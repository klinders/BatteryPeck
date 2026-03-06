# =====================================================================================================================
# heat_generation_SPMe.jl
#
# Runs a (dis)charge profile and plots individual heat generation terms
# =====================================================================================================================

# Import packages
using ModelingToolkit
using Plots
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using Plots.Measures

# Import module
using Revise
using BatteryPeck

Revise.revise()

# Load parameters
p = Chen2020()

# 0% SoC
p.n.c₀ = p.n.z_0 * p.n.c₊  # Empty negative electrode (Graphite)
p.p.c₀ = p.p.z_0 * p.p.c₊  # Fill positive electrode (NMC811)
p.Vmin = 2.4 # Lower event limit to prevent instant solver termination

# Build system
@mtkbuild sys = SingleCellPack(params=p, config=(1,1))

# Define experiment
exp = Experiment([
    RestStep(1),
    #PowerStep(30, 3600/1.5),

    CurrentStep(-5 * 0.2, 3600/0.15), # 0.2C charge
    #CurrentStep(-5 * 0.5, 3600/0.45), # 0.5C charge
    #CurrentStep(-5 * 1.0, 3600/1.0), # 1.0C charge
    #CurrentStep(-5 * 1.5, 3600/1.5), # 1.5C charge
    #CurrentStep(-5 * 2.0, 3600/2.0), # 2.0C charge

    #CurrentStep(5 * 0.2, 3600/0.15), # 0.2C discharge
    #CurrentStep(5 * 0.5, 3600/0.45), # 0.5C discharge
    #CurrentStep(5 * 1.0, 3600/1.0), # 1.0C discharge
    #CurrentStep(5 * 1.5, 3600/1.5), # 1.5C discharge
    #CurrentStep(5 * 2.0, 3600/2.0), # 2.0C discharge
    RestStep(60),
])

println("Starting Simulation...")
@time sol = simulate(sys, exp; saveat=1.0)

# Extract time and check end time
t = sol.t
println("Simulation ended at t = $(t[end]) seconds")
println("Solver Status: ", sol.retcode)

if t[end] < 10.0 # --- Debug for instant abort ---
    println("ERROR: Simulation aborted almost immediately. Check initial conditions or Vmin/Vmax events.")
else
    # Extract variables from solution
    v = sol[sys.cell.v]            
    pack_current = sol[sys.I]
    pack_power = sol[sys.P]

    q_rev = sol[sys.cell.Q_rev]
    q_irr = sol[sys.cell.Qᵢ]
    q_ohm_e = sol[sys.cell.Qₑ]
    q_ohm_s = sol[sys.cell.Qₛ]
    q_film = sol[sys.cell.Qf]
    q_total = sol[sys.cell.Q_total]

    # ---------------------------------------------------------
    # Bar chart data processing
    # ---------------------------------------------------------
    N = length(t)
    # Create 6 bin edges to make 5 equal sections across array length
    bin_edges = round.(Int, range(1, stop=N, length=6))
    
    labels_bar = String[]
    q_rev_bar = zeros(5); q_irr_bar = zeros(5); q_e_bar = zeros(5)
    q_s_bar = zeros(5); q_f_bar = zeros(5)
    
    for i in 1:5
        idx_start = bin_edges[i]
        idx_end = bin_edges[i+1]
        
        # Sum absolute values of heat
        sum_rev = sum(abs.(q_rev[idx_start:idx_end]))
        sum_irr = sum(abs.(q_irr[idx_start:idx_end]))
        sum_e   = sum(abs.(q_ohm_e[idx_start:idx_end]))
        sum_s   = sum(abs.(q_ohm_s[idx_start:idx_end]))
        sum_f   = sum(abs.(q_film[idx_start:idx_end]))
        
        # Calculate total absolute heat, and prevent division by zero during rest steps
        total_abs = sum_rev + sum_irr + sum_e + sum_s + sum_f
        total_abs = total_abs > 0 ? total_abs : 1.0 
        
        # Convert to percentages
        q_rev_bar[i] = 100 * sum_rev / total_abs
        q_irr_bar[i] = 100 * sum_irr / total_abs
        q_e_bar[i]   = 100 * sum_e / total_abs
        q_s_bar[i]   = 100 * sum_s / total_abs
        q_f_bar[i]   = 100 * sum_f / total_abs
        
        # Create x-axis time period labels
        t_start = round(Int, t[idx_start])
        t_end = round(Int, t[idx_end])
        push!(labels_bar, "$(t_start)s-\n$(t_end)s")
    end

    # ---------------------------------------------------------
    # Plotting code
    # ---------------------------------------------------------
    
    # Time
    xlim_range = (t[1], t[end]) 

    # Margins
    m = 7mm

    # MATLAB colour scheme
    matlab_colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]

    # Heat sources vs. time (top)
    p_heat = plot(t, q_rev, label="Reversible", ylabel="Volumetric heat (W/m³)", 
                  title="Heat generation sources", 
                  lw=2, legend=:bottomright, xlims=xlim_range, 
                  left_margin=m, right_margin=m, top_margin=m, bottom_margin=0mm, 
                  color=matlab_colors[1])
    plot!(p_heat, t, q_irr, label="Irreversible", lw=2, color=matlab_colors[2])
    plot!(p_heat, t, q_ohm_e, label="Electrolyte", lw=2, color=matlab_colors[3])
    plot!(p_heat, t, q_ohm_s, label="Solid phase", lw=2, color=matlab_colors[4])
    plot!(p_heat, t, q_film, label="Film", lw=2, color=matlab_colors[5])
    plot!(p_heat, t, q_total, label="Total", lw=2, linestyle=:dash, color=:black)

    # ---------------------------------------------------------
    # Voltage and current vs. time (middle)
    # ---------------------------------------------------------
    p_volt_curr = plot(t, v, label="Voltage", ylabel="Voltage (V)", lw=2, color=matlab_colors[1], 
                       legend=:bottomright, xlims=xlim_range,
                       left_margin=m, right_margin=m, top_margin=0mm, bottom_margin=0mm)
    # ylims=(2.5, 4.6)
    
    # Dummy trace for current so both labels share primary legend box
    plot!(p_volt_curr, [NaN], [NaN], label="Current", lw=2, color=matlab_colors[2])

    # Secondary axis (current)
    p_twin = twinx(p_volt_curr)
    
    # Add identical invisible legend to twin axis to squeeze it by exact same amount as plot beneath it
    plot!(p_twin, t, pack_current, label="", ylabel="Current (A)", lw=2, color=matlab_colors[2], 
          xlims=xlim_range, 
          left_margin=m, right_margin=m, top_margin=0mm, bottom_margin=0mm,
          legend=:bottomright, 
          foreground_color_legend=:transparent, 
          background_color_legend=:transparent, 
          legend_font_color=:transparent)
          
    # Feed identical strings into invisible legend so bounding boxes match perfectly
    plot!(p_twin, [NaN], [NaN], label="Voltage", color=:transparent)
    plot!(p_twin, [NaN], [NaN], label="Current", color=:transparent)

    # ---------------------------------------------------------
    # Power vs. time (bottom)
    # ---------------------------------------------------------
    p_pow = plot(t, pack_power, label="Power", xlabel="Time (s)", ylabel="Power (W)", lw=2, 
                 color=matlab_colors[4], legend=:bottomright, xlims=xlim_range, 
                 left_margin=m, right_margin=m, top_margin=0mm, bottom_margin=m)

    # ---------------------------------------------------------
    # Stacked bar chart (right column)
    # ---------------------------------------------------------
    stack_5 = q_rev_bar .+ q_irr_bar .+ q_e_bar .+ q_s_bar .+ q_f_bar 
    stack_4 = q_rev_bar .+ q_irr_bar .+ q_e_bar .+ q_s_bar
    stack_3 = q_rev_bar .+ q_irr_bar .+ q_e_bar
    stack_2 = q_rev_bar .+ q_irr_bar
    stack_1 = q_rev_bar

    p_bar = plot(title="Absolute contribution heat sources", ylabel="Percentage (%)", 
                 xlabel="Time period", legend=:outertopright, margin=m)
    
    bar!(p_bar, labels_bar, stack_5, label="Film", color=matlab_colors[5], lw=0)
    bar!(p_bar, labels_bar, stack_4, label="Solid phase", color=matlab_colors[4], lw=0)
    bar!(p_bar, labels_bar, stack_3, label="Electrolyte", color=matlab_colors[3], lw=0)
    bar!(p_bar, labels_bar, stack_2, label="Irreversible", color=matlab_colors[2], lw=0)
    bar!(p_bar, labels_bar, stack_1, label="Reversible", color=matlab_colors[1], lw=0)

    # Multi-figure layout
    l = @layout [
        [a{0.6h}; b; c] d{0.35w}
    ]
    
    p_combined = plot(p_heat, p_volt_curr, p_pow, p_bar, layout=l, size=(1400, 1000))
    
    display(p_combined)
end


# include("BatteryPeck/test/test_heat_cc.jl")