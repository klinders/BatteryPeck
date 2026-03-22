# =====================================================================================================================
# validation_chen.jl
# Compare SingleCellCoreShellPack simulation to experimental Chen dataset across 4 segments
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using CSV
using DataFrames
using Plots
using Plots.Measures
using DataInterpolations
using Logging
using Statistics
using Printf

using Revise
using BatteryToolkit

Revise.revise()

# Load experimental data
data_file = joinpath(@__DIR__, "..", "data", "Chen2020", "LGM50_cell03.csv")

# Skip metadata and load dataframe to extract plotting arrays
df = CSV.read(data_file, DataFrame, skipto=15, header=14)

t_exp = df[!, "Test Time [s]"]
v_exp = df[!, "Voltage [V]"]
T_surf_exp = df[!, "Temperature Cell [degC]"]
T_chamber_exp = df[!, "Temperature Chamber [degC]"]
T_amb_K = T_chamber_exp[1] + 273.15

# Read voltage SOC lookup table
lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
lut = CSV.read(lut_file, DataFrame)

# Create interpolator mapping voltage to SOC
v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)

# Process experimental current for plotting
raw_i_exp = df[!, "Current [A]"]
md_exp = df[!, "Md"]
i_exp = zeros(Float64, length(raw_i_exp))

for j in eachindex(raw_i_exp)
    mode_str = strip(String(md_exp[j]))
    if mode_str == "C"
        i_exp[j] = -raw_i_exp[j]
    elseif mode_str == "D"
        i_exp[j] = raw_i_exp[j]
    else
        i_exp[j] = 0.0
    end
end

# Identify segments
starts = Int[]
ends = Int[]

# Find start indices where mode changes from rest to discharge starting at final rest step
for i in 2:nrow(df)
    if strip(String(md_exp[i-1])) == "R" && strip(String(md_exp[i])) == "D"
        push!(starts, i - 1)
    end
end

# Find corresponding end indices where mode changes from rest to charge
for s in starts
    e = s
    while e < nrow(df)
        if strip(String(md_exp[e])) == "R" && strip(String(md_exp[e+1])) == "C"
            push!(ends, e)
            break
        end
        e += 1
    end
    # Safety catch if file ends without charging
    if e == nrow(df) && length(ends) < length(starts)
        push!(ends, e)
    end
end

# Discard first segment
starts = starts[2:end]
ends = ends[2:end]

# Mask experimental arrays with NaNs to create visual gaps
v_exp_masked = fill(NaN, length(v_exp))
T_surf_exp_masked = fill(NaN, length(T_surf_exp))
T_chamber_exp_masked = fill(NaN, length(T_chamber_exp))
i_exp_masked = fill(NaN, length(i_exp))

for (s, e) in zip(starts, ends)
    v_exp_masked[s:e] .= v_exp[s:e]
    T_surf_exp_masked[s:e] .= T_surf_exp[s:e]
    T_chamber_exp_masked[s:e] .= T_chamber_exp[s:e]
    i_exp_masked[s:e] .= i_exp[s:e]
end

# Initialise simulation loop
sim_results = []

println("Running validation simulations across $(length(starts)) segments...\n")

for (idx, (s, e)) in enumerate(zip(starts, ends))
    # Extract segment data
    t_seg = t_exp[s:e]
    v_seg = v_exp[s:e]
    T_surf_seg = T_surf_exp[s:e]
    i_seg = i_exp[s:e]
    T_chamber_seg = T_chamber_exp[s:e] .+ 273.15
    
    # Local time vector starting from zero for simulator
    t_sim_input = t_seg .- t_seg[1]
    
    # Initialise SOC based on segment starting voltage
    v_init = v_seg[1]
    soc_init = v_to_soc_interp(v_init) 
    
    # Load parameter set
    p = Chen2020()

    # Map SOC to stoichiometries
    z_n_init = p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)
    z_p_init = p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)

    # Set initial solid concentrations
    p.n.c₀ = z_n_init * p.n.c₊  
    p.p.c₀ = z_p_init * p.p.c₊

    p.Vmin = 1.0
    p.Vmax = 5.0

    # Build coupled core shell system
    @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=19.5, T_ambient=T_amb_K)

    # Format data for current drive step
    raw_dt = diff(t_sim_input)
    valid_indices = findall(x -> x > 0.0, raw_dt)
    clean_dt = raw_dt[valid_indices]
    clean_i = i_seg[valid_indices]
    clean_T_amb = T_chamber_seg[valid_indices .+ 1]
    
    # Define experiment for specific segment using time and current arrays
    exp_seg = Experiment([
        CurrentDriveStep(Any[clean_dt, clean_i, clean_T_amb], t_sim_input[end])
    ])

    # Run simulation
    local sol = with_logger(NullLogger()) do
        simulate(sys, exp_seg; saveat=10.0, verbose=false)
    end
    
    # Extract simulation results and shift time back to global timeframe
    t_sim_seg = sol.t .+ t_seg[1]
    v_sim_seg = sol[sys.cell.v]
    i_sim_seg = sol[sys.I]
    T_shell_sim_seg = sol[sys.thermal.shell_cap.T] .- 273.15
    
    # Interpolate for metrics calculation against experimental data points
    v_sim_interp = LinearInterpolation(v_sim_seg, t_sim_seg)
    T_shell_sim_interp = LinearInterpolation(T_shell_sim_seg, t_sim_seg)

    t_seg_clamped = clamp.(t_seg, t_sim_seg[1], t_sim_seg[end])
    v_sim_mapped = v_sim_interp.(t_seg_clamped)
    T_sim_mapped = T_shell_sim_interp.(t_seg_clamped)

    # Calculate RMSE
    rmse_v = sqrt(mean((v_sim_mapped .- v_seg).^2))
    rmse_T = sqrt(mean((T_sim_mapped .- T_surf_seg).^2))

    # Calculate R squared
    ss_res_v = sum((v_seg .- v_sim_mapped).^2)
    ss_tot_v = sum((v_seg .- mean(v_seg)).^2)
    r2_v = ss_tot_v == 0 ? 0.0 : 1.0 - (ss_res_v / ss_tot_v)

    ss_res_T = sum((T_surf_seg .- T_sim_mapped).^2)
    ss_tot_T = sum((T_surf_seg .- mean(T_surf_seg)).^2)
    r2_T = ss_tot_T == 0 ? 0.0 : 1.0 - (ss_res_T / ss_tot_T)

    # Print segment metrics
    println("Segment $idx (starts at t = $(round(t_seg[1], digits=1))s):")
    @printf("  Starting voltage | %.4f V\n", v_init)
    @printf("  Voltage          | RMSE: %7.2f mV | R²: %.4f\n", rmse_v * 1000, r2_v)
    @printf("  Temperature      | RMSE: %7.2f °C | R²: %.4f\n\n", rmse_T, r2_T)
    
    # Save for plotting overlay
    push!(sim_results, (t_sim_seg, v_sim_seg, T_shell_sim_seg, i_sim_seg))
end

# Plot layout
m = 7mm
xlim_range = (0, t_exp[end])
colours = ["#0072BD", "#D95319", "#EDB120"]

# Plot voltage
p_volt = plot(t_exp, v_exp_masked, label="Exp. voltage", ylabel="Voltage (V)", 
              lw=2, color=:black, linestyle=:dash, legend=:bottomleft, margin=m, xlims=xlim_range)
              
# Highlight segment start points
scatter!(p_volt, t_exp[starts], v_exp[starts], marker=:cross, color=:black, markersize=6, label="Experiment start")

# Overlay simulation results
for (i, res) in enumerate(sim_results)
    t_s, v_s, _, _ = res
    lbl = i == 1 ? "Sim. voltage" : ""
    plot!(p_volt, t_s, v_s, label=lbl, lw=2, color=colours[1])
end

# Plot temperature
p_temp = plot(t_exp, T_surf_exp_masked, label="Exp. surface temp", ylabel="Temperature (°C)", 
              lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)
plot!(p_temp, t_exp, T_chamber_exp_masked, label="Exp. chamber temp", lw=2, color=:gray, linestyle=:dot)

# Highlight segment start points
scatter!(p_temp, t_exp[starts], T_surf_exp[starts], marker=:cross, color=:black, markersize=6, label="")

# Overlay simulation results
for (i, res) in enumerate(sim_results)
    t_s, _, T_s, _ = res
    lbl = i == 1 ? "Sim. surface temp" : ""
    plot!(p_temp, t_s, T_s, label=lbl, lw=2, color=colours[2])
end

# Plot current
p_curr = plot(t_exp, i_exp_masked, label="Exp. current", xlabel="Time (s)", ylabel="Current (A)", 
              lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)

# Highlight segment start points
scatter!(p_curr, t_exp[starts], i_exp[starts], marker=:cross, color=:black, markersize=6, label="")

# Overlay simulation results
for (i, res) in enumerate(sim_results)
    t_s, _, _, i_s = res
    lbl = i == 1 ? "Sim. current" : ""
    plot!(p_curr, t_s, i_s, label=lbl, lw=2, color=colours[3])
end

# Combine layouts
l = @layout [a; b; c]
p_combined = plot(p_volt, p_temp, p_curr, layout=l, size=(1000, 1000), plot_title="Coupled core shell thermal SPMe validation")
display(p_combined)