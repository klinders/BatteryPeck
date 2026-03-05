# =====================================================================================================================
# test_validation_chen.jl
# Compares SingleCellCoreShellPack simulation to experimental Chen 2020 dataset
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using CSV
using DataFrames
using Plots
using Plots.Measures

using Revise
using BatteryPeck

Revise.revise()

# Load experimental data
data_file = joinpath(@__DIR__, "..", "data", "Chen2020", "LGM50_cell03.csv")

# Skip 13 rows of metadata from Maccor tester and load into DataFrame to extract plotting arrays
df = CSV.read(data_file, DataFrame, skipto=15, header=14)

t_exp = df[!, "Test Time [s]"]
v_exp = df[!, "Voltage [V]"]
T_surf_exp = df[!, "Temperature Cell [degC]"]
T_chamber_exp = df[!, "Temperature Chamber [degC]"]
T_amb_K = T_chamber_exp[1] + 273.15

# Process the experimental current for plotting
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

exp_duration = t_exp[end]

# Load parameter set
p = Chen2020()

soc_init = 0.063

# Map the SoC to the exact electrode stoichiometries
z_n_init = p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)
z_p_init = p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)

# Set the initial solid concentrations
p.n.c₀ = z_n_init * p.n.c₊  
p.p.c₀ = z_p_init * p.p.c₊

p.Vmin = 2.4
p.Vmax = 4.3

# Build coupled core-shell system
# Guess natural convection coefficient (e.g., h=15 W/m²K)
@mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=18.0)

# Run experiment
exp = Experiment([
    CurrentDriveStep(data_file, exp_duration)
])

println("Running Validation Simulation...")
@time sol = simulate(sys, exp, saveat=10.0)

# Extract simulation results
t_sim = sol.t
v_sim = sol[sys.cell.v]
i_sim = sol[sys.I]
T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15

# Plotting
m = 7mm
xlim_range = (0, exp_duration)
colors = ["#0072BD", "#D95319", "#EDB120"]

# Voltage plot
p_volt = plot(t_exp, v_exp, label="Exp. voltage", ylabel="Voltage (V)", 
              lw=2, color=:black, linestyle=:dash, legend=:bottomright, margin=m, xlims=xlim_range)
plot!(p_volt, t_sim, v_sim, label="Sim. voltage", lw=2, color=colors[1])

# Temperature plot
p_temp = plot(t_exp, T_surf_exp, label="Exp. surface temp", ylabel="Temperature (°C)", 
              lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)
plot!(p_temp, t_sim, T_shell_sim, label="Sim. surface temp", lw=2, color=colors[2])

# Current plot
p_curr = plot(t_exp, i_exp, label="Exp. current", xlabel="Time (s)", ylabel="Current (A)", 
              lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)
plot!(p_curr, t_sim, i_sim, label="Sim. current", lw=2, color=colors[3])

# Combine layouts
l = @layout [a; b; c]
p_combined = plot(p_volt, p_temp, p_curr, layout=l, size=(1000, 1000), plot_title="Coupled core-shell thermal SPMe validation (Chen 2020 dataset)")
display(p_combined)