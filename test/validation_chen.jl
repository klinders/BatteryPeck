# =====================================================================================================================
# validation_chen.jl
# Compare SingleCellCoreShellPack simulation to experimental Chen dataset
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

# Read voltage state of charge look up table
lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
lut = CSV.read(lut_file, DataFrame)

# Create interpolator mapping voltage to state of charge
v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)

# Calculate initial state of charge dynamically based on initial voltage
initial_measured_voltage = v_exp[1]
soc_init = v_to_soc_interp(initial_measured_voltage)

println("Dynamically initialised at SoC = $(round(soc_init*100, digits=2))% based on starting voltage of $(initial_measured_voltage)V")

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

exp_duration = t_exp[end]

# Load parameter set
p = Chen2020()

# Map state of charge to stoichiometries
z_n_init = p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)
z_p_init = p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)

# Set initial solid concentrations
p.n.c₀ = z_n_init * p.n.c₊  
p.p.c₀ = z_p_init * p.p.c₊

p.Vmin = 2.4
p.Vmax = 4.3

# Build coupled core shell system
@mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=19.5, T_ambient=T_amb_K)

# Define experiment
exp = Experiment([
    CurrentDriveStep(data_file, exp_duration)
])

println("Running validation simulation...")

# Run simulation while silencing internal package termination warnings
@time sol = with_logger(NullLogger()) do
    simulate(sys, exp; saveat=10.0, verbose=false)
end

# Extract simulation results
t_sim = sol.t
v_sim = sol[sys.cell.v]
i_sim = sol[sys.I]
T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15

# Interpolate simulation results to match experimental timestamps
v_sim_interp = LinearInterpolation(v_sim, t_sim)
T_shell_sim_interp = LinearInterpolation(T_shell_sim, t_sim)

t_exp_clamped = clamp.(t_exp, t_sim[1], t_sim[end])
v_sim_mapped = v_sim_interp.(t_exp_clamped)
T_sim_mapped = T_shell_sim_interp.(t_exp_clamped)

# Calculate root mean square error
rmse_v = sqrt(mean((v_sim_mapped .- v_exp).^2))
rmse_T = sqrt(mean((T_sim_mapped .- T_surf_exp).^2))

# Calculate R squared
ss_res_v = sum((v_exp .- v_sim_mapped).^2)
ss_tot_v = sum((v_exp .- mean(v_exp)).^2)
r2_v = 1.0 - (ss_res_v / ss_tot_v)

ss_res_T = sum((T_surf_exp .- T_sim_mapped).^2)
ss_tot_T = sum((T_surf_exp .- mean(T_surf_exp)).^2)
r2_T = 1.0 - (ss_res_T / ss_tot_T)

# Print validation metrics
println("\nValidation metrics:")
@printf("Voltage     | RMSE: %7.2f mV | R²: %.4f\n", rmse_v * 1000, r2_v)
@printf("Temperature | RMSE: %7.2f °C | R²: %.4f\n", rmse_T, r2_T)

# Plot results
m = 7mm
xlim_range = (0, exp_duration)
colors = ["#0072BD", "#D95319", "#EDB120"]

# Plot voltage
p_volt = plot(t_exp, v_exp, label="Exp. voltage", ylabel="Voltage (V)", 
              lw=2, color=:black, linestyle=:dash, legend=:bottomright, margin=m, xlims=xlim_range)
plot!(p_volt, t_sim, v_sim, label="Sim. voltage", lw=2, color=colors[1])

# Plot temperature
p_temp = plot(t_exp, T_surf_exp, label="Exp. surface temp", ylabel="Temperature (°C)", 
              lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)
plot!(p_temp, t_exp, T_chamber_exp, label="Exp. chamber temp", lw=2, color=:gray, linestyle=:dot)
plot!(p_temp, t_sim, T_shell_sim, label="Sim. surface temp", lw=2, color=colors[2])

# Plot current
p_curr = plot(t_exp, i_exp, label="Exp. current", xlabel="Time (s)", ylabel="Current (A)", 
              lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)
plot!(p_curr, t_sim, i_sim, label="Sim. current", lw=2, color=colors[3])

# Combine layouts
l = @layout [a; b; c]
p_combined = plot(p_volt, p_temp, p_curr, layout=l, size=(1000, 1000), plot_title="Coupled core shell thermal SPMe validation")
display(p_combined)