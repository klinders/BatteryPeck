# =====================================================================================================================
# coupled_thermal_SPMe.jl
# Runs dynamic charge profile and plots core vs shell temperatures
# =====================================================================================================================

using ModelingToolkit
using Plots
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using Plots.Measures

using Revise
using BatteryPeck

Revise.revise()

# Load parameters
p = Chen2020()
p.n.c₀ = p.n.z_0 * p.n.c₊  
p.p.c₀ = p.p.z_0 * p.p.c₊  
p.Vmin = 2.4 

# Build coupled system
# h_conv = 15 (roughly natural convection, higher values represent increased cooling capacity)
@mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=15.0)

# High C-rate experiment to force a thermal gradient
exp = Experiment([
    RestStep(10),
    CurrentStep(-5 * 1, 3600/1.5),
    RestStep(600)
])

println("Starting Coupled Simulation...")
@time sol = simulate(sys, exp; saveat=1.0)

t = sol.t
v = sol[sys.cell.v]            
pack_current = sol[sys.I]

# Extract thermal variables and convert to Celsius
T_core = sol[sys.thermal.core_cap.T] .- 273.15
T_shell = sol[sys.thermal.shell_cap.T] .- 273.15

# Plotting
m = 7mm
xlim_range = (t[1], t[end]) 
matlab_colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]

# Temperature Plot
p_temp = plot(t, T_core, label="Core (Jellyroll) Temp", ylabel="Temperature (°C)", 
              title="Coupled Thermal-Electrochemical Response", lw=2, color=matlab_colors[2],
              legend=:topleft, xlims=xlim_range, margin=m)
plot!(p_temp, t, T_shell, label="Shell (Casing) Temp", lw=2, color=matlab_colors[1])

# Voltage Plot
p_volt = plot(t, v, label="Voltage", ylabel="Voltage (V)", xlabel="Time (s)", 
              lw=2, color=matlab_colors[4], legend=:bottomright, xlims=xlim_range, margin=m)

l = @layout [a; b]
p_combined = plot(p_temp, p_volt, layout=l, size=(1000, 800))
display(p_combined)