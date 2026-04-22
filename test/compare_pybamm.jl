# =====================================================================================================================
# compare_pybamm.jl
#
# Validates the Julia SPMe implementation directly against a PyBaMM baseline CSV
# for a 1C constant current discharge using the Chen2020 parameter set.
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using CSV
using DataFrames
using Plots
using Plots.Measures
using Logging
using Statistics
using Printf

using Revise
using BatteryToolkit

Revise.revise()

# Load PyBaMM baseline data
pybamm_file = joinpath(@__DIR__, "..", "data", "PyBaMM", "pybamm_spme_1c.csv")
df_pybamm = CSV.read(pybamm_file, DataFrame)
t_pb = df_pybamm[!, "Time [s]"]
v_pb = df_pybamm[!, "Voltage [V]"]

# Load parameter set
p = Chen2020()

# Enforce standard voltage limits
p.Vmin = 2.5
p.Vmax = 4.3

# Initialise at exactly 100% SoC to match PyBaMM default state
p.n.c₀ = p.n.z_100 * p.n.c₊  
p.p.c₀ = p.p.z_100 * p.p.c₊

# Build basic isothermal pack for 1:1 comparison
@mtkbuild sys = SingleCellPack(params=p, config=(1,1))

# Define 1C discharge (5A for Chen2020)
exp = Experiment([
    RestStep(1),
    CurrentStep(5.0, 4000.0) 
])

println("Running Julia SPMe 1C discharge...")

# Run simulation
sol = with_logger(NullLogger()) do
    simulate(sys, exp; saveat=1.0, verbose=false)
end

# Extract Julia results
t_jl = sol.t
v_jl = sol[sys.cell.v]

# Plot formatting
m = 7mm
colours = ["#0072BD", "#D95319", "#EDB120"]

p_compare = plot(title="SPMe 1C Discharge: Julia vs PyBaMM", 
                 xlabel="Time (s)", ylabel="Voltage (V)", 
                 legend=:topright, margin=m, size=(800, 600))

# Plot PyBaMM baseline
plot!(p_compare, t_pb, v_pb, label="PyBaMM SPMe", lw=3, color=:black, linestyle=:dash)

# Plot Julia model
plot!(p_compare, t_jl, v_jl, label="Julia SPMe", lw=2, color=colours[1])

# Calculate and display simple error metric if lengths allow
min_len = min(length(v_jl), length(v_pb))
if min_len > 10
    rmse = sqrt(mean((v_jl[1:min_len] .- v_pb[1:min_len]).^2))
    @printf("Comparison RMSE: %7.2f mV\n", rmse * 1000)
end

display(p_compare)