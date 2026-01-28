using ModelingToolkit
using BatteryPeck
using Plots
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

p = Chen2020()

@mtkbuild sys = SingleCellPack()

exp = Experiment([
    PowerStep(20, 3600)
])

@time sol = simulate(sys, exp)

plot(sol[sys.t], sol[sys.source.v])
# plot(sol, vars=(sys.t, sys.cell.el.ηₑ), xlabel="Time (s)", ylabel="Electrolyte Potential Drop (V)", title="Electrolyte Potential Drop over Time")
# plot!(sol, vars=(sys.t, sys.cell.ηₑ2), xlabel="Time (s)", ylabel="Electrolyte Potential Drop (V)", title="Electrolyte Potential Drop over Time")


# plot(sol, vars=(sys.t, sys.V), xlabel="Time (s)", ylabel="Cell Voltage (V)", title="Single Cell Voltage Response")