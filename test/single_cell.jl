using ModelingToolkit
using BatteryPeck
using Plots

@mtkbuild sys = SingleCellPack(config=(96,2), Qcell=78);

e = Experiment([
    PowerStep(11000, 500),
    RestStep(8*3600),
]*7)

@time sol = simulate(sys,e, saveat=1);

plot(sol, vars=(sys.t, sys.V), xlabel="Time (s)", ylabel="Cell Voltage (V)", title="Single Cell Voltage Response")