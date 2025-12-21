# BatteryToolkit

a ModelingToolkit based battery pack inspired by PyBaMM

## How to run

Running a simulation can be as easy as

```julia
using ModelingToolkit
using BatteryPeck

@mtkbuild sys = SingleCellPack(config=(96,2), Qcell=78);

e = Experiment([
    PowerStep(11000, 500),
    RestStep(8*3600),
]*7)

@time sol = simulate(sys,e, saveat=1)
```

## Experiment

The experiment can be defined using steps like
```julia
Experiment([
    PowerStep(11000, 500),
    RestStep(8*3600),
    DriveStep("speed comparison/power.csv", 3600/2),
    RestStep(8*3600),
    DriveStep("speed comparison/power.csv", 3600/2),
    ChargeStep(0.8, 11000),
    RestStep(6*3600)
]*100)
```

## Plotting solutions
Solutions can be accessed like you're used to!

```julia
using Plots

plot(sol, idxs=[sys.cell.soc])
```
