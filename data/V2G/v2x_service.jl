using ModelingToolkit
using JSON
using GLMakie
using OrdinaryDiffEq
using Dates

using Revise
using BatteryToolkit

include("../../../src/services.jl")

@mtkbuild sys = SingleCellPack(config=(96,2), Qcell=78, params=OKane2022());

service = FCR

experiment = Experiment([
    DriveStep(joinpath(@__DIR__,"driving_power_wltp.csv"), 900.0),
    RestStep(8*3600),
    DriveStep(joinpath(@__DIR__,"driving_power_wltp.csv"), 900.0),
    service(16*3600-1800,dt=300)
]*90, DateTime(2020, 1, 1,9)) # Simulate from Jan 1, 2020 at 9 AM

# + [
#     RestStep(3*3600), 
#     DriveStep("data/power.csv", 1800.0),]
#     RestStep(3*3600),
#     DriveStep("data/power.csv", 1800.0),
#     service(17*3600)
# ]*2)

sol = simulate(sys,experiment, saveat=60);

# using JLD2
using DataFrames

# JLD2.@save "v0x_output.jld2" sol

df = DataFrame()

df.t = sol[sys.t]
df.v = sol[sys.cell.v]
df.i = sol[sys.cell.i]
df.q_plating = sol[sys.cell.plating.Q_loss]
df.q_sei = sol[sys.cell.sei.Q_loss]
df.q_loss = sol[sys.cell.Q_loss]

using CSV

CSV.write(joinpath(@__DIR__,"fcr_output.csv"), df)

# f2 = Figure()
# ax2 = Axis(f2[1,1], title="BatteryToolkit Plating and SEI", xlabel="Time [days]", ylabel="Concentration [mol/m³]")
# lines!(ax2, sol[sys.t]/3600/24, sol[sys.cell.plating.c_plating_x], label="Plated Li")
# lines!(ax2, sol[sys.t]/3600/24, sol[sys.cell.plating.c_dead_x], label="Dead Li")
# lines!(ax2, sol[sys.t]/3600/24, sol[sys.cell.sei.c_sei_x], label="SEI")
# axislegend()
# save(joinpath(@__DIR__,"btk_fcr_experiment_plating_sei.png"), f2)