module BatteryPeck

using DifferentialEquations

include("cell_model.jl")
include("pack_model.jl")
include("solvers.jl")

using .CellModel
using .PackModel
using .Solvers

export solve_ecm_battery

end  # module BatteryPeck