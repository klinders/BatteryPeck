module BatteryPeck

include("cell_model.jl")
include("pack_model.jl")
include("solvers.jl")

using .CellModel
using .PackModel
using .Solvers

ECM_Model = Model(battery_pack_ecm!,PackParameters(LG50T, 3), Model[])

end  # module BatteryPeck