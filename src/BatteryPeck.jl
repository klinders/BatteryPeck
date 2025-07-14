module BatteryPeck

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

#include("cell_model.jl")
#include("pack_model.jl")
include("Models/DynamicModels/SimpleECM.jl")
include("Models/DegradationModels/WangEtAl2014.jl")
#include("solvers.jl")

# using .CellModel
# using .PackModel
# using .Solvers

# ECM_Model = Model(battery_pack_ecm!,PackParameters(LG50T, 3), Model[])

export SimpleECM, WangEtAl2014

end  # module BatteryPeck