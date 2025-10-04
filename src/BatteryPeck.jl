module BatteryPeck

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SPMe/SPMe.jl")
include("ParameterSets/Base.jl")
include("ParameterSets/Chen2020.jl")
#include("solvers.jl")

export SPMe, Chen2020

end  # module BatteryPeck