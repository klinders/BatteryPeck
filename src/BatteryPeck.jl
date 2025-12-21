module BatteryPeck

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("ParameterSets/Base.jl")
include("helpers.jl")
include("ParameterSets/Chen2020.jl")
include("SPMe/SPMe.jl")
include("experiment.jl")
include("solvers.jl")
include("PackModels/SingleCellPack.jl")

export SPMe, Chen2020, BatteryParameters, PowerStep, DriveStep,RestStep,ChargeStep, SingleCellPack, Experiment, simulate

end  # module BatteryPeck