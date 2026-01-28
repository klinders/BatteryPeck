module BatteryPeck

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("ParameterSets/Base.jl")
include("helpers.jl")
include("ParameterSets/Chen2020.jl")
include("CellModels/SPMe/SPMe.jl")
include("experiment.jl")
include("solvers.jl")
include("PackModels/SingleCellPack.jl")

export SPMe, Chen2020, BatteryParameters, PowerStep,CurrentStep, DriveStep,RestStep,ChargeStep, SingleCellPack, Experiment, simulate

end  # module BatteryPeck