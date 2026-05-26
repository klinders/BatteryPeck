module BatteryToolkit

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("ParameterSets/Base.jl")
include("helpers.jl")
include("ParameterSets/Chen2020.jl")
include("ParameterSets/OKane2022.jl")
include("CellModels/SPMe/SPMe.jl")
include("experiment.jl")
include("solvers.jl")
include("PackModels/SingleCellPack.jl")
include("PackModels/MultiCellPack.jl")

export SPMe, Chen2020, OKane2022, BatteryParameters, PowerStep,CurrentStep, DriveStep,RestStep,ChargeStep, SingleCellPack, MultiCellPack, Experiment, simulate

end  # module BatteryPeck