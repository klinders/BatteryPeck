module BatteryToolkit

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal

include("ParameterSets/Base.jl")
include("helpers.jl")
include("ParameterSets/Chen2020.jl")
include("CellModels/SPMe/SPMe.jl")
include("CellModels/SPMe/LithiumPlating.jl")

include("ThermalModels/CellThermal.jl")
include("ThermalModels/PackGeometry.jl")
include("ThermalModels/PackParameters.jl")
include("ThermalModels/TMSComponents.jl")
include("ThermalModels/SystemBuilder.jl")

include("experiment.jl")
include("solvers.jl")
include("PackModels/SingleCellPack.jl")
include("PackModels/MultiCellPack.jl")
include("PackModels/SingleCellCoreShellPack.jl")


export SPMe, Chen2020, BatteryParameters, PowerStep, CurrentStep, DriveStep,
       RestStep, ChargeStep, CurrentDriveStep, SingleCellPack, MultiCellPack, Experiment, 
       simulate, CoreShellCell, TemperatureDependentJellyroll, SingleCellCoreShellPack,
       PackGeometry, build_pack_geometry,
       FluidProperties, SolidProperties, TMSGeometry, PackParameters,
       get_coolant_properties, get_solid_properties,
       FluidPort, PipeWallNode, FluidNode, TMSNode,
       build_pack_system, build_pack_parameters
end