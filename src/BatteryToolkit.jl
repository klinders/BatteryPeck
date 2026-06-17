module BatteryToolkit

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal
using DiffEqCallbacks

include("ParameterSets/Base.jl")
include("helpers.jl")
include("ParameterSets/Chen2020.jl")
include("ParameterSets/OKane2022.jl")

include("CellModels/SPMe/SPMe.jl")
include("ThermalModels/CellThermal.jl")
include("ThermalModels/PackGeometry.jl")
include("ThermalModels/PackParameters.jl")
include("ThermalModels/TMSComponents.jl")
include("ThermalModels/SystemBuilder.jl")
include("ThermalModels/TMSControl.jl")

include("experiment.jl")
include("solvers.jl")
include("PackModels/SingleCellPack.jl")
include("PackModels/MultiCellPack.jl")
include("PackModels/SingleCellCoreShellPack.jl")
include("PackModels/ExplicitPackSimulator.jl")

export SPMe, Chen2020, OKane2022, BatteryParameters, PowerStep, CurrentStep, TargetCurrentStep, DriveStep,
       RestStep, ChargeStep, CurrentDriveStep, SingleCellPack, MultiCellPack, Experiment, 
       simulate, CoreShellCell, TemperatureDependentJellyroll, SingleCellCoreShellPack,
       PackGeometry, build_pack_geometry,
       FluidProperties, SolidProperties, TMSGeometry, PackParameters,
       get_coolant_properties, get_solid_properties,
       FluidPort, PipeWallNode, FluidNode, TMSNode,
       build_pack_system, build_pack_parameters,
       velocity_to_mass_flow, 
       ExplicitPackSimulator, build_pack_simulator, simulate_pack!,
       TMSStrategy, ReactiveTMS, AnticipativeTMS, evaluate_tms_state, 
       get_future_load_avg, compile_experiment
end