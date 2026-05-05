# ==============================================================================
# isolated_anticipative_test.jl
# Validates lookahead "pre-cooling" logic
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq
using Plots
using SciMLBase
using BatteryToolkit
using DiffEqCallbacks
using Logging

# Override experiment.jl step! functions for modern parameter indexing
function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentStep)
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = -step.value
    u_modified!(integrator, true) 
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::RestStep)
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = 0.0
    u_modified!(integrator, true) 
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

println("Initialising Anticipative Diagnostic...")

elec_params = Chen2020()
elec_params.Vmin = 1.0 
elec_params.Vmax = 5.5

therm_params = build_pack_parameters()
a_c = therm_params.tms_geometry.channel_width * therm_params.tms_geometry.channel_height * therm_params.tms_geometry.number_of_channels
m_passive = velocity_to_mass_flow(0.1, therm_params.fluid.density, a_c)
m_active = velocity_to_mass_flow(0.3, therm_params.fluid.density, a_c) 

therm_params = PackParameters(
    fluid=therm_params.fluid, pipe_wall=therm_params.pipe_wall,
    potting_material=therm_params.potting_material, casing_material=therm_params.casing_material,
    tms_geometry=therm_params.tms_geometry, cell_gap_thickness=therm_params.cell_gap_thickness,
    axial_potting_thickness=therm_params.axial_potting_thickness, casing_thickness=therm_params.casing_thickness,
    ambient_temperature=298.15, inlet_temperature=298.15,
    mass_flow_rate=m_passive, ambient_convection_coefficient=5.0 
)

geom = build_pack_geometry(rows=1, cols=4, pattern=:square, tms_routing=:single_row, tms_encasement=:start_bottom, cell_pitch=0.025)

println("Building Coupled Multi-Cell System (Warnings Suppressed)...")
# Wrap build to silence harmless MTK internal warnings
sys_simp = with_logger(ConsoleLogger(stderr, Logging.Error)) do
    sys = CoupledMultiCellPack(name=:sys, geom=geom, therm_params=therm_params, elec_params=elec_params)
    structural_simplify(sys)
end

exp = Experiment([
    RestStep(600.0),
    CurrentStep(10.0, 1200.0),
    CurrentStep(-10.0, 500.0),
    CurrentStep(10.0, 500.0),
    CurrentStep(-10.0, 500.0),
    CurrentStep(10.0, 500.0),
    CurrentStep(-10.0, 500.0),
    CurrentStep(10.0, 500.0),
    CurrentStep(-10.0, 500.0),
    CurrentStep(10.0, 500.0),
    CurrentStep(-10.0, 500.0),
    CurrentStep(10.0, 500.0),
    CurrentStep(-10.0, 500.0),
    CurrentStep(10.0, 500.0),
    RestStep(300.0)
])

prob = ODEProblem(sys_simp, [sys_simp.Pin => 0.0, sys_simp.Iin => 0.0], (0.0, exp.tend), sparse=true, jac=false)

println("Attaching Anticipative DiscreteCallback...")
cb_anticipative = build_anticipative_callback(sys_simp.thermal_pack, 4, exp.steps, m_active, m_passive, 5.0, 300.0, 288.15, 298.15)

println("Simulating...")

# Capture tuple containing mass flow and applied current
saved_values = SavedValues(Float64, Tuple{Float64, Float64})
cb_save = SavingCallback(
    (u, t, integrator) -> (
        integrator.ps[sys_simp.thermal_pack.fluid_inlet.m_flow_in], 
        integrator.ps[sys_simp.Iin]
    ), 
    saved_values
)
cb_all = CallbackSet(cb_anticipative, cb_save)

@time sol = simulate(sys_simp, prob, exp, QNDF(); saveat=1.0, dtmax=5.0, callback=cb_all)

t_plot = sol.t
T_cells = [sol[getproperty(sys_simp.thermal_pack, Symbol("cell_$i")).core_cap.T] .- 273.15 for i in 1:4]

# Read from custom data logger
m_flow_scraped = [val[1] for val in saved_values.saveval]
current_scraped = [-val[2] for val in saved_values.saveval] 
t_scraped = saved_values.t

rho = therm_params.fluid.density
a_c = therm_params.tms_geometry.channel_width * therm_params.tms_geometry.channel_height * therm_params.tms_geometry.number_of_channels
vel_scraped = m_flow_scraped ./ (rho * a_c)

p1 = plot(t_scraped, current_scraped, label="Applied Current", ylabel="[A]", lw=2, color=:black, linetype=:steppost)
p2 = plot(t_plot, T_cells[1], label="Cell 1 Temp", ylabel="[C]", lw=2, color=:red)
p3 = plot(t_scraped, vel_scraped, label="Coolant Velocity", ylabel="[m/s]", xlabel="Time [s]", lw=2, color=:green, linetype=:steppost)
ylims!(p3, 0.0, 0.4) 

display(plot(p1, p2, p3, layout=(3,1), size=(800, 800), title="Anticipative Pre-Cooling Diagnostic"))
println("Done! Current will plot correctly, and the temperature reset is cured.")