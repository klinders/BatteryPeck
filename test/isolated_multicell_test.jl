# ==============================================================================
# isolated_multicell_test.jl
# Diagnoses the Virtual Series & Thermal Pack without Callback Interference
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq
using Plots
using SciMLBase
using BatteryToolkit

# experiment.jl overrides
function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentStep)
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = -step.value
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::RestStep)
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = 0.0
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

println("Initialising Multi-Cell Diagnostic...")

# Load parameters
elec_params = Chen2020()
elec_params.Vmin = 1.0 
elec_params.Vmax = 5.5

therm_params = build_pack_parameters()
geom = build_pack_geometry(rows=1, cols=4, pattern=:square, tms_routing=:single_row, tms_encasement=:start_bottom, cell_pitch=0.025)

println("Building Coupled Multi-Cell System...")
sys = CoupledMultiCellPack(name=:sys, geom=geom, therm_params=therm_params, elec_params=elec_params)
sys_simp = structural_simplify(sys)

# Define experiment: 60s Rest -> 3600s 1C Discharge (5A) -> 60s Rest
exp = Experiment([
    RestStep(60.0),
    CurrentStep(5.0, 3600.0), 
    RestStep(60.0)
])

# Modern MTK Problem Initialization mapping parameters
prob = ODEProblem(sys_simp, [sys_simp.Pin => 0.0, sys_simp.Iin => 0.0], (0.0, exp.tend), sparse=true, jac=false)

println("Simulating...")
@time sol = simulate(sys_simp, prob, exp, QNDF(); saveat=1.0)

t_plot = sol.t
v_plot = sol[sys_simp.V]
i_plot = sol[sys_simp.I]

# Extract all 4 cell temperatures
T_cells = [sol[getproperty(sys_simp.thermal_pack, Symbol("cell_$i")).core_cap.T] .- 273.15 for i in 1:4]

p1 = plot(t_plot, v_plot, label="Pack Voltage", ylabel="Voltage [V]", lw=2, color=:blue)
p2 = plot(t_plot, i_plot, label="Pack Current", ylabel="Current [A]", lw=2, color=:red)

p3 = plot(t_plot, T_cells[1], label="Cell 1", ylabel="Temperature [C]", xlabel="Time [s]", lw=2)
plot!(p3, t_plot, T_cells[2], label="Cell 2", lw=2)
plot!(p3, t_plot, T_cells[3], label="Cell 3", lw=2)
plot!(p3, t_plot, T_cells[4], label="Cell 4", lw=2)

display(plot(p1, p2, p3, layout=(3,1), size=(900, 900), title="Isolated 4-Cell Virtual Series Check"))
println("Done! Check if current holds steady and temperatures rise natively.")