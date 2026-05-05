# ==============================================================================
# isolated_test.jl
# Isolates the "Ghost Circuit" and "Comb Effect" on a single SPMe cell.
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots
using SciMLBase
using BatteryToolkit


# Define Pin and Iin explicitly as parameters
function SingleCellPack_Patched(;name, params=Chen2020(), Qcell=5)
    @parameters t
    
    @parameters begin
        Pin = 0.0
        Iin = 0.0
        Tin = 298.15
    end
    
    @variables begin
        V(t)
        I(t)
    end
    
    @named cell = SPMe(params=params, Q=Qcell)
    @named power = RealInput(guess=0)
    @named current = RealInput(guess=0)
    @named temp = RealInput(guess=298)
    @named source = Current()
    @named ground = Ground()
    
    eqs = [
        V ~ cell.v
        I ~ cell.i
        power.u ~ Pin
        current.u ~ Iin
        temp.u ~ Tin
        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        connect(temp, cell.T)
        source.I.u ~ power.u/cell.v + current.u
    ]
    return ODESystem(eqs, t; systems=[cell, power, temp, source, ground], name=name)
end

# Override the experiment.jl step! functions using integrator.ps
function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentStep)
    # Modern SciML parameter indexing
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = -step.value
    
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::RestStep)
    # Modern SciML parameter indexing
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = 0.0
    
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

println("Building patched SPMe system...")
sys = SingleCellPack_Patched(name=:sys)
sys_simp = structural_simplify(sys)

exp = Experiment([
    RestStep(60.0),
    CurrentStep(5.0, 3600.0), # 1C discharge for 1 hour
    RestStep(60.0)
])

println("Compiling ODE problem...")
# Modern MTK initialization: combine states and parameters into one map
prob = ODEProblem(sys_simp, [sys_simp.Pin => 0.0, sys_simp.Iin => 0.0], (0.0, exp.tend))

println("Simulating...")
@time sol = simulate(sys_simp, prob, exp, QNDF(); saveat=1.0)

t_plot = sol.t
v_plot = sol[sys_simp.V]
i_plot = sol[sys_simp.I]

p1 = plot(t_plot, v_plot, label="Cell Voltage", ylabel="Voltage [V]", lw=2, color=:blue)
p2 = plot(t_plot, i_plot, label="Cell Current", ylabel="Current [A]", xlabel="Time [s]", lw=2, color=:red)

display(plot(p1, p2, layout=(2,1), size=(800, 600), title="Isolated Current & Voltage Test"))
println("Done! Check the plot to ensure the current holds steady at 5A without resetting.")