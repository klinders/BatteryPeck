# =====================================================================================================================
# solvers.jl
#
# Main simulation function
# Inputs:
# - sys = model to be simulated, e.g. "SingleCellPack()"
# - experiment = struct defining the simulation protocol
# - alg = ODE algorithm (defaults to QNDF)
# - kwargs = keyword arguments for solver options of integrator
# =====================================================================================================================

using ModelingToolkit, OrdinaryDiffEq

"""
    simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, alg=QNDF(); reltol=1e-6, abstol=1e-7, callback=nothing, kwargs...)

Compile ODE problem every time for ease of use.
"""
function simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, alg=QNDF(); reltol=1e-5, abstol=1e-6, callback=nothing, kwargs...)
    
    # Proper ODEProblem syntax to avoid deprecation warnings.
    # Argument 2: Initial conditions (empty array forces reliance on kwargs `guesses`)
    # Argument 3: Timespan
    # Argument 4: Parameter map
    prob = ODEProblem(sys, [sys.Pin => experiment.p0, sys.Iin => 0.0], (0.0, experiment.tend))

    # Route to fast method passing callback explicitly
    return simulate(sys, prob, experiment, alg; reltol=reltol, abstol=abstol, callback=callback, kwargs...)
end

"""
    simulate(sys::ModelingToolkit.AbstractSystem, prob::ODEProblem, experiment::Experiment, alg=QNDF(); reltol=1e-6, abstol=1e-7, callback=nothing, kwargs...)

Accept pre-compiled ODE problem to bypass compilation overhead during benchmarks or repeated runs.
"""
function simulate(sys::ModelingToolkit.AbstractSystem, prob::ODEProblem, experiment::Experiment, alg=QNDF(); reltol=1e-5, abstol=1e-6, callback=nothing, kwargs...)
    
    # Stop simulation between every step when inputs are modified between experiment stages
    # Initialise integrator and apply optional callback for TMS control
    integrator = init(prob, alg; tstops=experiment.tstops, save_everystep=false, reltol=reltol, abstol=abstol, callback=callback, kwargs...)

    print("Simulating for: $(experiment.tend) seconds\n")
    
    # Update integrator after every experiment step and print time and voltage at each step
    for (i, step) in enumerate(experiment.steps)
        println("  [Step $i] Starting at t = $(integrator.t)s")
        step!(integrator, sys, step)
        
        # Dynamic single/multiple cell voltage readout
        v_val = try 
            integrator[sys.V] 
        catch 
            try integrator[sys.cell.v] catch; integrator[sys.cell1.v] end
        end
        
        println("  [Step $i] Finished at t = $(integrator.t)s | Voltage: $(round(v_val, digits=4))V")

        # Check if battery hit safety limit
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated
            println("  [!] Vmin/Vmax limit reached. Stopping experiment early.")
            break
        end
    end

    return integrator.sol
end