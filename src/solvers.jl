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

# Import packages
using ModelingToolkit, OrdinaryDiffEq

# Added alg=QNDF() as default, and explicit default tolerances
function simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, alg=QNDF(); reltol=1e-4, abstol=1e-7, kwargs...)
    
    # Convert symbolic ODE expressions into static code
    # Added sparse=true to automatically leverage Sparse AutoDiff for speed
    prob = ODEProblem(sys, [sys.P=>experiment.p0], (0.0,experiment.tend), sparse=true)

    # Stop simulation between every step, when inputs are modified between experiment stages
    # Pass reltol and abstol variables into init function
    integrator = init(prob, alg; tstops=experiment.tstops, save_everystep=false, reltol=reltol, abstol=abstol, kwargs...)

    print("Simulating for: $(experiment.tend) seconds\n")
    
    # Update integrator after every experiment step, and print time and voltage at each step
    for (i, step) in enumerate(experiment.steps)
        println("  [Step $i] Starting at t = $(integrator.t)s")
        step!(integrator, sys, step)
        println("  [Step $i] Finished at t = $(integrator.t)s | Current Voltage: $(integrator[sys.cell.v])V")

        # Check if battery hit safety limit
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated
            println("  [!] Vmin/Vmax limit reached. Stopping experiment early.")
            break
        end
    end

    # Force integrator to save last simulation point
    SciMLBase.savevalues!(integrator, true)
    
    return integrator.sol
end