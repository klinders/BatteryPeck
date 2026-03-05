# =====================================================================================================================
# solvers.jl
#
# Main simulation function
# Inputs:
# - sys = model to be simulated, e.g. "SingleCellPack()"
# - experiment = struct defining the simulation protocol
# - args = variable positional arguments
# - kwargs = keyword arguments for solver options of integrator, e.g. [saveat; abstol; reltol; dt]
#
# Creates static code from ODE expressions, handles multi-stage experiments, keeps track of simulation parameters
# =====================================================================================================================

# Import packages
using ModelingToolkit, OrdinaryDiffEq

function simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, args...; kwargs...)
    
    # Convert symbolic ODE expressions into static code
    # Inputs: ["equations"; "initial conditions"; "duration"]
    # [system equations from model; initialise P0 with experiment's value; set duration with experiment's duration] 
    prob = ODEProblem(sys, [sys.P=>experiment.p0], (0.0,experiment.tend))

    # Stop the simulation between every step, when inputs are modified between experiment stages, and keep track of simulation parameters
    # Inputs: ["static code of system equations"; "time of step transitions"; "saving configuration"; "solver options"]
    integrator = init(prob,args...; tstops=experiment.tstops, save_everystep=false, kwargs...)

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

    # without unit test: for step in experiment.steps (also remove print statements)

    # Force integrator to save last simulation point
    SciMLBase.savevalues!(integrator)

    # Return full history of simulation parameters
    return integrator.sol

end
