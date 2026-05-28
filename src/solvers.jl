
using ModelingToolkit, OrdinaryDiffEq, ProgressBars


function simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, args...; kwargs...)
    
    prob = ODEProblem(sys, [sys.Pin=>experiment.p0], (0.0,experiment.tend))
    integrator = init(prob,args...; tstops=experiment.tstops, save_everystep=false, kwargs...)

    time_scale, time_unit, time_symbol = format_time(experiment.tend)

    print("Simulating for: $(round(experiment.tend*time_scale,digits=2)) $(time_unit)\n")
    
    # Update integrator after every experiment step and print time and voltage at each step
    for (i, step) in ProgressBar(enumerate(experiment.steps))

        # add the start time if the step wants it
        if hasproperty(step, :t_start)
            t_start = experiment.start_time + Second(integrator.t)
            step!(integrator, sys, step, t_start)
        else
            step!(integrator, sys, step)
        end

                
        # Check if battery hit safety limit
        if integrator.sol.retcode != SciMLBase.ReturnCode.Success
            println("[!] Simulation aborted at t = $(integrator.t)s. Stopping experiment early.")
            break
        end

    end

    return integrator.sol
end

function format_time(seconds::Real)
    if seconds < 60
        return 1,"seconds", "s"
    elseif seconds < 3600
        return 1/60,"minutes", "m"
    elseif seconds < 3600*24
        1/3600, "hours", "h"
    else
        return 1/(3600*24), "days", "d"
    end
end