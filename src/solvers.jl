
using ModelingToolkit, OrdinaryDiffEq, ProgressBars


function simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, args...; kwargs...)
    
    prob = ODEProblem(sys, [sys.Pin=>experiment.p0], (0.0,experiment.tend))
    integrator = init(prob,args...; tstops=experiment.tstops, save_everystep=false, kwargs...)

    time_scale, time_unit, time_symbol = format_time(experiment.tend)

    print("Simulating for: $(experiment.tend*time_scale) $(time_unit)\n")
    
    for step in ProgressBar(experiment.steps)
        step!(integrator, sys, step)
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