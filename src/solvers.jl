
using ModelingToolkit, OrdinaryDiffEq


function simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, args...; kwargs...)
    
    prob = ODEProblem(sys, [sys.Pin=>experiment.p0], (0.0,experiment.tend))
    integrator = init(prob,args...; tstops=experiment.tstops, save_everystep=false, kwargs...)

    print("Simulating for: $(experiment.tend) seconds\n")
    
    for step in experiment.steps
        step!(integrator, sys, step)
    end

    return integrator.sol
end
