
using ModelingToolkit, OrdinaryDiffEq, ProgressBars

"""
    simulate(sys::AbstractSystem, experiment::Experiment, args...; kwargs...)

Run a battery simulation using the provided ModelingToolkit system and experiment profile.

Simulates the battery system through all steps defined in the experiment, updating the solver
after each step and checking safety limits. Returns the solution object with time series data.

# Arguments
- `sys::AbstractSystem`: ModelingToolkit system (typically from `SPMe()` or pack models)
- `experiment::Experiment`: Experiment profile containing steps and parameters
- `args...`: Positional arguments passed to ODE solver (e.g., solver algorithm)
- `kwargs...`: Keyword arguments passed to ODE solver (e.g., `abstol`, `reltol`)

# Returns
- Solution object with fields `t` (time) and `u` (state variables) for plotting/analysis

# Example
```julia
params = OKane2022()
sys = SPMe(params=params)
exp = Experiment([PowerStep(1000, 3600)])  # 1000W for 1 hour
sol = simulate(sys, exp, Rodas4())
```
"""
function simulate(sys::ModelingToolkit.AbstractSystem, experiment::Experiment, args...; parameters=nothing, kwargs...)
    u0 = parameters!==nothing ? [[sys.Pin=>experiment.p0, sys.Iin=>0]; parameters] : [sys.Pin=>experiment.p0, sys.Iin=>0]
    prob = ODEProblem(sys, u0, (0.0,experiment.tend))
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
            println(integrator.sol.retcode)
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