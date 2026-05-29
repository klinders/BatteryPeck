using CSV, Tables
using SciMLBase
using ModelingToolkit
using Dates

"""
Abstract type for all step types in the experiment.
"""
abstract type AbstractStep end

"""
Rest for a period
"""
struct RestStep <: AbstractStep
    period::Real
end

"""
Charge up to specified SoC using the given power for the given period.
"""
struct ChargeStep <: AbstractStep
    soc::Real
    period::Real
    power::Real
    ChargeStep(soc::Real, period::Real=0, power::Real=11000) = new(soc,period,power)
end

"""
Apply a given  `power` for a given `period`
"""
struct PowerStep <: AbstractStep
    value::Real
    period::Real
end

"""
Apply a given  `current` for a given `period`
"""
struct CurrentStep <: AbstractStep
    value::Real
    period::Real
end

"""
Apply a drivecycle from the given csv

**Arguments**
- `csv ::Vector{Any}` Path to the csv driving cycle
- `period ::Real` (optional) time to apply the cycle in seconds

At the moment, the period can be up to the lenght of the csv. 

TODO: Repeat the cycle when period is longer then csv
"""
struct DriveStep <: AbstractStep
    csv::Vector{Any}
    period::Real
    DriveStep(file::String, period::Real=nothing) = begin
        f = CSV.File(file) |> Tables.matrix
        t = f[:,1]
        p = f[:,2]
        dt = diff(t)
        tend = Int64(t[end]-1)
        if !isnothing(period) && period < tend
            tend = findfirst(t.>=period)
        end

        return new([dt[1:tend], p[1:tend]], tend)
    end
end

"""
Get the initial value for the first step of the experiment. This is used to set the initial conditions for the simulation. The behavior depends on the specific step type (PowerStep, RestStep, ChargeStep, CurrentStep, DriveStep).
"""
function get_p0(s::AbstractStep)
    error("get_p0 not implemented for step type $(typeof(s))")
end

function get_p0(s::PowerStep)
    return -s.value
end

function get_p0(s::CurrentStep)
    return -s.value*4.2
end

function get_p0(s::DriveStep)
    return -s.csv[2][1]
end

function get_p0(s::ChargeStep)
    return s.power
end

function get_p0(s::RestStep)
    return 0
end



"""
    Experiment(steps::Vector{<:AbstractStep}, start_time::DateTime=DateTime(2020, 1, 1))

Create an experimental profile composing multiple battery operation steps.

Combines a sequence of operation steps (power, current, rest, charge, drive cycle) into
a single experiment. Automatically calculates step timing and prepares parameters for
simulation with the `simulate()` function.

# Arguments
- `steps::Vector{<:AbstractStep}`: Vector of step objects (RestStep, PowerStep, CurrentStep, ChargeStep, DriveStep)
- `start_time::DateTime`: Real-world timestamp for first step (default: 2020-01-01)

# Fields (automatically calculated)
- `steps::Vector`: Original step vector
- `tstops::Vector{Float64}`: Cumulative time at end of each step except the last (s)
- `tend::Float64}`: Total simulation duration (s)
- `step_count::Int64`: Number of steps
- `p0::Float64`: Initial power/current value (W or A)
- `start_time::DateTime`: Experiment start timestamp

# Example
```julia
steps = [
    PowerStep(1000, 1800),      # 1000W for 30 min
    RestStep(300),               # 5 min rest
    PowerStep(-500, 3600)        # -500W (discharge) for 1 hour
]
exp = Experiment(steps)
sol = simulate(sys, exp, Rodas4())
```
"""
struct Experiment
    steps::Array{AbstractStep}
    tstops::Array{Float64}
    tend::Float64
    step_count::Int64
    p0::Float64
    start_time::DateTime
    Experiment(steps::Vector{T} where T<:AbstractStep, start_time::DateTime=DateTime(2020, 1, 1)) = begin
        tstops = cumsum([s.period for s in steps])
        tend = tstops[end]
        # Remove the last Tstop since it is the end of the simulation
        pop!(tstops)
        step_count = length(steps)
        p0 = get_p0(steps[1])
        return new(steps, tstops, tend, step_count, p0, start_time)
    end
end

function Base.:*(a::AbstractVector{<:AbstractStep}, n::Integer)
    return repeat(a,n)
end

function Base.:+(a::AbstractVector{<:AbstractStep}, b::AbstractVector{<:AbstractStep})
    return vcat(a,b)
end

"""
    step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::AbstractStep)
Apply a single step of the experiment to the DE integrator. This function modifies the integrator's input parameters according to the step type
and advances the simulation by the step's period. The behavior depends on the specific step type (PowerStep, RestStep, ChargeStep, CurrentStep, DriveStep).

# Arguments
- `integrator::SciMLBase.DEIntegrator`: The DE integrator to modify and step
- `sys::ModelingToolkit.AbstractSystem`: The system being simulated, used to access input variables
- `step::AbstractStep`: The step to apply, which determines how the integrator's inputs are modified and how long to step the simulation


"""
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::AbstractStep)
    error("step! not implemented for step type $(typeof(step))")
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::PowerStep)
    set_u!(integrator, sys.Pin, -step.value)
    set_u!(integrator, sys.Iin, 0)
    u_modified!(integrator, true)
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::RestStep)
    set_u!(integrator, sys.Pin, 0)
    set_u!(integrator, sys.Iin, 0)
    u_modified!(integrator, true)
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::ChargeStep)
    soc = integrator.sol[sys.cell.soc][end]
    end_soc = step.soc
    t_start = integrator.t

    while soc < end_soc && integrator.t - t_start < step.period
        set_u!(integrator, sys.Pin, step.power)
        set_u!(integrator, sys.Iin, 0)
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, 60, true)
        soc = integrator.sol[sys.cell.soc][end]
    end
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentStep)
    t_start = integrator.t

    while integrator.t - t_start < step.period
        v = integrator.sol[sys.V][end]
        set_u!(integrator, sys.Pin, 0)
        set_u!(integrator, sys.Iin, -step.value)
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, 1, true)
    end
end

function step!(integrator::SciMLBase.DEIntegrator,sys::ModelingToolkit.AbstractSystem, step::DriveStep)
    # print("Stepping $(length(step.csv[1])) steps\n")
    for (dt, value) in zip(step.csv[1], step.csv[2])
        set_u!(integrator, sys.Pin, -value)
        set_u!(integrator, sys.Iin, 0)
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, dt, true)
    end
end
