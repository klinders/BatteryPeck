using CSV, Tables
using SciMLBase
using ModelingToolkit
# using Dates

abstract type Step end

"""
Rest for a period
"""
struct RestStep <: Step
    period::Real
end

"""
Charge up to specified SoC using the given power for the given period.
"""
struct ChargeStep <: Step
    soc::Real
    period::Real
    power::Real
    ChargeStep(soc::Real, period::Real=0, power::Real=11000) = new(soc,period,power)
end

"""
Apply a given  `power` for a given `period`
"""
struct PowerStep <: Step
    value::Real
    period::Real
end

"""
Apply a given  `current` for a given `period`
"""
struct CurrentStep <: Step
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
struct DriveStep <: Step
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

"Get the initial value from a power step"
function get_p0(s::PowerStep)
    return -s.value
end

"Get the initial value from a current step"
function get_p0(s::CurrentStep)
    return -s.value*4.2
end

"Get the initial value from a drive step"
function get_p0(s::DriveStep)
    return -s.csv[2][1]
end

"Get the initial value from a charge step"
function get_p0(s::ChargeStep)
    return s.power
end

"Get the initial value from a rest step"
function get_p0(s::RestStep)
    return 0
end


struct Experiment
    steps::Array{Step}
    tstops::Array{Float64}
    tend::Float64
    step_count::Int64
    p0::Float64
    # start_time::Date
    Experiment(steps::Vector{T} where T<:Step) = begin #, start_time::Date=Date(2020, 1, 1)
        tstops = cumsum([s.period for s in steps])
        tend = tstops[end]
        # Remove the last Tstop since it is the end of the simulation
        pop!(tstops)
        step_count = length(steps)
        p0 = get_p0(steps[1])
        return new(steps, tstops, tend, step_count, p0, start_time)
    end
end

function Base.:*(a::AbstractVector{<:Step}, n::Integer)
    return repeat(a,n)
end

function Base.:+(a::AbstractVector{<:Step}, b::AbstractVector{<:Step})
    return vcat(a,b)
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
