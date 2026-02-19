# =====================================================================================================================
# experiment.jl
#
# Defines experimental protocol structures and execution logic for simulations.
# Provides [Rest; Power; Charge; Drive]Step primitives to construct load profiles,
# and updates/operates integrator solver accordingly using step! function.
# =====================================================================================================================

# Import packages
using CSV, Tables
using SciMLBase
using ModelingToolkit

# Define parent category "Step"
abstract type Step end

# Rest = battery does nothing
# Inputs: ["time"]
struct RestStep <: Step
    period::Real
end

# Charge with set power value until set SoC is reached
# Inputs: ["SoC", "time", "power"]
struct ChargeStep <: Step
    soc::Real
    period::Real
    power::Real

    # Provide default values (allows function to work with <3 input arguments)
    ChargeStep(soc::Real, period::Real=0, power::Real=11000) = new(soc,period,power)
end

# Constant power (dis)Charge
# Inputs: ["power", "time"]
struct PowerStep <: Step
    value::Real
    period::Real
end

# Constant current (dis)Charge
# Inputs: ["current", "time"]
struct CurrentStep <: Step
    value::Real
    period::Real
end

# Drive cycle (dis)charge profile
# Inputs: ["power vs. time csv"; "time"]
struct DriveStep <: Step
    csv::Vector{Any}
    period::Real

    DriveStep(file::String, period::Real=nothing) = begin
        # Read csv and convert to matrix
        f = CSV.File(file) |> Tables.matrix
        # Time = col1
        t = f[:,1]
        # Power = col2
        p = f[:,2]
        # Delta time
        dt = diff(t)
        # Total time
        tend = Int64(t[end]-1)
        # Abort (dis)charge profile at first time step after input period if [period != 0] and [period < tend]
        if !isnothing(period) && period < tend
            tend = findfirst(t.>=period)
        end

        # Only store data up to the cutoff point
        return new([dt[1:tend], p[1:tend]], tend)
    end
end

# Extract initial power at t=0 for solver
get_p0(s::PowerStep)  = -s.value
get_p0(s::ChargeStep) = s.power
get_p0(s::CurrentStep)= -s.value*4      # Assume V = 4V at t=0
get_p0(s::DriveStep)  = -s.csv[2][1]
get_p0(s::RestStep)   = 0.0

# Concatenate multiple steps into a single instruction list
# Inputs: ["list of steps"]
struct Experiment
    # (Dis)charge instructions
    steps::Array{Step}
    # Timestamps between steps
    tstops::Array{Float64}
    # Total time
    tend::Float64
    # Number of steps
    step_count::Int64
    # Initial power at t=0
    p0::Float64

    # Constructor able to accept any combination of steps
    Experiment(steps::Vector{T} where T<:Step) = begin
        # Timestamps between steps > tstop[n] = t[1] + ... + t[n-1] + t[n]
        tstops = cumsum([s.period for s in steps])
        tend = tstops[end]
        # Remove the last tstop since it is the end of the simulation
        pop!(tstops)
        step_count = length(steps)
        p0 = get_p0(steps[1])

        return new(steps, tstops, tend, step_count, p0)
    end
end

# Duplicates step sequence "a", "n" times
function Base.:*(a::AbstractVector{<:Step}, n::Integer)
    return repeat(a,n)
end

# Apply set power value for specific period during runtime
function apply_power!(integrator, sys, power, dt)
    # Force update sys.P inside solver memory
    set_u!(integrator, sys.P, power)
    # Notify solver about update
    u_modified!(integrator, true)
    # Run simulation for period "dt"
    OrdinaryDiffEq.step!(integrator, dt, true)
end

# Update sys.P in integrator solver memory, using power value and period from Rest/PowerStep
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::Union{RestStep, PowerStep})
    apply_power!(integrator, sys, get_p0(step), step.period)
end

# Keep updating sys.P in the integrator solver memory, using power values per timeframe from DriveStep profile, until end of csv is reached
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::DriveStep)
    for (dt, value) in zip(step.csv[1], step.csv[2])
        apply_power!(integrator, sys, -value, dt)

        # Stop applying power if solver terminates (e.g., due to event limit)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated
            break
        end
    end
end

# Update sys.P in integrator solver memory, and run until target SoC or target period from ChargeStep is reached
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::ChargeStep)
    soc = integrator.sol[sys.cell.soc][end]
    t_start = integrator.t

    while (soc < step.soc) && (integrator.t - t_start > step.period)
        apply_power!(integrator, sys, step.power, 60)
        soc = integrator.sol[sys.cell.soc][end]
        
        # Stop applying power if solver terminates (e.g., due to event limit)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated
            break
        end
    end
end

# Update sys.P in integrator solver memory, and run until target period from CurrentStep is reached
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentStep)
    t_start = integrator.t

    while (integrator.t - t_start) < step.period
        v = integrator.sol[sys.V][end]
        apply_power!(integrator, sys, -step.value*v, 1)

        # Stop applying power if solver terminates (e.g., due to event limit)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated
            break
        end
    end
end
    

