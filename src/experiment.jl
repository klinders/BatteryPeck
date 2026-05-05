# =====================================================================================================================
# experiment.jl
#
# Defines experimental protocol structures and execution logic for simulations.
# Provides [Rest; Power; Charge; Drive]Step primitives to construct load profiles,
# and updates/operates integrator solver accordingly using step! function.
# =====================================================================================================================

# Import packages
using CSV, Tables, DataFrames
using SciMLBase
using ModelingToolkit

# Define parent category "Step"
abstract type Step end

"""
Rest for a given `period`.
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
Apply a given `power` for a given `period`.
"""
struct PowerStep <: Step
    value::Real
    period::Real
end

"""
Apply a given `current` for a given `period`.
"""
struct CurrentStep <: Step
    value::Real
    period::Real
end

"""
Apply a drive cycle from the given csv.

**Arguments**
- `csv ::Vector{Any}` Path to the csv driving cycle
- `period ::Real` (optional) time to apply the cycle in seconds

At the moment, the period can be up to the length of the csv. 

TODO: Repeat the cycle when period is longer then csv.
"""
struct DriveStep <: Step
    csv::Vector{Any}
    period::Real

    # Constructor for pre-parsed data arrays
    DriveStep(csv::Vector{Any}, period::Real) = new(csv, period)

    DriveStep(file::String, period::Real=nothing) = begin
        # Read CSV and convert to matrix
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

"""
Drive cycle (dis)charge profile, current-based instead of power-based for usage of Chen2020 data.
Inputs: ["current vs. time csv"; "time limit"]
"""
struct CurrentDriveStep <: Step
    csv::Vector{Any}
    period::Real

    CurrentDriveStep(csv::Vector{Any}, period::Real) = new(csv, period)

    CurrentDriveStep(file::String, period=nothing) = begin
        f = CSV.read(file, DataFrame, skipto=15, header=14)
        
        t = f[!, "Test Time [s]"]
        raw_i = f[!, "Current [A]"]
        md = f[!, "Md"] # Mode column: "C" (charge), "D" (discharge), "R" (rest)

        # Extract chamber temperature and convert to Kelvin
        raw_T_amb = f[!, "Temperature Chamber [degC]"] .+ 273.15
        
        # Apply sign convention: charge = negative current, discharge = positive current
        signed_i = zeros(Float64, length(raw_i))
        for j in eachindex(raw_i)
            # Use strip() to remove any accidental whitespace from the CSV strings
            mode_str = strip(String(md[j])) 
            if mode_str == "C"
                signed_i[j] = -raw_i[j]
            elseif mode_str == "D"
                signed_i[j] = raw_i[j]
            else
                signed_i[j] = 0.0 # Rest
            end
        end

        # Calculate raw time steps
        raw_dt = diff(t)
        
        # Find valid rows where time actually moves forward (dt > 0)
        valid_indices = findall(x -> x > 0.0, raw_dt)
        
        if length(valid_indices) < length(raw_dt)
            bad_count = length(raw_dt) - length(valid_indices)
            println("  [Info] Cleaned $bad_count invalid CSV rows where dt <= 0")
        end
        
        clean_dt = raw_dt[valid_indices]
        clean_i  = signed_i[valid_indices]
        clean_t  = t[valid_indices .+ 1]            # Match shifted time indices
        clean_T_amb = raw_T_amb[valid_indices .+ 1] # ^

        tend = length(clean_dt)
        end_time = t[end]
        
        if !isnothing(period)
            cutoff_idx = findfirst(clean_t .>= period)
            if !isnothing(cutoff_idx)
                tend = cutoff_idx - 1
                end_time = clean_t[tend]
            end
        end

        # Pass cleaned current arrays to solver
        return new(Any[clean_dt[1:tend], clean_i[1:tend], clean_T_amb[1:tend]], end_time)
    end
end

"""
Get the initial value from a power step.
"""
function get_p0(s::PowerStep)
    return -s.value
end

"""
Get the initial value from a current step.
"""
function get_p0(s::CurrentStep)
    return -s.value * 4.2
end

"""
Get the initial value from a drive step.
"""
function get_p0(s::DriveStep)
    return -s.csv[2][1]
end

"""
Get the initial value from a currentdrive step.
"""
function get_p0(s::CurrentDriveStep)
    return -s.csv[2][1] * 4.2
end

"""
Get the initial value from a charge step.
"""
function get_p0(s::ChargeStep)
    return s.power
end

"""
Get the initial value from a rest step.
"""
function get_p0(s::RestStep)
    return 0
end

"""
Concatenate multiple steps into single instruction list.
Inputs: ["list of steps"]
"""
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
        # Remove last tstop; end of simulation
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

# Concatenate step sequences "a" and "b"
function Base.:+(a::AbstractVector{<:Step}, b::AbstractVector{<:Step})
    return [a;b]
end

"""
Update power in integrator solver memory and set current to zero.
"""
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::PowerStep)
    integrator.ps[sys.Pin] = -step.value
    integrator.ps[sys.Iin] = 0.0
    u_modified!(integrator, true)
    
    # Place a stop sign at the end of the period and step adaptively
    t_target = integrator.t + step.period
    SciMLBase.add_tstop!(integrator, t_target)
    while integrator.t < t_target
        SciMLBase.step!(integrator)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

"""
Clear power and current in integrator solver memory.
"""
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::RestStep)
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = 0.0
    u_modified!(integrator, true)
    
    t_target = integrator.t + step.period
    SciMLBase.add_tstop!(integrator, t_target)
    while integrator.t < t_target
        SciMLBase.step!(integrator)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

"""
Apply power in integrator solver memory until target SoC is reached.
"""
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::ChargeStep)
    soc = integrator.sol[sys.cell.soc][end]
    end_soc = step.soc
    t_start = integrator.t

    while (soc < end_soc) && (integrator.t - t_start < step.period)
        integrator.ps[sys.Pin] = step.power
        integrator.ps[sys.Iin] = 0.0
        u_modified!(integrator, true)
        
        # Advance adaptively in 60-second chunks to check SoC
        t_target = min(integrator.t + 60.0, t_start + step.period)
        SciMLBase.add_tstop!(integrator, t_target)
        while integrator.t < t_target
            SciMLBase.step!(integrator)
            if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
        end
        
        soc = integrator.sol[sys.cell.soc][end]
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

"""
Apply current in integrator solver memory until target period is reached.
"""
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentStep)
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = -step.value
    u_modified!(integrator, true)
    
    t_target = integrator.t + step.period
    SciMLBase.add_tstop!(integrator, t_target)
    while integrator.t < t_target
        SciMLBase.step!(integrator)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

"""
Apply power values from CSV profile until end is reached.
"""
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::DriveStep)
    has_Tamb = length(step.csv) >= 3

    for idx in 1:length(step.csv[1])
        dt = step.csv[1][idx]
        
        integrator.ps[sys.Pin] = step.csv[2][idx]
        integrator.ps[sys.Iin] = 0.0
        
        if has_Tamb && hasproperty(sys, :T_amb)
            integrator.ps[sys.T_amb] = step.csv[3][idx]
        end
        u_modified!(integrator, true)

        t_target = integrator.t + dt
        SciMLBase.add_tstop!(integrator, t_target)
        while integrator.t < t_target
            SciMLBase.step!(integrator)
            if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
        end
        
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

"""
Apply current values from CSV profile until end is reached.
"""
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentDriveStep)
    has_Tamb = length(step.csv) >= 3

    for idx in 1:length(step.csv[1])
        dt = step.csv[1][idx]
        
        integrator.ps[sys.Pin] = 0.0
        integrator.ps[sys.Iin] = -step.csv[2][idx]
        
        if has_Tamb && hasproperty(sys, :T_amb)
            integrator.ps[sys.T_amb] = step.csv[3][idx]
        end
        u_modified!(integrator, true)

        t_target = integrator.t + dt
        SciMLBase.add_tstop!(integrator, t_target)
        while integrator.t < t_target
            SciMLBase.step!(integrator)
            if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
        end
        
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end