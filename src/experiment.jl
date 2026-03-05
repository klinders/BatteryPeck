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

# Drive cycle (dis)charge profile, current-based instead of power-based for usage of Chen2020 data
# Inputs: ["current vs. time csv"; "time limit"]
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

# Extract initial power at t=0 for solver
get_p0(s::PowerStep)  = -s.value
get_p0(s::ChargeStep) = s.power
get_p0(s::CurrentStep)= -s.value*4                 # Assume V = 4V at t=0
get_p0(s::CurrentDriveStep) = -s.csv[2][1] * 4.0   # Assume V = 4V at t=0
get_p0(s::DriveStep)  = -s.csv[2][1]
get_p0(s::RestStep)   = 0.0

# Concatenate multiple steps into single instruction list
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

# Apply set power value (and optionally ambient temperature) for specific period during runtime
function apply_power!(integrator, sys, power, dt, T_amb=nothing)
    # Force update sys.P inside solver memory
    set_u!(integrator, sys.P, power)
    
    # Update ambient temperature if provided (e.g., Chen2020 dataset)
    if !isnothing(T_amb) && hasproperty(sys, :T_amb)
        set_u!(integrator, sys.T_amb, T_amb)
    end
    
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
    
# Update sys.P in integrator solver memory by converting CSV current into power using dynamic cell voltage
function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentDriveStep)
    # Check if ambient temperature data was loaded from CSV
    has_Tamb = length(step.csv) >= 3

    for idx in 1:length(step.csv[1])
        dt = step.csv[1][idx]
        current_val = step.csv[2][idx]
        T_amb_val = has_Tamb ? step.csv[3][idx] : nothing
        
        v = first(integrator[sys.V])
        
        # Also pass dynamic ambient temperature into apply_power
        apply_power!(integrator, sys, -current_val * v, dt, T_amb_val)
        
        # Stop applying power if solver terminates (e.g., due to event limit)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated
            break
        end
    end
end
