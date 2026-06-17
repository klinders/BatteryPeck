# =====================================================================================================================
# experiment.jl
#
# Defines experimental protocol structures and execution logic for simulations.
# Provides [Rest; Power; Charge; Drive]Step primitives to construct load profiles.
# Includes both DAE step! assignments and an explicit state machine compiler.
# =====================================================================================================================

# Import packages
using CSV, Tables, DataFrames
using SciMLBase
using ModelingToolkit

export compile_experiment, Step, RestStep, ChargeStep, PowerStep, CurrentStep, DriveStep, CurrentDriveStep, Experiment

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
Apply a given `current` until a target SoC or Voltage is reached.
If limits are reached before `period` elapses, the system rests for the remainder of the period.
"""
struct TargetCurrentStep <: Step
    value::Real
    target_soc::Real
    target_v::Real
    period::Real
end

"""
Apply a drive cycle from the given csv.
"""
struct DriveStep <: Step
    csv::Vector{Any}
    period::Real

    DriveStep(csv::Vector{Any}, period::Real) = new(csv, period)

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
Drive cycle (dis)charge profile, current-based.
"""
struct CurrentDriveStep <: Step
    csv::Vector{Any}
    period::Real

    CurrentDriveStep(csv::Vector{Any}, period::Real) = new(csv, period)

    CurrentDriveStep(file::String, period=nothing) = begin
        f = CSV.read(file, DataFrame, skipto=15, header=14)
        
        t = f[!, "Test Time [s]"]
        raw_i = f[!, "Current [A]"]
        md = f[!, "Md"]

        raw_T_amb = f[!, "Temperature Chamber [degC]"] .+ 273.15
        signed_i = zeros(Float64, length(raw_i))
        for j in eachindex(raw_i)
            mode_str = strip(String(md[j])) 
            if mode_str == "C"
                signed_i[j] = -raw_i[j]
            elseif mode_str == "D"
                signed_i[j] = raw_i[j]
            else
                signed_i[j] = 0.0 
            end
        end

        raw_dt = diff(t)
        valid_indices = findall(x -> x > 0.0, raw_dt)
        
        clean_dt = raw_dt[valid_indices]
        clean_i  = signed_i[valid_indices]
        clean_t  = t[valid_indices .+ 1]            
        clean_T_amb = raw_T_amb[valid_indices .+ 1] 

        tend = length(clean_dt)
        end_time = t[end]
        
        if !isnothing(period)
            cutoff_idx = findfirst(clean_t .>= period)
            if !isnothing(cutoff_idx)
                tend = cutoff_idx - 1
                end_time = clean_t[tend]
            end
        end

        return new(Any[clean_dt[1:tend], clean_i[1:tend], clean_T_amb[1:tend]], end_time)
    end
end

function get_p0(s::PowerStep) return -s.value end
function get_p0(s::CurrentStep) return -s.value * 4.2 end
function get_p0(s::TargetCurrentStep) return -s.value * 4.2 end
function get_p0(s::DriveStep) return -s.csv[2][1] end
function get_p0(s::CurrentDriveStep) return -s.csv[2][1] * 4.2 end
function get_p0(s::ChargeStep) return s.power end
function get_p0(s::RestStep) return 0 end

struct Experiment
    steps::Array{Step}
    tstops::Array{Float64}
    tend::Float64
    step_count::Int64
    p0::Float64

    Experiment(steps::Vector{T} where T<:Step) = begin
        # For arbitrary duration ChargeSteps, tstops is just an estimate
        tstops = cumsum([s.period for s in steps])
        tend = tstops[end]
        pop!(tstops)
        step_count = length(steps)
        p0 = get_p0(steps[1])
        return new(steps, tstops, tend, step_count, p0)
    end
end

function Base.:*(a::AbstractVector{<:Step}, n::Integer) return repeat(a,n) end
function Base.:+(a::AbstractVector{<:Step}, b::AbstractVector{<:Step}) return [a;b] end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::PowerStep)
    integrator.ps[sys.Iin] = 0.0
    integrator.ps[sys.Pin] = -step.value
    u_modified!(integrator, true)
    
    t_target = integrator.t + step.period
    SciMLBase.add_tstop!(integrator, t_target)
    while integrator.t < t_target
        SciMLBase.step!(integrator)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

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

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::ChargeStep)
    get_soc() = try integrator.sol[sys.soc][end] catch; try integrator.sol[sys.cell.soc][end] catch; integrator.sol[sys.cell1.soc][end] end end
    soc = get_soc()
    t_start = integrator.t

    while (soc < step.soc) && (step.period == 0 || integrator.t - t_start < step.period)
        integrator.ps[sys.Pin] = step.power
        integrator.ps[sys.Iin] = 0.0
        u_modified!(integrator, true)
        
        # Advance adaptively in small chunks to monitor SoC crossing
        t_target = min(integrator.t + 10.0, step.period > 0 ? t_start + step.period : Inf)
        SciMLBase.add_tstop!(integrator, t_target)
        while integrator.t < t_target
            SciMLBase.step!(integrator)
            if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
        end
        
        soc = get_soc()
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

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

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::TargetCurrentStep)
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

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::DriveStep)
    has_Tamb = length(step.csv) >= 3
    for idx in 1:length(step.csv[1])
        dt = step.csv[1][idx]
        integrator.ps[sys.Pin] = step.csv[2][idx]
        integrator.ps[sys.Iin] = 0.0
        if has_Tamb && hasproperty(sys, :T_amb) integrator.ps[sys.T_amb] = step.csv[3][idx] end
        
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

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::CurrentDriveStep)
    has_Tamb = length(step.csv) >= 3
    for idx in 1:length(step.csv[1])
        dt = step.csv[1][idx]
        integrator.ps[sys.Pin] = 0.0
        integrator.ps[sys.Iin] = -step.csv[2][idx]
        if has_Tamb && hasproperty(sys, :T_amb) integrator.ps[sys.T_amb] = step.csv[3][idx] end
        
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
    compile_experiment(exp::Experiment)

Generates a dynamic state machine closure for explicit co-simulators, as oposed to monolithic ODE solvers.
Returns `get_load(t, v_pack, soc)`, which outputs `(current_load_A, is_finished_boolean)`.
"""
function compile_experiment(exp::Experiment)
    # State tracking variables
    step_idx = 1
    t_step_start = 0.0
    drive_sub_idx = 1
    step_limit_hit = false
    
    function get_load(t::Float64, v_pack::Float64, soc::Float64)
        if step_idx > length(exp.steps)
            return 0.0, true # Simulation complete
        end
        
        step = exp.steps[step_idx]
        dt = t - t_step_start
        
        # RestStep
        if step isa RestStep
            if dt >= step.period
                step_idx += 1; t_step_start = t; step_limit_hit = false
                return get_load(t, v_pack, soc)
            end
            return 0.0, false
            
        # CurrenStep
        elseif step isa CurrentStep
            if dt >= step.period
                step_idx += 1; t_step_start = t; step_limit_hit = false
                return get_load(t, v_pack, soc)
            end
            return -step.value, false

        # TargetCurrentStep
        elseif step isa TargetCurrentStep
            if dt >= step.period
                step_idx += 1; t_step_start = t; step_limit_hit = false
                return get_load(t, v_pack, soc)
            end
            
            # If limit is hit earlier in this step, force a Rest
            if step_limit_hit
                return 0.0, false 
            end
            
            I_total = -step.value
            is_charging = I_total > 0
            
            # Check limits
            hit_soc = is_charging ? (soc >= step.target_soc) : (soc <= step.target_soc)
            hit_v   = is_charging ? (v_pack >= step.target_v) : (v_pack <= step.target_v)
            
            if hit_soc || hit_v
                step_limit_hit = true # Lock into Rest mode for remainder of this period
                return 0.0, false
            end
            
            return I_total, false
            
        # PowerStep
        elseif step isa PowerStep
            if dt >= step.period
                step_idx += 1; t_step_start = t; step_limit_hit = false
                return get_load(t, v_pack, soc)
            end
            # Convention: step.value is discharge power
            return (-step.value / v_pack), false
            
        # ChargeStep
        elseif step isa ChargeStep
            if soc >= step.soc || (step.period > 0 && dt >= step.period)
                step_idx += 1; t_step_start = t; step_limit_hit = false
                return get_load(t, v_pack, soc)
            end
            # Charging uses negative current (I = -P / V)
            return -(step.power / v_pack), false
            
        # CurrentDriveStep
        elseif step isa CurrentDriveStep
            # Advance internal CSV index to match elapsed time
            while drive_sub_idx <= length(step.csv[1]) && dt >= step.csv[1][drive_sub_idx]
                dt -= step.csv[1][drive_sub_idx]
                t_step_start += step.csv[1][drive_sub_idx]
                drive_sub_idx += 1
            end
            
            if drive_sub_idx > length(step.csv[1]) || (!isnothing(step.period) && (t - (t_step_start - dt)) >= step.period)
                step_idx += 1; t_step_start = t; drive_sub_idx = 1; step_limit_hit = false
                return get_load(t, v_pack, soc)
            end
            return -step.csv[2][drive_sub_idx], false
            
        # Drive step
        elseif step isa DriveStep
            while drive_sub_idx <= length(step.csv[1]) && dt >= step.csv[1][drive_sub_idx]
                dt -= step.csv[1][drive_sub_idx]
                t_step_start += step.csv[1][drive_sub_idx]
                drive_sub_idx += 1
            end
            
            if drive_sub_idx > length(step.csv[1]) || (!isnothing(step.period) && (t - (t_step_start - dt)) >= step.period)
                step_idx += 1; t_step_start = t; drive_sub_idx = 1; step_limit_hit = false
                return get_load(t, v_pack, soc)
            end
            return (step.csv[2][drive_sub_idx] / v_pack), false
        end
        
        return 0.0, false
    end
    
    return get_load
end