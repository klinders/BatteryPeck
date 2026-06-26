# =====================================================================================================================
# experiment.jl
#
# Defines experimental protocol structures and execution logic for simulations.
# Provides [Rest; Power; Charge; Drive; FCR; SmartSoC]Step primitives to construct load profiles.
# Includes both DAE step! assignments and an explicit state machine compiler.
# Contains:
# 1. get_p0: Extract initial power or current value from step
# 2. step!: Execute single experimental step via ModelingToolkit integrator
# 3. compile_experiment: Generate dynamic state machine closure for explicit co-simulators
# =====================================================================================================================

using CSV, Tables, DataFrames
using SciMLBase
using ModelingToolkit

export compile_experiment, Step, RestStep, ChargeStep, PowerStep, CurrentStep, TargetCurrentStep, SmartSoCStep, DriveStep, CurrentDriveStep, FCRStep, Experiment

const GLOBAL_RAMP_DURATION = 1.0

abstract type Step end

"""
    RestStep

Rest for given period.
"""
struct RestStep <: Step
    period::Real
end

"""
    ChargeStep

Charge up to specified SoC using given power for given period.
"""
struct ChargeStep <: Step
    soc::Real
    period::Real
    power::Real
    ChargeStep(soc::Real, period::Real=0, power::Real=11000) = new(soc,period,power)
end

"""
    PowerStep

Apply given power for given period.
"""
struct PowerStep <: Step
    value::Real
    period::Real
end

"""
    CurrentStep

Apply given current for given period.
"""
struct CurrentStep <: Step
    value::Real
    period::Real
end

"""
    TargetCurrentStep

Apply given current until target SoC or voltage is reached (uses CV exponential decay).
"""
struct TargetCurrentStep <: Step
    value::Real
    target_soc::Real
    target_v::Real
    period::Real
end

"""
    SmartSoCStep

Dynamically charges or discharges at `max_current` to hit `target_soc`.
Smoothly tapers in the last 2%, and permanently snaps to 0.0A (Rest) once reached.
"""
struct SmartSoCStep <: Step
    max_current::Real
    target_soc::Real
    period::Real
end

"""
    DriveStep

Apply drive cycle from given csv.
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
    CurrentDriveStep

Drive cycle charge profile based on current array.
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
                signed_i[j] = raw_i[j]
            elseif mode_str == "D"
                signed_i[j] = -raw_i[j]
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

"""
    FCRStep

Apply vehicle to grid response based on frequency deviation.
"""
struct FCRStep <: Step
    t_data::Vector{Float64}
    f_data::Vector{Float64}
    max_power_w::Float64
    period::Real
end

"""
    get_p0(s::Step)

Extract initial power or current value from step.
"""
function get_p0(s::PowerStep) return s.value end
function get_p0(s::CurrentStep) return s.value * 4.2 end
function get_p0(s::TargetCurrentStep) return s.value * 4.2 end
function get_p0(s::SmartSoCStep) return s.max_current * 4.2 end
function get_p0(s::DriveStep) return s.csv[2][1] end
function get_p0(s::CurrentDriveStep) return s.csv[2][1] * 4.2 end
function get_p0(s::ChargeStep) return s.power end
function get_p0(s::RestStep) return 0.0 end
function get_p0(s::FCRStep) return 0.0 end

struct Experiment
    steps::Array{Step}
    tstops::Array{Float64}
    tend::Float64
    step_count::Int64
    p0::Float64

    Experiment(steps::Vector{T} where T<:Step) = begin
        tstops = cumsum([s.period for s in steps])
        tend = tstops[end]
        pop!(tstops)
        step_count = length(steps)
        p0 = get_p0(steps[1])
        return new(steps, tstops, tend, step_count, p0)
    end
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::PowerStep)
    integrator.ps[sys.Iin] = 0.0
    integrator.ps[sys.Pin] = step.value
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
    integrator.ps[sys.Iin] = step.value
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
    integrator.ps[sys.Iin] = step.value
    u_modified!(integrator, true)
    
    t_target = integrator.t + step.period
    SciMLBase.add_tstop!(integrator, t_target)
    while integrator.t < t_target
        SciMLBase.step!(integrator)
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
    end
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::SmartSoCStep)
    integrator.ps[sys.Pin] = 0.0
    integrator.ps[sys.Iin] = step.max_current
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
        integrator.ps[sys.Iin] = step.csv[2][idx]
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
    compile_experiment(protocol::Experiment)

Generate dynamic state machine closure for explicit co-simulators.
"""
function compile_experiment(protocol::Experiment)
    step_idx = 1
    t_step_start = 0.0
    drive_sub_idx = 1
    step_limit_hit = false
    t_limit_hit = 0.0
    last_step_val = 0.0
    
    function get_load(t::Float64, v_pack::Float64, soc::Float64)
        if step_idx > length(protocol.steps)
            return 0.0, true, :rest
        end
        
        step = protocol.steps[step_idx]
        dt = t - t_step_start
        
        raw_target = 0.0
        segment_period = step.period
        transition = false
        regime = :smooth
        
        if step isa RestStep
            raw_target = 0.0
            regime = :rest
            if dt >= step.period transition = true end
            
        elseif step isa CurrentStep
            raw_target = step.value
            if dt >= step.period transition = true end
            
        elseif step isa TargetCurrentStep
            raw_target = step.value
            
            if !step_limit_hit
                hit_soc = (raw_target > 0) ? (soc >= step.target_soc) : (soc <= step.target_soc)
                hit_v   = (raw_target > 0) ? (v_pack >= step.target_v) : (v_pack <= step.target_v)
                if hit_soc || hit_v 
                    step_limit_hit = true 
                    t_limit_hit = dt 
                end
            end
            
            if step_limit_hit 
                raw_target = raw_target * Base.exp(-(dt - t_limit_hit) / 600.0)
            end
            if dt >= step.period transition = true end

        elseif step isa SmartSoCStep
            error_soc = step.target_soc - soc
            
            # Snap strictly to Rest once hit
            if step_limit_hit
                raw_target = 0.0
                regime = :rest
            else
                if abs(error_soc) <= 0.0005 # 0.05% tolerance to perfectly nail the target
                    step_limit_hit = true
                    raw_target = 0.0
                    regime = :rest
                else
                    # Positive current = charge, Negative = discharge
                    raw_target = sign(error_soc) * step.max_current
                    
                    # Proportional taper for the last 2% of SoC to land smoothly without overshooting
                    if abs(error_soc) < 0.02
                        taper = abs(error_soc) / 0.02
                        raw_target *= clamp(taper, 0.05, 1.0) # Down to 5% of max current before snapping
                    end
                end
            end
            if dt >= step.period transition = true end
            
        elseif step isa PowerStep
            raw_target = step.value / v_pack
            if dt >= step.period transition = true end
            
        elseif step isa ChargeStep
            raw_target = (step.power / v_pack)
            if soc >= step.soc || (step.period > 0 && dt >= step.period) transition = true end
            
        elseif step isa FCRStep
            idx = searchsortedlast(step.t_data, dt)
            idx = clamp(idx, 1, length(step.f_data))
            freq = step.f_data[idx]
            
            delta_f = 50.0 - freq
            power_fraction = clamp(delta_f / 0.2, -1.0, 1.0)
            p_req = power_fraction * step.max_power_w
            
            raw_target = -p_req / v_pack
            segment_period = 10.0 
            regime = :grid
            if dt >= step.period transition = true end

        elseif step isa CurrentDriveStep || step isa DriveStep
            regime = :drive
            while drive_sub_idx <= length(step.csv[1]) && dt >= step.csv[1][drive_sub_idx]
                dt -= step.csv[1][drive_sub_idx]
                t_step_start += step.csv[1][drive_sub_idx]
                last_step_val = (step isa CurrentDriveStep) ? step.csv[2][drive_sub_idx] : (step.csv[2][drive_sub_idx] / v_pack)
                drive_sub_idx += 1
            end
            
            if drive_sub_idx > length(step.csv[1]) || (!isnothing(step.period) && (t - (t_step_start - dt)) >= step.period)
                transition = true
            else
                segment_period = step.csv[1][drive_sub_idx]
                raw_target = (step isa CurrentDriveStep) ? step.csv[2][drive_sub_idx] : (step.csv[2][drive_sub_idx] / v_pack)
            end
        end
        
        if transition
            last_step_val = raw_target
            step_idx += 1
            t_step_start = t
            drive_sub_idx = 1
            step_limit_hit = false
            t_limit_hit = 0.0
            return get_load(t, v_pack, soc)
        end
        
        actual_ramp = min(GLOBAL_RAMP_DURATION, segment_period / 2.0)
        if dt < actual_ramp && actual_ramp > 0.0
            ramped_target = last_step_val + (raw_target - last_step_val) * (dt / actual_ramp)
            return ramped_target, false, regime
        else
            return raw_target, false, regime
        end
    end
    
    return get_load
end