
using BatteryToolkit
import BatteryToolkit: step!

struct UncontrolledCharging <: BatteryToolkit.AbstractStep
    period::Real
    soc::Real
    power::Real
    UncontrolledCharging(period::Real=0, end_soc::Real=0.8, power::Real=11000) = new(period,end_soc,power)
end

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::UncontrolledCharging)
    soc = integrator.sol[sys.cell.soc][end]
    end_soc = step.soc
    t_start = integrator.t

    state = "Charging"
    set_u!(integrator, sys.Pin, step.power)
    set_u!(integrator, sys.Iin, 0)
    u_modified!(integrator, true)
    
    while integrator.t - t_start < step.period 
        soc = integrator.sol[sys.cell.soc][end]

        if state=="Charging" && soc >= end_soc
            state = "Idle"
            set_u!(integrator, sys.Pin, 0)
            set_u!(integrator, sys.Iin, 0)
            u_modified!(integrator, true)
        end

        OrdinaryDiffEq.step!(integrator, 60, true)

        if integrator.sol.retcode != SciMLBase.ReturnCode.Success
            return integrator.sol.retcode
        end
    end
end

struct DelayedCharging <: BatteryToolkit.AbstractStep
    period::Real
    soc::Real
    power::Real
    DelayedCharging(period::Real=0, end_soc::Real=0.8, power::Real=11000) = new(period,end_soc,power)
end

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::DelayedCharging)
    end_soc = step.soc
    t_start = integrator.t
    vehicle_capacity = 58 # kWh

    state = "Idle"
    set_u!(integrator, sys.Pin, 0)
    set_u!(integrator, sys.Iin, 0)
    u_modified!(integrator, true)


    while integrator.t - t_start < step.period 

        current_soc = integrator.sol[sys.cell.soc][end]

        if state != "End" && current_soc >= end_soc
            state = "End"
            @info "[$(integrator.t)] Reached end_soc, switching to idle mode"

            set_u!(integrator, sys.Pin, 0)
            set_u!(integrator, sys.Iin, 0)
            u_modified!(integrator, true)
        elseif state == "Idle" 
            # Check if we need to switch to charging
            t_left = step.period - (integrator.t - t_start)

            # time to charge to end_soc at max power
            soc_dt = (end_soc - current_soc)/(step.power/1000/3600/vehicle_capacity)
            
            if soc_dt > t_left - 3600 # Add some buffer to make sure we are charged in the end
                state = "Charging"
                @info "[$(integrator.t)] Switching to charging mode to reach end_soc"
                set_u!(integrator, sys.Pin, step.power)
                set_u!(integrator, sys.Iin, 0)
                u_modified!(integrator, true)
            end
        end

        OrdinaryDiffEq.step!(integrator, 60, true)

        if integrator.sol.retcode != SciMLBase.ReturnCode.Success
            return integrator.sol.retcode
        end
    end
end

struct FCR <: BatteryToolkit.AbstractStep
    period::Real
    soc::Real
    power::Real
    f0::Real
    dt::Int64
    droop::Real
    deadband::Real
    t_start::DateTime
    FCR(
        period::Real;
        end_soc::Real=0.8, 
        power::Real=11000, 
        f0::Real=50.0,
        dt::Int64=60,
        droop::Real=1/0.09, 
        deadband::Real=0.01
    ) = new(period, end_soc, power, f0, dt, droop, deadband, DateTime(0))
end

include("../data/DataSources.jl")
using .DataSources

global fcr_data = DataSources.GridFrequencySource()

function resample_to_dt(data::AbstractArray, step::AbstractStep)
    # step time in the data (in seconds)
    data_dt = ceil(Int,step.period / length(data))
    
    # Average over the dt to reduce the number of steps and make it more realistic
    n_windows = step.period ÷ step.dt
    
    if data_dt < step.dt
        # down sample
        window_size = step.dt ÷ data_dt
        return [sum(data[(i-1)*window_size+1 : i*window_size])/step.dt for i in 1:n_windows]
    elseif data_dt > step.dt
        window_size = data_dt ÷ step.dt
        return repeat(data, inner=window_size)
    else
        return data
    end
end

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::FCR, start_time::DateTime)
    end_soc = step.soc
    t_start = integrator.t
    vehicle_capacity = 58 # kWh

    # @info "Get FCR droop curve"
    global fcr_data
    df = DataSources.get_period(fcr_data, start_time, step.period)
    
    # Low frequency, give power
    @. df.power = step.droop * sign(df.frequency - step.f0) * max(abs(df.frequency - step.f0) - step.deadband, 0.0)*step.power
    # limit to the max
    df.power = clamp.(df.power, -step.power, step.power)
    
    power_dt = resample_to_dt(df.power, step)
    
    # @info "[$(integrator.t)] Start FCR service"
    state = "FCR"
    i = 1

    while integrator.t - t_start < step.period 

        if state == "FCR" 
            if i > length(power_dt)
                @warn "FCR step out of range ($(i)/$(length(power_dt)))"
                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)
            else
                set_u!(integrator, sys.Pin, power_dt[i])
                set_u!(integrator, sys.Iin, 0)
            end

            i += 1

            # Check if we need to switch to charging
            t_left = step.period - (integrator.t - t_start)
            current_soc = integrator.sol[sys.cell.soc][end]

            # time to charge to end_soc at max power
            soc_dt = (end_soc - current_soc)/(step.power/1000/3600/vehicle_capacity)
            
            if soc_dt > t_left - 3600 # Add some buffer to make sure we are charged in the end
                state = "Charging"
                # @info "[$(integrator.t)] Switching to charging mode to reach end_soc"
            end

        elseif state == "Charging"
            current_soc = integrator.sol[sys.cell.soc][end]

            if current_soc >= end_soc
                state = "Idle"
                # @info "[$(integrator.t)] Reached end_soc, switching to idle mode"

                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)
            else
                set_u!(integrator, sys.Pin, step.power)
                set_u!(integrator, sys.Iin, 0)
            end
        elseif state == "Idle"
            set_u!(integrator, sys.Pin, 0)
            set_u!(integrator, sys.Iin, 0)
        end
        
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, step.dt, true)

        if integrator.sol.retcode != SciMLBase.ReturnCode.Success
            return integrator.sol.retcode
        end
    end
end

struct aFRR <: BatteryToolkit.AbstractStep
    period::Real
    soc::Real
    power::Real
    dt::Int64
    soc_max::Real
    soc_min::Real
    activation::Real
    t_start::DateTime
    aFRR(
        period::Real;
        end_soc::Real=0.8, 
        power::Real=11000,
        dt::Int64=900, 
        soc_max::Real=0.8,
        soc_min::Real=0.2,
        activation::Real=0.05
    ) = new(period, end_soc, power, dt,soc_max,soc_min, activation, DateTime(0))
end

global afrr_data = DataSources.ActivatedBalancingSource()

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::aFRR, start_time::DateTime)
    end_soc = step.soc
    t_start = integrator.t
    vehicle_capacity = 58 # kWh

    # @info "Get aFRR activated and available"
    global afrr_data
    df = DataSources.get_period(afrr_data, start_time, step.period)
    
    # Low frequency, give power
    df.up_pct   = df.up_activated   ./ df.up_accepted
    df.down_pct = df.down_activated ./ df.down_accepted

    df.power = ifelse.(
        (df.up_pct .> df.down_pct) .& (df.up_pct .> step.activation),   -step.power,
    ifelse.(
        (df.down_pct .> df.up_pct) .& (df.down_pct .> step.activation), step.power,
        0.0))

    # limit to the max
    df.power = clamp.(df.power, -step.power, step.power)
    
    # @info "[$(integrator.t)] Start FCR service"
    state = "aFRR"
    i = 1

    while integrator.t - t_start < step.period 

        if state == "aFRR" 
            current_soc = integrator.sol[sys.cell.soc][end]

            if i > length(df.power)
                @warn "aFRR step out of range ($(i)/$(length(df.power)))"
                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)

            elseif current_soc > step.soc_max && df.power[i] > 0 || current_soc < step.soc_min && df.power[i] < 0
                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)
            else
                set_u!(integrator, sys.Pin, df.power[i])
                set_u!(integrator, sys.Iin, 0)
            end

            i += 1

            # Check if we need to switch to charging
            t_left = step.period - (integrator.t - t_start)

            # time to charge to end_soc at max power
            soc_dt = (end_soc - current_soc)/(step.power/1000/3600/vehicle_capacity)
            
            if soc_dt > t_left - 3600 # Add some buffer to make sure we are charged in the end
                state = "Charging"
                # @info "[$(integrator.t)] Switching to charging mode to reach end_soc"
            end

        elseif state == "Charging"
            current_soc = integrator.sol[sys.cell.soc][end]

            if current_soc >= end_soc
                state = "Idle"
                # @info "[$(integrator.t)] Reached end_soc, switching to idle mode"

                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)
            else
                set_u!(integrator, sys.Pin, step.power)
                set_u!(integrator, sys.Iin, 0)
            end
        elseif state == "Idle"
            set_u!(integrator, sys.Pin, 0)
            set_u!(integrator, sys.Iin, 0)
        end
        
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, step.dt, true)

        if integrator.sol.retcode != SciMLBase.ReturnCode.Success
            return integrator.sol.retcode
        end
    end
end

struct SolarArbitrage <: BatteryToolkit.AbstractStep
    period::Real
    soc::Real
    power::Real
    dt::Int64
    soc_max::Real
    soc_min::Real
    t_start::DateTime
    SolarArbitrage(
        period::Real;
        end_soc::Real=0.8, 
        power::Real=11000,
        dt::Int64=60, 
        soc_max::Real=0.8,
        soc_min::Real=0.2,
    ) = new(period, end_soc, power, dt,soc_max,soc_min, DateTime(0))
end

global solar_data = DataSources.SolarPowerSource()
global day_ahead_data = DataSources.DayAheadSource()

using JuMP
using Statistics
using DataFrames

function BatteryToolkit.step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::SolarArbitrage, start_time::DateTime)
    end_soc = step.soc
    t_start = integrator.t
    vehicle_capacity = 58 # kWh

    # @info "Get aFRR activated and available"
    global solar_data
    global day_ahead_prices

    df_solar = DataSources.get_period(solar_data, start_time, step.period)
    df_dayahead = DataSources.get_period(day_ahead_data, start_time, step.period)

    power = DataFrame([])

    # if there is solar available: charge
    # During low prices

    q75 = quantile(df_dayahead.price, .75)
    q25 = quantile(df_dayahead.price, .25)

    # High price start discharging
    power = ifelse.(df_dayahead.price .> q75, -step.power, ifelse.(df_dayahead.price .< q25, step.power, 0))

    # limit to the max
    power = clamp.(power, -step.power, step.power)

    power = resample_to_dt(power, step)

    # @info "[$(integrator.t)] Start FCR service"
    state = "Arbitrage"
    i = 1

    while integrator.t - t_start < step.period 

        if state == "Arbitrage" 
            current_soc = integrator.sol[sys.cell.soc][end]

            if i > length(power)
                @warn "SolarArbitrage step out of range ($(i)/$(length(power)))"
                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)

            elseif current_soc > step.soc_max && power[i] > 0 || current_soc < step.soc_min && power[i] < 0
                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)
            else
                set_u!(integrator, sys.Pin, power[i])
                set_u!(integrator, sys.Iin, 0)
            end

            i += 1

            # Check if we need to switch to charging
            t_left = step.period - (integrator.t - t_start)

            # time to charge to end_soc at max power
            soc_dt = (end_soc - current_soc)/(step.power/1000/3600/vehicle_capacity)
            
            if soc_dt > t_left - 3600 # Add some buffer to make sure we are charged in the end
                state = "Charging"
                # @info "[$(integrator.t/3600)] Switching to charging mode to reach end_soc"
            end

        elseif state == "Charging"
            current_soc = integrator.sol[sys.cell.soc][end]

            if current_soc >= end_soc
                state = "Idle"
                # @info "[$(integrator.t)] Reached end_soc, switching to idle mode"

                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)
            else
                set_u!(integrator, sys.Pin, step.power)
                set_u!(integrator, sys.Iin, 0)
            end
        elseif state == "Idle"
            set_u!(integrator, sys.Pin, 0)
            set_u!(integrator, sys.Iin, 0)
        end
        
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, step.dt, true)

        if integrator.sol.retcode != SciMLBase.ReturnCode.Success
            return integrator.sol.retcode
        end
    end
end