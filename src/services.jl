struct UncontrolledCharging <: Step
    period::Real
    soc::Real
    power::Real
    UncontrolledCharging(period::Real=0, end_soc::Real=0.8, power::Real=11000) = new(period,end_soc,power)
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::UncontrolledCharging)
    soc = integrator.sol[sys.cell.soc][end]
    end_soc = step.soc
    t_start = integrator.t

    current = integrator.sol[sys.cell.i][end]

    while integrator.t - t_start < step.period 
        soc = integrator.sol[sys.cell.soc][end]

        # Add hysteresis to prevent rapid switching
        if soc < end_soc*0.99
            new_current = round(integrator.sol[sys.cell.i][end])
            if new_current != current
                current = new_current
            end
            set_u!(integrator, sys.Pin, step.power)
            set_u!(integrator, sys.Iin, 0)
        elseif soc > end_soc*1.01
            set_u!(integrator, sys.Pin, 0)
            set_u!(integrator, sys.Iin, 0)
        end
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, 60, true)

        if integrator.sol.retcode != SciMLBase.ReturnCode.Success
            return integrator.sol.retcode
        end
    end
end

struct FCR <: Step
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

include("../data/GridFrequency/GridFrequency.jl")

function average_power(df, dt)
    n_windows = nrow(df) ÷ dt
    times  = [df.t[i*dt] for i in 1:n_windows]          # timestamp at end of each window
    freqs = [sum(df.f[(i-1)*dt+1 : i*dt])/dt for i in 1:n_windows] # average frequency in each window
    powers = [sum(df.p[(i-1)*dt+1 : i*dt])/dt for i in 1:n_windows]
    return DataFrame(t=times, f=freqs, p=powers)
end

function step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::FCR, start_time::DateTime)
    soc = integrator.sol[sys.cell.soc][end]
    end_soc = step.soc
    t_start = integrator.t
    vehicle_capacity = 58 # kWh
    current = integrator.sol[sys.cell.i][end]

    @info "Get FCR droop curve"
    df = get_frequency(start_time, step.period)
    
    # Low frequency, give power
    @. df.p = step.droop * sign(df.f - step.f0) * max(abs(df.f - step.f0) - step.deadband, 0.0)*step.power
    
    # limit to the max
    df.p = clamp.(df.p, -step.power, step.power)
    
    # Average over the dt to reduce the number of steps and make it more realistic
    df = average_power(df, step.dt)
    
    @info "[$(integrator.t)] Start FCR service"
    state = "FCR"
    i = 1

    while integrator.t - t_start < step.period 

        if state == "FCR" 
            if i > nrow(df)
                set_u!(integrator, sys.Pin, 0)
                set_u!(integrator, sys.Iin, 0)
            else
                set_u!(integrator, sys.Pin, df.p[i])
                set_u!(integrator, sys.Iin, 0)
            end

            i += 1

            # Check if we need to switch to charging
            t_left = step.period - (integrator.t - t_start)
            current_soc = integrator.sol[sys.cell.soc][end]

            # time to charge to end_soc at max power
            dt = (end_soc - current_soc)/(step.power/1000/3600/vehicle_capacity)
            
            if dt > t_left - 3600 # Add some buffer to make sure we are charged in the end
                state = "Charging"
                @info "[$(integrator.t)] Switching to charging mode to reach end_soc"
            end

        elseif state == "Charging"
            current_soc = integrator.sol[sys.cell.soc][end]

            if current_soc >= end_soc
                state = "Idle"
                @info "[$(integrator.t)] Reached end_soc, switching to idle mode"

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