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

    while integrator.t - t_start < step.period
        soc = integrator.sol[sys.cell.soc][end]

        # Add hysteresis to prevent rapid switching
        if soc < end_soc*0.99
            set_u!(integrator, sys.Pin, step.power)
            set_u!(integrator, sys.Iin, 0)
        elseif soc > end_soc*1.01
            set_u!(integrator, sys.Pin, 0)
            set_u!(integrator, sys.Iin, 0)
        end
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, 60, true)

    end
end