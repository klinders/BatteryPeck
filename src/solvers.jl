
using ModelingToolkit, OrdinaryDiffEq, CSV, Tables
# using ..PackModel
# using ..CellModel
# export solve, Model, Experiment, Simulation, Step, CurrentStep

# function constant_current(t)
#     return 2.0  # Constant discharge of 2A
# end

abstract type Step end

struct PowerStep <: Step
    value::Float64
    period::Float64
end

struct DriveStep <: Step
    csv::Vector{Any}
    period::Float64
    DriveStep(file::String, period::Float64=nothing) = begin
        f = CSV.File(file) |> Tables.matrix
        t = f[:,1]
        p = f[:,2]
        dt = diff(t)
        tend = t[end]
        if !isnothing(period) && period < tend
            tend = findfirst(t.>=period)
        end

        return new([dt[1:tend], p[1:tend]], tend)
    end
end

function step!(sys, integrator::SciMLBase.DEIntegrator, step::PowerStep)
    set_u!(integrator, sys.power.u, step.value)
    u_modified!(integrator, true)
    OrdinaryDiffEq.step!(integrator, step.period, true)
end

function step!(sys, integrator::SciMLBase.DEIntegrator, step::DriveStep)
    print("Stepping $(length(step.csv[1])) steps\n")
    for (dt, value) in zip(step.csv[1], step.csv[2])
        set_u!(integrator, sys.power.u, -value)
        u_modified!(integrator, true)
        OrdinaryDiffEq.step!(integrator, dt, true)
    end
end

function solve(sys::ModelingToolkit.System, experiment::Vector{Step})
    tend = sum(s.period for s in experiment)
    print("Simulating for: $tend seconds\n")

    prob = ODEProblem(sys, [sys.power.u=>0], (0.0, tend))
    
    integrator = init(prob, reltol = 1e-4, abstol = 1e-6, saveat=60.0)

    for s in experiment
        step!(sys, integrator, s)
    end

    return integrator.sol
end
