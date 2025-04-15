module Solvers

using DifferentialEquations
using ..PackModel
using ..CellModel
export solve, Model, Experiment, Simulation, Step, CurrentStep

function constant_current(t)
    return 2.0  # Constant discharge of 2A
end

abstract type Step end

struct CurrentStep <: Step
    current::Float64
    period::Float64
end

struct CrateStep <: Step
    crate::Float64
    period::Float64
end

struct PowerStep <: Step
    power::Float64
    period::Float64
end

Base.@kwdef mutable struct Model
    objective::Function
    submodels::Vector{Model}
    params::Parameters
    is_valid::Bool
    function Model(objective::Function, submodels::Vector{Model}, params::Parameters)

        # Do some validation
        is_valid = true
        
        # Validate submodels
        for m in submodels
            if !m.is_valid
                is_valid = false
            end
        end

        # Set the values
        new(objective, submodels, params, is_valid)
    end
end

Base.@kwdef mutable struct Experiment
    steps::Vector{Step}
    is_valid::Bool
    function Experiment(a::Vector{Step})
        # Do some validation
        is_valid = true

        # Set the values
        x = new(a)
        x.is_valid = is_valid
    end
end

Base.@kwdef mutable struct Simulation
    model::Model
    experiment::Experiment
end

function solve(sim::Simulation, tspan=(0,100))
      
    f = ODEFunction(sim.model.objective)

    prob = ODEProblem(f, sim.model._u0, tspan, sim.model.params)
    
    return solve(prob)
end

end  # module Solvers