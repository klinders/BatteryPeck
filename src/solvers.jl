module Solvers

using DifferentialEquations
using ..PackModel
using ..CellModel
export solve, Model, Experiment, Simulation, Step, CurrentStep

function constant_current(t)
    return 2.0  # Constant discharge of 2A
end


struct Step
    current::Function
    period::Int
    function Step(c::Function,p::Int)
        new(c,p)
    end
    function Step(c::Number, p::Int)
        new(x->c,p)
    end
end

# struct CrateStep <: Step
#     crate::Float64
#     period::Float64
# end

# struct PowerStep <: Step
#     power::Float64
#     period::Float64
# end

Base.@kwdef mutable struct Model
    objective::Function
    params::Parameters
    submodels::Vector{Model}=nothing
    is_valid::Bool=false
    function Model(objective::Function, params::Parameters, submodels::Vector{Model}=nothing)

        # Do some validation
        is_valid = true
        
        # Validate submodels
        for m in submodels
            if !m.is_valid
                is_valid = false
            end
        end

        # Set the values
        x = new(objective, params, submodels)
        x.is_valid = is_valid

        x
    end
end

Base.@kwdef mutable struct Experiment
    steps::Vector{Step}
    is_valid::Bool=false
    function Experiment(a::Vector{Step})
        # Do some validation
        is_valid = true

        # Set the values
        x = new(a)
        x.is_valid = is_valid

        x
    end
end

Base.@kwdef mutable struct Simulation
    model::Model
    experiment::Experiment
end

function solve(sim::Simulation, tspan=(0,100))
      
    f = ODEFunction(sim.model.objective)

    sim.model.params.I = sim.experiment.steps[1].current

    prob = ODEProblem(f, sim.model.params.u0, tspan, sim.model.params)
    
    return DifferentialEquations.solve(prob)
end

end  # module Solvers