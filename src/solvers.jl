module Solvers

# using DifferentialEquations
# using ..PackModel
# using ..CellModel
# export solve, Model, Experiment, Simulation, Step, CurrentStep

# function constant_current(t)
#     return 2.0  # Constant discharge of 2A
# end


# struct Step
#     current::Function
#     period::Int
#     function Step(current::Function,period::Int)
#         new(current,period)
#     end
#     function Step(current::Number, period::Int)
#         new(x->current,period)
#     end
# end

# # struct CrateStep <: Step
# #     crate::Float64
# #     period::Float64
# # end

# # struct PowerStep <: Step
# #     power::Float64
# #     period::Float64
# # end

# function solve_model(du, u, p, t)
#     # Unpack the model parameters
#     model = p.model
#     if model.is_valid
#         # Call the model's objective function
#         model.objective(du, u, model.params, t)
#     else
#         error("Model is not valid")
#     end
# end



# Base.@kwdef mutable struct Experiment
#     steps::Vector{Step}
#     is_valid::Bool=false
#     function Experiment(a::Vector{Step})
#         # Do some validation
#         is_valid = true

#         # Set the values
#         x = new(a)
#         x.is_valid = is_valid

#         x
#     end
# end

# Base.@kwdef mutable struct Simulation
#     model::Model
#     experiment::Experiment
# end

# function solve(sim::Simulation, tspan=(0,100))
      
#     f = ODEFunction(sim.model.objective)

#     sim.model.params.I = constant_current# experiment_current(sim.experiment)

#     prob = ODEProblem(f, sim.model.params.u0, tspan, sim.model.params)
    
#     return DifferentialEquations.solve(prob)
# end

end  # module Solvers