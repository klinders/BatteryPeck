module Solvers

using DifferentialEquations
using ..PackModel
export solve_ecm_battery

function constant_current(t)
    return -2.0  # Constant discharge of 2A
end

function solve_ecm_battery()
    return 0
    
    params = Dict(
        :num_cells => 3,
        :Q => 3600.0,   # 1Ah battery (3600 Coulombs)
        :R0 => 0.01,    # Internal resistance
        :R1 => 0.02,    # RC resistance
        :C1 => 500.0,   # RC capacitance
        :I => constant_current  # Current function
    )

    print(u0)
    u0 = vcat(fill([1.0, 0.0], params[:num_cells])...)  # Initial SOC=1, V1=0
    tspan = (0.0, 100.0)

    #f = ODEFunction(battery_pack_ecm!)

    #prob = ODEProblem(f, u0, tspan, params)
    
    return solve(prob)
end

end  # module Solvers