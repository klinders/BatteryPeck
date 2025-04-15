module CellModel

export battery_cell_ecm!, Parameters, CellParameters, LG50T
abstract type Parameters end

mutable struct CellParameters <: Parameters
    Q::Float64   # 1Ah battery (3600 Coulombs)
    R0::Float64    # Internal resistance
    R1::Float64    # RC resistance
    C1::Float64   # RC capacitance

    function CellParameters(params::Dict)

        set = new()
        for (key, value) in params
            # set field  of the struct with the name key to the corresponding value
            setproperty!(set, Symbol(key), value)
        end
        set
    end
end

LG50T = Dict{String, Any}(
    "Q"=>3600,
    "R0" => 0.2,
    "R1" => 0.1,
    "C1" => 350
)

function battery_cell_ecm!(du, u, p::CellParameters, t)
    SOC, V1 = u  # State variables: State of charge and voltage across R1C1
    I = p.I(t)  # Current input function (can be constant or time-varying)

    # Open circuit voltage as a simple linear function (can be a lookup table)
    Voc = 3.0 + 0.5 * SOC  # Example: linear approximation

    # Differential equations
    du[1] = -I / p.Q  # dSOC/dt
    du[2] = (I * p.R1 - V1) / (p.R1 * p.C1)  # dV1/dt

    return nothing
end

end  # module CellModel