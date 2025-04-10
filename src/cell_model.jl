module CellModel

export battery_cell_ecm!

function battery_cell_ecm!(du, u, p, t)
    SOC, V1 = u  # State variables: State of charge and voltage across R1C1
    I = p[:I](t)  # Current input function (can be constant or time-varying)
    Q, R0, R1, C1 = p[:Q], p[:R0], p[:R1], p[:C1]

    # Open circuit voltage as a simple linear function (can be a lookup table)
    Voc = 3.0 + 0.5 * SOC  # Example: linear approximation

    # Differential equations
    du[1] = -I / Q  # dSOC/dt
    du[2] = (I * R1 - V1) / (R1 * C1)  # dV1/dt

    return nothing
end

end  # module CellModel