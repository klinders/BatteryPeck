# using ModelingToolkit

function ConstantTemperature(;name, T)
    @parameters t
    @named output = RealOutput()
    
    eqs = [
        output.u ~ T
    ]
    return System(eqs, t; systems=[output], name)
end