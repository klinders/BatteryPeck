using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

using BatteryPeck

function SingleCellPack(;name, params=Chen2020(), config=(96,3), Qcell=5)

    @parameters
        t # Time variable
    end

    D = Differential(t)

    @variables begin 
        P(t)=0, [input=true]
        T(t)=298, [input=true]
        V(t)
        I(t)
    end

    @named cell = SPMe(params=params, Q=Qcell) # Battery cell model
    @named power = RealInput(guess=0)
    @named temp = RealInput(guess=298)
    @named source = Current()
    @named ground = Ground()

    eqs = [
        V ~ config[1]*cell.v
        I ~ config[2]*cell.i
        D(P) ~ 0
        D(T) ~ 0
        power.u ~ P
        temp.u ~ T

        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        connect(temp, cell.T)

        source.I.u ~ power.u/config[1]/config[2]/cell.v
    ]

    return System(eqs, t; systems=[cell, power, temp, source, ground], name=name)
end