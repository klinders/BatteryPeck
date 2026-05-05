using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
Construct single cell pack model.
"""
function SingleCellPack(;name, params=Chen2020(), config=(96,3), Qcell=5)

    @parameters begin
        t # Time variable
        
        # Define boundary inputs as parameters
        Pin = 0.0
        Iin = 0.0
        Tin = 298.15
    end

    @variables begin 
        V(t)
        I(t)
    end

    @named cell = SPMe(params=params, Q=Qcell) # Battery cell model
    @named power = RealInput(guess=0)
    @named current = RealInput(guess=0)
    @named temp = RealInput(guess=298)
    @named source = Current()
    @named ground = Ground()

    eqs = [
        V ~ config[1]*cell.v
        I ~ config[2]*cell.i
        power.u ~ Pin
        current.u ~ Iin
        temp.u ~ Tin

        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        connect(temp, cell.T)

        source.I.u ~ power.u/config[1]/config[2]/cell.v + current.u/config[2]
    ]

    return System(eqs, t; systems=[cell, power, temp, source, ground], name=name)
end