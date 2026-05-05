using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
Construct multi-cell pack in series.
"""
function MultiCellPack(;name, params=Chen2020(), config=(12,1), Qcell=5)

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

    Ncell = config[1]*config[2]

    cell = [SPMe(name=Symbol("cell_$i"), params=params) for i in 1:Ncell] # Battery model
    @named power = RealInput(guess=0)
    @named current = RealInput(guess=0)
    @named temp = RealInput(guess=298)
    @named source = Current()
    @named ground = Ground()

    eqs = [
        V ~ cell[1].p.v - cell[end].n.v
        I ~ cell[1].i
        power.u ~ Pin
        current.u ~ Iin
        temp.u ~ Tin

        connect(source.n, cell[end].n)
        connect(source.p, cell[1].p)
        connect(ground.g, source.n)     
        [connect(temp, cell[i].T) for i in 1:Ncell]...
        [connect(cell[i].n, cell[i+1].p) for i in 1:Ncell-1]...

        source.I.u ~ power.u/V + current.u
    ]

    return System(eqs, t; systems=[cell..., power, current, temp, source, ground], name=name)
end