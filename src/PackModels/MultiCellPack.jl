using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

function MultiCellPack(;name, params=Chen2020(), config=(12,1), Qcell=5)

    @parameters begin
        t # Time variable
    end

    D = Differential(t)

    @variables begin 
        P(t)=0, [input=true]
        T(t)=298.15, [input=true]
        V(t)
        I(t)
    end

    Ncell = config[1]*config[2]

    cell = [SPMe(name=Symbol("cell_$i"), params=params) for i in 1:Ncell] # Battery model
    @named power = RealInput(guess=0)
    @named temp = RealInput(guess=298.15)
    @named source = Current()
    @named ground = Ground()

    eqs = [
        V ~ cell[1].p.v - cell[end].n.v
        I ~ cell[1].i
        D(P) ~ 0
        D(T) ~ 0
        power.u ~ P
        temp.u ~ T

        connect(source.n, cell[end].n)
        connect(source.p, cell[1].p)
        connect(ground.g, source.n)     
        [connect(temp, cell[i].T) for i in 1:Ncell]...
        [connect(cell[i].n, cell[i+1].p) for i in 1:Ncell-1]...

        source.I.u ~ power.u/V
    ]

    return System(eqs, t; systems=[cell..., power, temp, source, ground], name=name)
end