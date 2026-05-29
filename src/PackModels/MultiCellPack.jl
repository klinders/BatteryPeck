using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    MultiCellPack(; name, params=Chen2020(), config=(12,1), Qcell=5)

Create a multi-cell battery pack with flexible series/parallel configuration.

Builds a ModelingToolkit system representing a large battery pack with multiple cells in
series and parallel. Each cell is modeled using SPMe, and they are electrically connected
according to the specified configuration.

# Keyword Arguments
- `name`: System name for ModelingToolkit (required)
- `params::BatteryParameters`: Parameter set for all cells (default: Chen2020 parameters)
- `config::Tuple{Int,Int}`: Series and parallel configuration, `(n_series, n_parallel)` (default: 12 series, 1 parallel)
- `Qcell::Float64`: Nominal cell capacity in Ah (default: 5 Ah)

# Returns
- ModelingToolkit system with inputs (P, T) and outputs (V, I) for the entire pack

# Input/Output Signals
- Inputs: `P` (power), `T` (temperature)
- Outputs: `V` (pack voltage = sum of series cell voltages), `I` (pack current)

# Notes
- All cells share the same temperature input
- Total pack voltage is the sum of series-connected cell voltages
- Parallel cells share equal current

# Example
```julia
params = Chen2020()
pack = MultiCellPack(name=:battery, params=params, config=(96,3), Qcell=5)
```
"""
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