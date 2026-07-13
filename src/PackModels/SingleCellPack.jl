using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    SingleCellPack(; name, params=Chen2020(), config=(96,3), Qcell=5)

Create a single battery cell pack with configurable series/parallel arrangement.

Builds a ModelingToolkit system representing a battery pack of series-parallel connected cells.
Each cell is modeled using the SPMe (Single Particle Model with Electrolyte) framework.

# Keyword Arguments
- `name`: System name for ModelingToolkit (required)
- `params::BatteryParameters`: Parameter set for the cell model (default: Chen2020 parameters)
- `config::Tuple{Int,Int}`: Series and parallel configuration, `(n_series, n_parallel)` (default: 96 series, 3 parallel)
- `Qcell::Float64`: Nominal cell capacity in Ah (default: 5 Ah)

# Returns
- ModelingToolkit system with inputs (Pin, Iin, Tin) and outputs (V, I)

# Input/Output Signals
- Inputs: `Pin` (power), `Iin` (current), `Tin` (temperature)
- Outputs: `V` (voltage), `I` (current)

# Example
```julia
params = OKane2022()
pack = SingleCellPack(name=:pack, params=params, config=(96,3), Qcell=5)
```
"""
function SingleCellPack(;name, params=Chen2020(), config=(96,3), Qcell=5, T=298.15, kargs...)

    @parameters begin
        t # Time variable
    end

    D = Differential(t)

    @variables begin 
        Pin(t)=0, [input=true]
        Iin(t)=0, [input=true]
        Tin(t)=T, [input=true]
        V(t)
        I(t)
    end

    @named cell = SPMe(params=params, Q=Qcell) # Battery cell model
    @named power = RealInput(guess=0)
    @named current = RealInput(guess=0)
    @named temp = RealInput(guess=T)
    @named source = Current()
    @named ground = Ground()

    eqs = [
        V ~ config[1]*cell.v
        I ~ config[2]*cell.i
        D(Pin) ~ 0
        D(Iin) ~ 0
        D(Tin) ~ 0
        power.u ~ Pin
        current.u ~ Iin
        temp.u ~ Tin

        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        connect(temp, cell.T)

        source.I.u ~ power.u/config[1]/config[2]/cell.v + current.u/config[2]
    ]

    return System(eqs, t; systems=[cell, power, temp, source, ground], name=name, kargs...)
end