# =====================================================================================================================
# SingleCellPack.jl
#
# Pack model for a single cell
# =====================================================================================================================

# Import packages
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

# Import module
using BatteryPeck

# Single cell model
# Inputs: ["component name"; "parameter set"; "series-parallel configuration"; "cell capacity"]
function SingleCellPack(; name, params=Chen2020(), config=(1,1), Qcell=5)
    # Independent variables
    @parameters begin
        t # Time variable
    end

    # Time derivative
    D = Differential(t)

    # Time-dependent I/O variables
    @variables begin 
        P(t) = 0, [input=true]
        T(t) = 298, [input=true]
        V(t), [guess=4.19*config[1]]
        I(t), [guess=0.0]
        
        Q_rev(t)
        Qᵢ(t)
        Qₑ(t)
        Qₛ(t)
        Qf(t)
        Q_total(t)
    end

    # Components
    @named cell = SPMe(params=params, Q=Qcell)
    @named power = RealInput(guess=0)
    @named temp = RealInput(guess=298)
    @named source = Current()
    @named ground = Ground()

    # Equations and connections
    eqs = [
        # Pack voltage
        V ~ config[1]*cell.v
        # Pack current
        I ~ config[2]*cell.i
        # Power and time are not affected by equations, only by user input
        D(P) ~ 0
        D(T) ~ 0
        # Link component to input variables
        power.u ~ P
        temp.u ~ T

        # Heat sources
        Q_rev ~ cell.Q_rev
        Qᵢ ~ cell.Qᵢ
        Qₑ ~ cell.Qₑ
        Qₛ ~ cell.Qₛ
        Qf ~ cell.Qf
        Q_total ~ cell.Q_total
        
        # Source- > cell-
        connect(source.n, cell.n)
        # Source+ > cell+
        connect(source.p, cell.p)
        # Ground > source- (reference)
        connect(ground.g, source.n)
        # Cell temp = ambient
        connect(temp, cell.T)

        # Applied current [I_pack / N_par = P_pack / (N_ser * V_cell)]
        source.I.u ~ power.u/config[1]/config[2]/cell.v
    ]

    return System(eqs, t; systems=[cell, power, temp, source, ground], name=name)
end