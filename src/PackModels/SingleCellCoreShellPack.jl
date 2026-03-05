# =====================================================================================================================
# SingleCellCoreShellPack.jl
# Pack model for a single cell dynamically coupled to a core-shell LPTN
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal
using BatteryPeck

function SingleCellCoreShellPack(; name, params=Chen2020(), config=(1,1), Qcell=5, h_conv=15.0, T_ambient=298.15)
    # Surface area LGM50
    A_cell = 0.0053 
    
    @parameters begin
        t
    end
    D = Differential(t)

    @variables begin 
        P(t) = 0, [input=true]
        T_amb(t) = T_ambient, [input=true] # Ambient boundary temperature
        V(t), [guess=4.19*config[1]]
        I(t), [guess=0.0]
    end

    # Instantiate both domains
    @named cell = SPMe(params=params, Q=Qcell)
    @named thermal = CoreShellCell(T_start=T_ambient)
    
    # Instantiate electrical and thermal boundaries
    @named power = RealInput(guess=0)
    @named source = Current()
    @named ground = Ground()
    
    @named R_conv = ThermalResistor(R = 1 / (h_conv * A_cell))
    @named amb_temp_in = RealInput(guess=T_ambient)
    @named amb_source = PrescribedTemperature()

    eqs = [
        V ~ config[1]*cell.v
        I ~ config[2]*cell.i
        D(P) ~ 0
        D(T_amb) ~ 0
        
        power.u ~ P
        amb_temp_in.u ~ T_amb
        
        # Electrical connections
        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        source.I.u ~ power.u/config[1]/config[2]/cell.v
        
        # Thermal boundary connections (shell to ambient)
        connect(amb_temp_in, amb_source.T)
        connect(thermal.port_shell, R_conv.port_a)
        connect(R_conv.port_b, amb_source.port)
        
        # Multi-physics coupling
        # Feed Q_total (W/m^3) from SPMe to LPTN
        thermal.Q_volumetric_in.u ~ cell.Q_total
        # Feed T_core (K) from LPTN back to SPMe
        cell.T.u ~ thermal.core_cap.T
    ]

    return System(eqs, t; systems=[cell, thermal, power, source, ground, R_conv, amb_temp_in, amb_source], name=name)
end