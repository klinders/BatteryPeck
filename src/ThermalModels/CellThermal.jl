# =====================================================================================================================
# CellThermal.jl
#
# 1D Core-Shell Lumped Parameter Thermal Network (LPTN) for the LGM50
# Based on experimental dataset by O'Regan et al. (2022)[1]:
# https://doi.org/10.1016/j.electacta.2022.140700
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Thermal
using ModelingToolkitStandardLibrary.Blocks

"""
    TemperatureDependentJellyroll(; name, m_jelly=0.05708, T_start=300)

Temperature-dependent jellyroll component containing heat source, thermal resistances, and thermal capacitance.

Equations:
Specific heat (C_p) uses third-order polynomials for constituent materials.
Heat capacity dynamics: dT/dt = Q_flow / (mass * C_p)

Editable values:
`m_jelly`: Alters total mass of jellyroll, affecting thermal inertia.
`T_start`: Sets initial thermal state.
"""
@component function TemperatureDependentJellyroll(; name, m_jelly=0.05708, T_start=300)
    @named port = HeatPort()
    
    @parameters begin
        t
        mass = m_jelly
    end
    
    @variables begin
        T(t) = T_start
        C_p(t)
    end
    
    D = Differential(t)
    
    # Specific heat equations
    # Mass fractions calculated from (Tab. 4)[1]
    # Third-order polynomials from (Tab. S4)[1]
    eqs = [
        T ~ port.T,
        
        C_p ~ (
            0.4571 * (-0.0008414*T^3 + 0.7892*T^2 - 241.3*T + 2.508e4) +  # Positive electrode
            0.2892 * ( 0.0004932*T^3 - 0.4910*T^2 + 169.4*T - 1.897e4) +  # Negative electrode
            0.0343 * ( 0.0014940*T^3 - 1.4440*T^2 + 475.5*T - 5.130e4) +  # Separator
            0.0503 * ( 4.503e-6*T^3  - 0.006256*T^2 + 3.281*T + 355.7) +  # Aluminium foil
            0.1025 * ( 1.445e-6*T^3  - 0.001946*T^2 + 0.9633*T + 236.0) + # Copper foil
            0.0666 * 229.0                                                # Electrolyte (Tab. 7)[1]
        ),
        
        # Heat capacity dynamics (Q = m * c_p * dT/dt)
        # MTK thermal convention: port.Q is heat flowing INTO component
        D(T) ~ port.Q_flow / (mass * C_p)
    ]
    
    return ODESystem(eqs, t, [T, C_p], [mass]; name=name, systems=[port])
end

"""
    CoreShellCell(; name, T_start=298.15)

Core-shell LPTN architecture component.
Models internal radial and axial heat transfer pathways.

Equations:
Solid cylinder with internal generation: R = 1 / (4 * pi * L * k)
Lumped average-mass node: R = 1 / (8 * pi * L * k)

Editable values:
`T_start`: Sets initial boundary temperature.
"""
@component function CoreShellCell(; name, T_start=298.15)
    @parameters begin
        t
        # Jellyroll volume (Tab. 3)[1]
        V_jellyroll = 2.13e-5 
        
        # Harmonic mean of layers yields k_rad ≈ 1.13 W/mK
        # Solid cylinder with internal generation: R = 1 / (4 * pi * L * k)
         #R_rad_val = 1.07 

        # Halved from 1.07 to represent T_avg instead of from centre to edge
         R_rad_val = 0.535 
        
        # Volumetric mean yields k_ax ≈ 42 W/mK
        # Cell cylinder resistance = 1.21 K/W
        # Plastic insulator caps (0.22mm bottom, 0.2mm top) add ~1.62 K/W bottleneck
        R_ax_val = 2.83 
        
        # Mass = 10.64g
        # Specific heat (SS type 304) = 477 J/kgK (Tab. 7)[1]
        # C_th = 0.01064 * 477
        C_shell_val = 5.075 
    end
    
    # Instantiate nodes
    @named core_cap = TemperatureDependentJellyroll(T_start=T_start)
    @named shell_cap = HeatCapacitor(C=C_shell_val, T=T_start)
    
    # Instantiate internal pathways
    @named R_rad = ThermalResistor(R=R_rad_val)
    @named R_ax = ThermalResistor(R=R_ax_val)
    
    # Instantiate heat source interface
    @named heat_source = PrescribedHeatFlow()
    @named Q_volumetric_in = RealInput()
    
    # Absolute outer surface of cell
    @named port_shell = HeatPort()
    
    eqs = [
        # Convert volumetric heat (W/m^3) to absolute heat (W)
        heat_source.Q_flow.u ~ Q_volumetric_in.u * V_jellyroll,
        
        # Inject generated heat directly into core node
        connect(heat_source.port, core_cap.port),
        
        # Connect core to shell via radial pathway
        connect(core_cap.port, R_rad.port_a),
        connect(R_rad.port_b, port_shell),
        
        # Connect core to shell via axial pathway (parallel to radial)
        connect(core_cap.port, R_ax.port_a),
        connect(R_ax.port_b, port_shell),
        
        # Attach shell thermal mass to outer node 
        # Implicitly connects it to thermal ground
        connect(port_shell, shell_cap.port)
    ]
    
    subsystems = [core_cap, shell_cap, R_rad, R_ax, heat_source, Q_volumetric_in, port_shell]
    
    return ODESystem(eqs, t, [], [V_jellyroll, R_rad_val, R_ax_val, C_shell_val]; name=name, systems=subsystems)
end