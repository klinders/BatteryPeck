# =====================================================================================================================
# CellThermal.jl
#
# 1D Core-Shell Lumped Parameter Thermal Network (LPTN) for the LGM50
# Based on experimental dataset by O'Regan et al. (2022)[1]:
# https://doi.org/10.1016/j.electacta.2022.140700
# Contains:
# 1. TemperatureDependentJellyroll: Generate temperature-dependent jellyroll component
# 2. CoreShellCell: Generate core-shell LPTN architecture component
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Thermal
using ModelingToolkitStandardLibrary.Blocks

"""
    TemperatureDependentJellyroll(; name, m_jelly=0.05708, T_start=300)

Generate temperature-dependent jellyroll component.

Contains heat source, thermal resistances, and thermal capacitance.

Equations:
Specific heat (C_p) uses third-order polynomials for constituent materials.
Heat capacity dynamics: dT/dt = Q_flow / (mass * C_p)

Editable values:
`m_jelly`: Alters total mass of jellyroll, affecting thermal inertia.
`T_start`: Sets initial thermal state.

# Arguments
- `name`: Component name for ModelingToolkit
- `m_jelly`: Total mass of jellyroll
- `T_start`: Initial temperature

# Returns
- Instantiated TemperatureDependentJellyroll component system
"""
@component function TemperatureDependentJellyroll(; name, m_jelly=0.05708, T_start=300)
    # Instantiate thermal port and declare required parameters and variables
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
    
    # Formulate specific heat capacity using constituent material polynomials and apply heat capacity dynamics
    # Mass fractions calculated from (Tab. 4)[1]
    # Third-order polynomials from (Tab. S4)[1]
    # MTK thermal convention: port.Q is heat flowing INTO component
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
        
        D(T) ~ port.Q_flow / (mass * C_p)
    ]
    
    return ODESystem(eqs, t, [T, C_p], [mass]; name=name, systems=[port])
end

"""
    CoreShellCell(; name, T_start=298.15)

Generate core-shell LPTN architecture component.

Models internal radial and axial heat transfer pathways.

Equations:
Solid cylinder with internal generation: R = 1 / (4 * pi * L * k)
Lumped average-mass node: R = 1 / (8 * pi * L * k)

Editable values:
`T_start`: Sets initial boundary temperature.

# Arguments
- `name`: Component name for ModelingToolkit
- `T_start`: Initial boundary temperature

# Returns
- Instantiated CoreShellCell component system
"""
@component function CoreShellCell(; name, T_start=298.15)
    # Define physical parameters and thermal properties for core and shell
    @parameters begin
        t
        
        # Jellyroll volume (Tab. 3)[1]
        V_jellyroll = 2.13e-5 
        
        # Harmonic mean of layers yields k_rad ≈ 1.13 W/mK
        # Halved from 1.07 to represent T_avg instead of from centre to edge
        R_rad_val = 0.535 
        
        # Volumetric mean yields k_ax ≈ 42 W/mK
        # Cell cylinder internal resistance = 1.21 K/W
        # Volume-averaged internal resistance = 1.21 / 2 = 0.605 K/W
        # Plastic insulator caps (0.22mm bottom, 0.2mm top) add ~1.62 K/W boundary bottleneck
        # Total axial resistance = 0.605 + 1.62 = 2.225 K/W
        R_ax_val = 2.225 
        
        # Shell thermal capacity
        # Mass = 10.64g
        # Specific heat (SS type 304) = 477 J/kgK (Tab. 7)[1]
        # C_th = 0.01064 * 477
        C_shell_val = 5.075 
    end
    
    # Instantiate internal thermal nodes, pathways, and external interfaces
    @named core_cap = TemperatureDependentJellyroll(T_start=T_start)
    @named shell_cap = HeatCapacitor(C=C_shell_val, T=T_start)
    
    @named R_rad = ThermalResistor(R=R_rad_val)
    @named R_ax = ThermalResistor(R=R_ax_val)
    
    @named heat_source = PrescribedHeatFlow()
    @named Q_volumetric_in = RealInput()
    
    @named port_shell = HeatPort()
    
    # Formulate routing connections between heat source, internal nodes, and shell boundary
    eqs = [
        heat_source.Q_flow.u ~ Q_volumetric_in.u * V_jellyroll,
        
        connect(heat_source.port, core_cap.port),
        
        connect(core_cap.port, R_rad.port_a),
        connect(R_rad.port_b, port_shell),
        
        connect(core_cap.port, R_ax.port_a),
        connect(R_ax.port_b, port_shell),
        
        connect(port_shell, shell_cap.port)
    ]
    
    subsystems = [core_cap, shell_cap, R_rad, R_ax, heat_source, Q_volumetric_in, port_shell]
    
    return ODESystem(eqs, t, [], [V_jellyroll, R_rad_val, R_ax_val, C_shell_val]; name=name, systems=subsystems)
end