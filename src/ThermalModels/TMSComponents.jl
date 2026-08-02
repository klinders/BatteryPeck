# ==============================================================================
# TMSComponents.jl
# 1D fluid and pipe wall elements using MTK
# Reference: Kakac, S., Shah, R.K. and Aung, W., 1987. Handbook of single-phase convective heat transfer. Wiley.
# Contains:
# 1. FluidPort: Define connector for 1D fluid network
# 2. PipeWallNode: Model lumped thermal mass of aluminium cooling channel wall
# 3. FluidNode: Model 1D advection and pressure drop of coolant as ordinary differential equation
# 4. ConvectionModel: Act as thermal bridge calculating dynamic convective heat flux into fluid
# 5. TMSNode: Couple PipeWallNode, ConvectionModel, and FluidNode into container component
# ==============================================================================

export FluidPort, PipeWallNode, FluidNode, ConvectionModel, MinorLoss, TMSNode

using ModelingToolkit
using ModelingToolkitStandardLibrary.Thermal: HeatPort 
using IfElse

@parameters t

"""
    FluidPort(; name)

Define connector for 1D fluid network.

Contains pressure (p), mass flow rate (m_flow), and temperature (T).

# Arguments
- `name`: Component name for ModelingToolkit

# Returns
- Instantiated FluidPort connector system
"""
@connector function FluidPort(; name)
    # Declare state variables for fluid port
    sts = @variables begin
        p(t) 
        m_flow(t), [connect = Flow] 
        T(t) 
    end
    
    return ODESystem(Equation[], t, sts, []; name=name)
end

"""
    PipeWallNode(; name, params::PackParameters, length::Float64)

Model lumped thermal mass of aluminium cooling channel wall.

Equations:
Heat capacity dynamics: m * c_p * dT/dt = Q_in + Q_out

Editable values:
`T(t)`: Initial temperature of wall. Alters starting boundary condition.

# Arguments
- `name`: Component name for ModelingToolkit
- `params::PackParameters`: Struct containing physical and geometric properties
- `length::Float64`: Length of pipe segment

# Returns
- Instantiated PipeWallNode component system
"""
@component function PipeWallNode(; name, params::PackParameters, length::Float64)
    # Create external thermal ports and extract geometric dimensions
    @named port = HeatPort()
    @named fluid_port = HeatPort() 
    
    W = params.tms_geometry.channel_width
    H = params.tms_geometry.channel_height
    th = params.tms_geometry.wall_thickness
    
    # Calculate physical volumes and compute lumped thermal capacitance
    outer_area = (W + 2*th) * (H + 2*th)
    inner_area = W * H
    wall_volume = (outer_area - inner_area) * length
    C_wall = wall_volume * params.pipe_wall.density * params.pipe_wall.specific_heat
    
    # Declare initial state variables
    sts = @variables T(t)=params.ambient_temperature
    
    # Formulate governing energy balance equations for thermal mass
    eqs = [
        T ~ port.T,
        T ~ fluid_port.T,
        C_wall * Differential(t)(T) ~ port.Q_flow + fluid_port.Q_flow
    ]
    
    return ODESystem(eqs, t, sts, []; systems=[port, fluid_port], name=name)
end

"""
    FluidNode(; name, params::PackParameters, length::Float64)

Model 1D advection and pressure drop of coolant as ordinary differential equation.

Fast compilation ensured by tracking fluid mass dynamics explicitly.

Equations:
Pressure drop (Darcy-Weisbach): Δp = f * (L / D_h) * (ρ * v² / 2)
Advection dynamics: m * c_p * dT/dt = Q_conv + m_dot * c_p * (T_in - T_out)
Curvature stabilisation (Kakac, p. 5.9): Re_crit = 2100 * (1 + 12 * (R_c / a)^-0.5)
Laminar friction multiplier (Kakac, p. 5.7): f_c / f_s = 0.1125 * De^0.5 
Turbulent friction (Kakac, p. 5.22): f_c * (R_c / a)^0.5 = 0.00725 + 0.076 * (Re * (R_c / a)^-2)^-0.25

Editable values:
`T(t)`: Initial temperature of fluid. Alters starting thermal state.

# Arguments
- `name`: Component name for ModelingToolkit
- `params::PackParameters`: Struct containing physical and geometric properties
- `length::Float64`: Length of pipe segment

# Returns
- Instantiated FluidNode component system
"""
@component function FluidNode(; name, params::PackParameters, length::Float64)
    # Create required ports and extract thermophysical fluid properties
    @named port_a = FluidPort()
    @named port_b = FluidPort()
    @named heat_port = HeatPort() 
    
    rho = params.fluid.density
    cp  = params.fluid.specific_heat
    mu  = params.fluid.dynamic_viscosity
    k   = params.fluid.thermal_conductivity
    
    # Calculate hydraulic parameters and radius of curvature for serpentine channel
    W = params.tms_geometry.channel_width
    H = params.tms_geometry.channel_height
    N = params.tms_geometry.number_of_channels
    
    A_flow = (W * H) * N
    P_wet = 2 * (W + H) * N
    D_h = 4 * A_flow / P_wet
    a_hyd = D_h / 2.0 
    R_c = 0.0125 
    fluid_mass = A_flow * length * rho
    
    # Declare fluid state variables
    sts = @variables begin
        T(t)=params.ambient_temperature
        p(t) 
        m_flow(t) 
        v(t) 
        Re(t) 
        De(t) 
        Re_crit(t) 
        f_curved_ratio(t) 
        f_major(t) 
        dp(t) 
    end
    
    # Formulate mass conservation, curvature stabilisation, hydraulic friction, and convective energy balance equations
    eqs = [
        0 ~ port_a.m_flow + port_b.m_flow,
        m_flow ~ port_a.m_flow, 
        
        v ~ m_flow / (rho * A_flow),
        Re ~ (rho * max(v, 1e-6) * D_h) / mu,
        De ~ Re * sqrt(a_hyd / R_c),
        Re_crit ~ 2100.0 * (1.0 + 12.0 * (R_c / a_hyd)^-0.5),
        
        f_curved_ratio ~ IfElse.ifelse(De < 30.0, 
                            1.0, 
                            IfElse.ifelse(De < 300.0, 
                                0.419 * De^0.275, 
                                0.1125 * sqrt(De)
                            )
                         ),
                         
        f_major ~ IfElse.ifelse(Re < Re_crit, 
            (1.0 / Re) * f_curved_ratio,
            4.0 * (sqrt(a_hyd / R_c) * (0.00725 + 0.076 * (Re * (a_hyd / R_c)^2)^-0.25))
        ),
        
        dp ~ (f_major * (length / D_h) * (rho * v^2 / 2.0)) + ((42 * length) * (rho * v^2 / 2.0)),
                           
        port_a.p - port_b.p ~ dp,
        p ~ (port_a.p + port_b.p) / 2.0, 
        
        heat_port.T ~ T,
        port_b.T ~ T, 
        
        fluid_mass * cp * Differential(t)(T) ~ heat_port.Q_flow + m_flow * cp * (port_a.T - T)
    ]
    
    return ODESystem(eqs, t, sts, []; systems=[port_a, port_b, heat_port], name=name)
end

"""
    ConvectionModel(; name, params::PackParameters, length::Float64)

Act as thermal bridge calculating dynamic convective heat flux into fluid.

Equations:
Convective heat transfer (Gnielinski correlation):
Nu = ((f/8) * (Re - 1000) * Pr) / (1 + 12.7 * √(f/8) * (Pr^(2/3) - 1))

# Arguments
- `name`: Component name for ModelingToolkit
- `params::PackParameters`: Struct containing physical and geometric properties
- `length::Float64`: Length of pipe segment

# Returns
- Instantiated ConvectionModel component system
"""
@component function ConvectionModel(; name, params::PackParameters, length::Float64)
    # Create thermal interaction ports and extract physical fluid properties
    @named solid_port = HeatPort()
    @named fluid_port = HeatPort()
    
    k  = params.fluid.thermal_conductivity
    cp = params.fluid.specific_heat
    mu = params.fluid.dynamic_viscosity
    
    # Compute hydraulic geometry parameters and convective surface area
    W = params.tms_geometry.channel_width
    H = params.tms_geometry.channel_height
    N = params.tms_geometry.number_of_channels
    
    A_flow = (W * H) * N
    P_wet = 2 * (W + H) * N
    D_h = 4 * A_flow / P_wet
    
    A_surface = P_wet * length
    Pr = (cp * mu) / k
    
    # Declare state variables for heat transfer calculations
    sts = @variables begin
        Re(t) 
        Nu(t) 
        h(t) 
        Q_flow(t) 
    end
    
    # Formulate nusselt number correlations and calculate total convective heat flux
    eqs = [
        Nu ~ IfElse.ifelse(Re < 2300, 
                           8.23, 
                           let f_turb = (0.79 * log(max(Re, 2300.0)) - 1.64)^-2
                            ((f_turb/8.0) * (Re - 1000.0) * Pr) / (1.0 + 12.7 * sqrt(f_turb/8.0) * (Pr^(2/3) - 1.0))
                           end),
                           
        h ~ (Nu * k) / D_h,
        Q_flow ~ h * A_surface * (solid_port.T - fluid_port.T),
        
        solid_port.Q_flow ~ Q_flow,
        fluid_port.Q_flow ~ -Q_flow
    ]
    
    return ODESystem(eqs, t, sts, []; systems=[solid_port, fluid_port], name=name)
end

"""
    TMSNode(; name, params::PackParameters, length::Float64)

Couple PipeWallNode, ConvectionModel, and FluidNode into container component.

# Arguments
- `name`: Component name for ModelingToolkit
- `params::PackParameters`: Struct containing physical and geometric properties
- `length::Float64`: Length of pipe segment

# Returns
- Instantiated TMSNode coupled component system
"""
@component function TMSNode(; name, params::PackParameters, length::Float64)
    # Initialise internal subcomponents and external routing ports
    @named wall = PipeWallNode(params=params, length=length)
    @named fluid = FluidNode(params=params, length=length)
    @named convection = ConvectionModel(params=params, length=length)
    
    @named port_a = FluidPort()
    @named port_b = FluidPort()
    @named heat_port = HeatPort() 
    
    # Define routing connections between fluid domains, pipe wall, and convection bridge
    eqs = [
        connect(port_a, fluid.port_a),
        connect(port_b, fluid.port_b),
        connect(heat_port, wall.port),
        
        connect(wall.fluid_port, convection.solid_port),
        connect(fluid.heat_port, convection.fluid_port),
        
        convection.Re ~ fluid.Re
    ]
    
    return ODESystem(eqs, t, [], []; systems=[wall, fluid, convection, port_a, port_b, heat_port], name=name)
end