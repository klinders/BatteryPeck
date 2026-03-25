# ==============================================================================
# TMSComponents.jl
# 1D fluid and pipe wall elements using MTK
# ==============================================================================

export FluidPort, PipeWallNode, FluidNode, ConvectionModel, MinorLoss, TMSNode

using ModelingToolkit
using ModelingToolkitStandardLibrary.Thermal: HeatPort 
using IfElse

@parameters t

"""
    FluidPort(; name)

Defines connector for 1D fluid network.
Contains pressure (p), mass flow rate (m_flow), and temperature (T).
"""
@connector function FluidPort(; name)
    # Declare state variables
    sts = @variables begin
        p(t) 
        m_flow(t), [connect = Flow] 
        T(t) 
    end
    
    # Return ODE system for fluid port
    ODESystem(Equation[], t, sts, []; name=name)
end

"""
    PipeWallNode(; name, params::PackParameters, length::Float64)

Models lumped thermal mass of aluminium cooling channel wall.

Equations:
Heat capacity dynamics: m * c_p * dT/dt = Q_in + Q_out

Editable values:
`T(t)`: Initial temperature of wall. Alters starting boundary condition.
"""
@component function PipeWallNode(; name, params::PackParameters, length::Float64)
    # Create external thermal ports
    @named port = HeatPort()
    @named fluid_port = HeatPort() 
    
    # Extract geometric dimensions
    W = params.tms_geometry.channel_width
    H = params.tms_geometry.channel_height
    th = params.tms_geometry.wall_thickness
    
    # Calculate cross sectional areas
    outer_area = (W + 2*th) * (H + 2*th)
    inner_area = W * H
    
    # Calculate total wall volume
    wall_volume = (outer_area - inner_area) * length
    
    # Compute lumped thermal capacitance
    C_wall = wall_volume * params.pipe_wall.density * params.pipe_wall.specific_heat
    
    # Set initial wall temperature
    sts = @variables T(t)=298.15 
    
    # Define governing equations
    eqs = [
        T ~ port.T,
        T ~ fluid_port.T,
        
        # Apply energy balance for thermal mass
        C_wall * Differential(t)(T) ~ port.Q_flow + fluid_port.Q_flow
    ]
    
    # Return ODE system for pipe wall
    ODESystem(eqs, t, sts, []; systems=[port, fluid_port], name=name)
end

"""
    FluidNode(; name, params::PackParameters, length::Float64)

Models 1D advection and pressure drop of coolant as ordinary differential equation.
Fast compilation ensured by tracking fluid mass dynamics explicitly.

Equations:
Pressure drop (Darcy-Weisbach): Δp = f * (L / D_h) * (ρ * v² / 2)
Advection dynamics: m * c_p * dT/dt = Q_conv + m_dot * c_p * (T_in - T_out)

Editable values:
`T(t)`: Initial temperature of fluid. Alters starting thermal state.
"""
@component function FluidNode(; name, params::PackParameters, length::Float64)
    # Create fluid and thermal ports
    @named port_a = FluidPort()
    @named port_b = FluidPort()
    @named heat_port = HeatPort() 
    
    # Extract thermophysical properties
    rho = params.fluid.density
    cp  = params.fluid.specific_heat
    mu  = params.fluid.dynamic_viscosity
    k   = params.fluid.thermal_conductivity
    
    # Extract channel geometry
    W = params.tms_geometry.channel_width
    H = params.tms_geometry.channel_height
    N = params.tms_geometry.number_of_channels
    
    # Calculate hydraulic parameters
    A_flow = (W * H) * N
    P_wet = 2 * (W + H) * N
    D_h = 4 * A_flow / P_wet
    
    # Compute total fluid mass
    fluid_mass = A_flow * length * rho
    
    # Declare state variables
    sts = @variables begin
        T(t) = 298.15 
        p(t) 
        m_flow(t) 
        v(t) 
        Re(t) 
        dp(t) 
    end
    
    # Define fluid dynamics equations
    eqs = [
        # Apply mass conservation
        0 ~ port_a.m_flow + port_b.m_flow,
        m_flow ~ port_a.m_flow, 
        
        # Calculate fluid velocity and Reynolds number
        v ~ m_flow / (rho * A_flow),
        Re ~ (rho * v * D_h) / mu,
        
        # Calculate major and minor pressure losses
        dp ~ IfElse.ifelse(Re < 2300, 
                           (48.0 * mu * length / D_h^2) * v, 
                           let f_turb = (0.79 * log(max(Re, 2300.0)) - 1.64)^-2
                               f_turb * (length / D_h) * (rho * v^2 / 2.0)
                           end) + (33.0 * length) * (rho * v^2 / 2.0),
                           
        # Apply pressure drop across ports
        port_a.p - port_b.p ~ dp,
        p ~ (port_a.p + port_b.p) / 2.0, 
        
        # Equate thermal port temperatures
        heat_port.T ~ T,
        port_b.T ~ T, 
        
        # Apply convective energy balance
        fluid_mass * cp * Differential(t)(T) ~ heat_port.Q_flow + m_flow * cp * (port_a.T - T)
    ]
    
    # Return ODE system for fluid node
    ODESystem(eqs, t, sts, []; systems=[port_a, port_b, heat_port], name=name)
end

"""
    ConvectionModel(; name, params::PackParameters, length::Float64)

Acts as thermal bridge calculating dynamic convective heat flux into fluid.

Equations:
Convective heat transfer (Gnielinski correlation):
Nu = ((f/8) * (Re - 1000) * Pr) / (1 + 12.7 * √(f/8) * (Pr^(2/3) - 1))
"""
@component function ConvectionModel(; name, params::PackParameters, length::Float64)
    # Create thermal interaction ports
    @named solid_port = HeatPort()
    @named fluid_port = HeatPort()
    
    # Extract fluid properties
    k  = params.fluid.thermal_conductivity
    cp = params.fluid.specific_heat
    mu = params.fluid.dynamic_viscosity
    
    # Extract geometric dimensions
    W = params.tms_geometry.channel_width
    H = params.tms_geometry.channel_height
    N = params.tms_geometry.number_of_channels
    
    # Compute hydraulic parameters
    A_flow = (W * H) * N
    P_wet = 2 * (W + H) * N
    D_h = 4 * A_flow / P_wet
    
    # Calculate convective surface area and Prandtl number
    A_surface = P_wet * length
    Pr = (cp * mu) / k
    
    # Declare state variables
    sts = @variables begin
        Re(t) 
        Nu(t) 
        h(t) 
        Q_flow(t) 
    end
    
    # Define heat transfer equations
    eqs = [
        # Calculate Nusselt number for laminar or turbulent flow
        Nu ~ IfElse.ifelse(Re < 2300, 
                           8.23, 
                           let f_turb = (0.79 * log(max(Re, 2300.0)) - 1.64)^-2
                               ((f_turb/8.0) * (Re - 1000.0) * Pr) / (1.0 + 12.7 * sqrt(f_turb/8.0) * (Pr^(2/3) - 1.0))
                           end),
                           
        # Compute convective heat transfer coefficient
        h ~ (Nu * k) / D_h,
        
        # Calculate total convective heat flux
        Q_flow ~ h * A_surface * (solid_port.T - fluid_port.T),
        
        # Apply heat flux to boundaries
        solid_port.Q_flow ~ Q_flow,
        fluid_port.Q_flow ~ -Q_flow
    ]
    
    # Return ODE system for convection model
    ODESystem(eqs, t, sts, []; systems=[solid_port, fluid_port], name=name)
end

"""
    TMSNode(; name, params::PackParameters, length::Float64)

Container component coupling PipeWallNode, ConvectionModel, and FluidNode.
"""
@component function TMSNode(; name, params::PackParameters, length::Float64)
    # Initialise internal subcomponents
    @named wall = PipeWallNode(params=params, length=length)
    @named fluid = FluidNode(params=params, length=length)
    @named convection = ConvectionModel(params=params, length=length)
    
    # Create external routing ports
    @named port_a = FluidPort()
    @named port_b = FluidPort()
    @named heat_port = HeatPort() 
    
    # Define component connections
    eqs = [
        # Connect external fluid ports to internal fluid node
        connect(port_a, fluid.port_a),
        connect(port_b, fluid.port_b),
        
        # Connect external heat port to pipe wall
        connect(heat_port, wall.port),
        
        # Connect pipe wall and fluid node to convection bridge
        connect(wall.fluid_port, convection.solid_port),
        connect(fluid.heat_port, convection.fluid_port),
        
        # Pass Reynolds number to convection model
        convection.Re ~ fluid.Re
    ]
    
    # Return coupled ODE system for TMS node
    ODESystem(eqs, t, [], []; systems=[wall, fluid, convection, port_a, port_b, heat_port], name=name)
end