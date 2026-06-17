# ==============================================================================
# SystemBuilder.jl
# Assembles the full LPTN by wiring cells, TMS nodes, and boundaries together
# Contains:
# 1. FluidSource: Provide dynamic mass flow rate and inlet temperature boundary condition
# 2. FluidSink: Provide constant pressure boundary condition at outlet
# 3. build_pack_system: Construct complete ODESystem containing cells, casing, fluid nodes, and thermal connections
# ==============================================================================

export build_pack_system, FluidSource, FluidSink

using ModelingToolkit
using ModelingToolkitStandardLibrary.Thermal

@parameters t

"""
    FluidSource(; name, m_flow_val, T_val)

Provide dynamic mass flow rate and inlet temperature boundary condition.

# Arguments
- `name`: Component name for ModelingToolkit
- `m_flow_val`: Initial mass flow rate constraint
- `T_val`: Initial inlet temperature constraint

# Returns
- Instantiated FluidSource component system
"""
@component function FluidSource(; name, m_flow_val, T_val)
    @named port = FluidPort()
    
    # Define boundaries as parameters to allow callback intervention
    @parameters begin
        m_flow_in = m_flow_val
        T_inlet = T_val
    end
    
    # Apply mass flow and temperature port constraints
    eqs = [
        port.m_flow ~ -m_flow_in, 
        port.T ~ T_inlet
    ]
    
    return ODESystem(eqs, t, [], [m_flow_in, T_inlet]; systems=[port], name=name)
end

"""
    FluidSink(; name, p_val)

Provide constant pressure boundary condition at outlet.

Editable values:
`p_val`: Adjusts absolute baseline pressure at exit node.

# Arguments
- `name`: Component name for ModelingToolkit
- `p_val`: Constant exit pressure value

# Returns
- Instantiated FluidSink component system
"""
@component function FluidSink(; name, p_val)
    @named port = FluidPort()
    
    # Set pressure at exit port
    eqs = [
        port.p ~ p_val
    ]
    
    return ODESystem(eqs, t, [], []; systems=[port], name=name)
end

"""
    build_pack_system(name::Symbol, geom, params::PackParameters)

Construct complete ODESystem containing cells, casing, fluid nodes, and thermal connections.

Implements advanced casing heat rejection and intercellular gap resistances.

Editable values:
`R_radial_val`:  Transverse conduction override. Cell to cell conduction.
`R_contact_val`: Thermal resistance between cell and cooling ribbon.

# Arguments
- `name::Symbol`: Base name for constructed system
- `geom`: Pack geometry structure containing coordinates and edges
- `params::PackParameters`: Struct containing physical and thermal properties

# Returns
- Complete coupled ODE system representing battery pack
"""
function build_pack_system(name::Symbol, geom, params::PackParameters)
    
    # Count elements and calculate total cooling channel length to determine individual node lengths
    num_cells = length(geom.cell_coords)
    num_tms_nodes = length(geom.tms_coords)
    
    total_length = 0.0
    for edge in geom.flow_edges
        t1, t2 = edge
        c1 = geom.tms_coords[t1]
        c2 = geom.tms_coords[t2]
        total_length += sqrt((c1[1] - c2[1])^2 + (c1[2] - c2[2])^2)
    end
    
    node_length = total_length / num_tms_nodes
    
    # Initialise cell components, TMS nodes, and ambient temperature boundary
    cells = [CoreShellCell(name=Symbol("cell_$i"), T_start=params.ambient_temperature) 
             for i in 1:num_cells]
             
    tms_nodes = [TMSNode(name=Symbol("tms_$i"), params=params, length=node_length) 
                 for i in 1:num_tms_nodes]
                 
    @named ambient_temp = FixedTemperature(T = params.ambient_temperature)
    
    # Calculate true 3D surface area of pack bounding box to derive casing thermal mass
    L_pack = 7 * 0.025 
    W_pack = 4 * 0.025 * 0.866 
    H_pack = 0.070 
    pack_surface_area = 2 * (L_pack * W_pack + L_pack * H_pack + W_pack * H_pack)
    
    casing_vol = pack_surface_area * params.casing_thickness
    C_casing = casing_vol * params.casing_material.density * params.casing_material.specific_heat
    @named casing_mass = HeatCapacitor(C = C_casing, T = params.ambient_temperature)
    
    # Compute ambient convection resistance and connect casing to environment
    R_conv_val = 1.0 / (params.ambient_convection_coefficient * pack_surface_area)
    @named casing_convection = ThermalResistor(R = R_conv_val)
    
    eqs = Equation[]
    
    push!(eqs, connect(casing_mass.port, casing_convection.port_a))
    push!(eqs, connect(casing_convection.port_b, ambient_temp.port))
    
    # Retain potting bottleneck for active scenario as cells are physically embedded
    cell_diameter = 0.021
    cell_face_area = pi * (cell_diameter / 2.0)^2
    # R_axial_val = params.axial_potting_thickness / (params.potting_material.thermal_conductivity * (2 * cell_face_area))
    R_axial_val = 0.25

    # Connect cells to casing via axial thermal resistors
    axial_resistors = [ThermalResistor(name=Symbol("R_ax_$i"), R=R_axial_val) for i in 1:num_cells]
    
    for i in 1:num_cells
        push!(eqs, connect(cells[i].port_shell, axial_resistors[i].port_a))
        push!(eqs, connect(axial_resistors[i].port_b, casing_mass.port))
    end
    
    # Override transverse conduction
    # SSCC is continuous 1mm thick aluminium sheet weaving through pack
    # Override silicone potting resistance to simulate highly conductive metal highway
    # Allows heat to more easily short-circuit between rows and prevents artificial downstream bottlenecking
    R_radial_val = 0.2
    # 0.027 old
    
    # Connect adjacent cells via gap resistors
    gap_resistors = [ThermalResistor(name=Symbol("R_gap_$idx"), R=R_radial_val) for idx in 1:length(geom.cell_edges)]
    
    for (idx, edge) in enumerate(geom.cell_edges)
        c1, c2 = edge
        push!(eqs, connect(cells[c1].port_shell, gap_resistors[idx].port_a))
        push!(eqs, connect(gap_resistors[idx].port_b, cells[c2].port_shell))
    end
    
    # Tunable contact resistance between cell and cooling ribbon
    R_contact_val = 4
    # 3.5 old
    
    # Connect cells to TMS nodes via contact resistors
    contact_resistors = [ThermalResistor(name=Symbol("R_contact_$idx"), R=R_contact_val) 
                         for idx in 1:length(geom.convection_edges)]
    
    for (idx, edge) in enumerate(geom.convection_edges)
        c_idx, t_idx = edge
        push!(eqs, connect(cells[c_idx].port_shell, contact_resistors[idx].port_a))
        push!(eqs, connect(contact_resistors[idx].port_b, tms_nodes[t_idx].heat_port))
    end
    
    # Route fluid nodes in series and define inlet and outlet boundary conditions
    for edge in geom.flow_edges
        t1, t2 = edge
        push!(eqs, connect(tms_nodes[t1].port_b, tms_nodes[t2].port_a))
    end
    
    first_tms = tms_nodes[geom.flow_edges[1][1]]
    last_tms  = tms_nodes[geom.flow_edges[end][2]]
    
    @named fluid_inlet = FluidSource(m_flow_val = params.mass_flow_rate, T_val = params.inlet_temperature)
    @named fluid_outlet = FluidSink(p_val = 101325.0)
    
    push!(eqs, connect(fluid_inlet.port, first_tms.port_a))
    push!(eqs, connect(last_tms.port_b, fluid_outlet.port))
    
    # Aggregate all system components and return ODE system
    systems = Any[ambient_temp, casing_mass, casing_convection, fluid_inlet, fluid_outlet]
    append!(systems, cells)
    append!(systems, tms_nodes)
    append!(systems, axial_resistors)
    append!(systems, gap_resistors)
    append!(systems, contact_resistors)

    return ODESystem(eqs, t, [], []; systems=systems, name=name)
end