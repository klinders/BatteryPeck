# ==============================================================================
# SystemBuilder.jl
# Assembles the full LPTN by wiring cells, TMS nodes, and boundaries together
# ==============================================================================

export build_pack_system, FluidSource, FluidSink

using ModelingToolkit
using ModelingToolkitStandardLibrary.Thermal

@parameters t

"""
    FluidSource(; name, m_flow_val, T_val)

Provides constant mass flow rate and inlet temperature boundary condition.

Editable values:
`m_flow_val`: Alters total coolant mass flow entering system.
`T_val`: Alters inlet temperature of fluid.
"""
@component function FluidSource(; name, m_flow_val, T_val)
    @named port = FluidPort()
    
    # Indicate negative value for mass leaving source to enter pipe
    eqs = [
        port.m_flow ~ -m_flow_val, 
        port.T ~ T_val
    ]
    
    # Return ODE system for fluid source
    ODESystem(eqs, t, [], []; systems=[port], name=name)
end

"""
    FluidSink(; name, p_val)

Provides constant pressure boundary condition at outlet.

Editable values:
`p_val`: Adjusts absolute baseline pressure at exit node.
"""
@component function FluidSink(; name, p_val)
    @named port = FluidPort()
    
    # Set pressure at exit port
    eqs = [
        port.p ~ p_val
    ]
    
    # Return ODE system for fluid sink
    ODESystem(eqs, t, [], []; systems=[port], name=name)
end

"""
    build_pack_system(name::Symbol, geom, params::PackParameters)

Constructs complete ODESystem containing cells, casing, fluid nodes, and thermal connections.
Implements advanced casing heat rejection and intercellular gap resistances.

Editable values:
`R_radial_val`:  Transverse conduction override. Cell to cell conduction.
`R_contact_val`: Thermal resistance between cell and cooling ribbon.
"""
function build_pack_system(name::Symbol, geom, params::PackParameters)
    
    # Count total cells and TMS nodes
    num_cells = length(geom.cell_coords)
    num_tms_nodes = length(geom.tms_coords)
    
    # Calculate total cooling channel length
    total_length = 0.0
    for edge in geom.flow_edges
        t1, t2 = edge
        c1 = geom.tms_coords[t1]
        c2 = geom.tms_coords[t2]
        total_length += sqrt((c1[1] - c2[1])^2 + (c1[2] - c2[2])^2)
    end
    
    # Determine length per fluid node
    node_length = total_length / num_tms_nodes
    
    # Initialise cell components
    cells = [CoreShellCell(name=Symbol("cell_$i"), T_start=params.ambient_temperature) 
             for i in 1:num_cells]
             
    # Initialise TMS nodes
    tms_nodes = [TMSNode(name=Symbol("tms_$i"), params=params, length=node_length) 
                 for i in 1:num_tms_nodes]
                 
    # Define ambient temperature boundary
    @named ambient_temp = FixedTemperature(T = params.ambient_temperature)
    
    # Calculate true 3D surface area of pack bounding box
    L_pack = 7 * 0.025 
    W_pack = 4 * 0.025 * 0.866 
    H_pack = 0.070 
    pack_surface_area = 2 * (L_pack * W_pack + L_pack * H_pack + W_pack * H_pack)
    
    # Calculate casing volume and thermal mass
    casing_vol = pack_surface_area * params.casing_thickness
    C_casing = casing_vol * params.casing_material.density * params.casing_material.specific_heat
    @named casing_mass = HeatCapacitor(C = C_casing, T = params.ambient_temperature)
    
    # Calculate ambient convection resistance
    R_conv_val = 1.0 / (params.ambient_convection_coefficient * pack_surface_area)
    @named casing_convection = ThermalResistor(R = R_conv_val)
    
    # Initialise equation array
    eqs = Equation[]
    
    # Connect casing to ambient environment
    push!(eqs, connect(casing_mass.port, casing_convection.port_a))
    push!(eqs, connect(casing_convection.port_b, ambient_temp.port))
    
    # Retain potting bottleneck for active scenario as cells are physically embedded
    cell_diameter = 0.021
    cell_face_area = pi * (cell_diameter / 2.0)^2
    #R_axial_val = params.axial_potting_thickness / (params.potting_material.thermal_conductivity * (2 * cell_face_area))
    R_axial_val = 0.00001

    # Create axial thermal resistors
    axial_resistors = [ThermalResistor(name=Symbol("R_ax_$i"), R=R_axial_val) for i in 1:num_cells]
    
    # Connect cells to casing via axial resistors
    for i in 1:num_cells
        push!(eqs, connect(cells[i].port_shell, axial_resistors[i].port_a))
        push!(eqs, connect(axial_resistors[i].port_b, casing_mass.port))
    end
    
    # Override transverse conduction
    # SSCC is continuous 1mm thick aluminium sheet weaving through pack
    # Override silicone potting resistance to simulate highly conductive metal highway
    # Allows heat to more easily short-circuit between rows and prevents artificial downstream bottlenecking
    R_radial_val = 0.000027 
    # 0.027 old
    
    # Create intercellular gap resistors
    gap_resistors = [ThermalResistor(name=Symbol("R_gap_$idx"), R=R_radial_val) for idx in 1:length(geom.cell_edges)]
    
    # Connect adjacent cells via gap resistors
    for (idx, edge) in enumerate(geom.cell_edges)
        c1, c2 = edge
        push!(eqs, connect(cells[c1].port_shell, gap_resistors[idx].port_a))
        push!(eqs, connect(gap_resistors[idx].port_b, cells[c2].port_shell))
    end
    
    # Tunable contact resistance between cell and cooling ribbon
    R_contact_val = 3.5
    # 3.5 old
    
    # Create contact resistors
    contact_resistors = [ThermalResistor(name=Symbol("R_contact_$idx"), R=R_contact_val) 
                         for idx in 1:length(geom.convection_edges)]
    
    # Connect cells to TMS nodes
    for (idx, edge) in enumerate(geom.convection_edges)
        c_idx, t_idx = edge
        push!(eqs, connect(cells[c_idx].port_shell, contact_resistors[idx].port_a))
        push!(eqs, connect(contact_resistors[idx].port_b, tms_nodes[t_idx].heat_port))
    end
    
    # Connect fluid nodes in series
    for edge in geom.flow_edges
        t1, t2 = edge
        push!(eqs, connect(tms_nodes[t1].port_b, tms_nodes[t2].port_a))
    end
    
    # Identify inlet and outlet nodes
    first_tms = tms_nodes[geom.flow_edges[1][1]]
    last_tms  = tms_nodes[geom.flow_edges[end][2]]
    
    # Define fluid inlet and outlet boundaries
    @named fluid_inlet = FluidSource(m_flow_val = params.mass_flow_rate, T_val = params.inlet_temperature)
    @named fluid_outlet = FluidSink(p_val = 101325.0)
    
    # Connect boundaries to fluid network
    push!(eqs, connect(fluid_inlet.port, first_tms.port_a))
    push!(eqs, connect(last_tms.port_b, fluid_outlet.port))
    
    # Aggregate all system components
    systems = Any[ambient_temp, casing_mass, casing_convection, fluid_inlet, fluid_outlet]
    append!(systems, cells)
    append!(systems, tms_nodes)
    append!(systems, axial_resistors)
    append!(systems, gap_resistors)
    append!(systems, contact_resistors)

    # Return complete ODE system
    return ODESystem(eqs, t, [], []; systems=systems, name=name)
end