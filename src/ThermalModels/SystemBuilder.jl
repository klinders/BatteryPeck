# ==============================================================================
# SystemBuilder.jl
# Assembles the full LPTN by wiring Cells, TMS Nodes, and Boundaries together
# ==============================================================================

export build_pack_system

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
    eqs = [
        # Negative value indicates mass leaves source to enter pipe
        port.m_flow ~ -m_flow_val, 
        port.T ~ T_val
    ]
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
    eqs = [
        port.p ~ p_val
    ]
    ODESystem(eqs, t, [], []; systems=[port], name=name)
end

"""
    build_pack_system(name::Symbol, geom, params::PackParameters)

Constructs complete ODESystem containing cells, casing, fluid nodes, and thermal connections.
Implements advanced casing heat rejection and intercellular gap resistances.

Editable values:
`R_contact_val`: Thermal resistance between cell and cooling ribbon. Adjusts right side of temperature curve.
"""
function build_pack_system(name::Symbol, geom, params::PackParameters)
    
    # Extracts counts from geometry arrays safely
    num_cells = length(geom.cell_coords)
    num_tms_nodes = length(geom.tms_coords)
    
    # Calculates physical length of cooling pipe
    total_length = 0.0
    for edge in geom.flow_edges
        t1, t2 = edge
        c1 = geom.tms_coords[t1]
        c2 = geom.tms_coords[t2]
        total_length += sqrt((c1[1] - c2[1])^2 + (c1[2] - c2[2])^2)
    end
    node_length = total_length / num_tms_nodes
    
    # Instantiates core components
    cells = [CoreShellCell(name=Symbol("cell_$i"), T_start=params.ambient_temperature) 
             for i in 1:num_cells]
             
    tms_nodes = [TMSNode(name=Symbol("tms_$i"), params=params, length=node_length) 
                 for i in 1:num_tms_nodes]
                 
    # Instantiates global boundaries and casing
    @named ambient_temp = FixedTemperature(T = params.ambient_temperature)
    
    casing_area = num_cells * (0.021 * 0.021) 
    casing_vol = casing_area * params.casing_thickness
    C_casing = casing_vol * params.casing_material.density * params.casing_material.specific_heat
    @named casing_mass = HeatCapacitor(C = C_casing, T = params.ambient_temperature)
    
    # Calculates convective resistance
    R_conv_val = 1.0 / (params.ambient_convection_coefficient * casing_area)
    @named casing_convection = ThermalResistor(R = R_conv_val)
    
    # Establishes connections
    eqs = Equation[]
    
    push!(eqs, connect(casing_mass.port, casing_convection.port_a))
    push!(eqs, connect(casing_convection.port_b, ambient_temp.port))
    
    # Multiplies cell face area by two because heat leaves both top and bottom faces to casing
    cell_face_area = pi * (params.tms_geometry.channel_width / 2.0)^2
    R_axial_val = params.axial_potting_thickness / (params.potting_material.thermal_conductivity * (2 * cell_face_area))
    axial_resistors = [ThermalResistor(name=Symbol("R_ax_$i"), R=R_axial_val) for i in 1:num_cells]
    
    for i in 1:num_cells
        push!(eqs, connect(cells[i].port_shell, axial_resistors[i].port_a))
        push!(eqs, connect(axial_resistors[i].port_b, casing_mass.port))
    end
    
    R_radial_val = params.cell_gap_thickness / (params.potting_material.thermal_conductivity * (0.065 * 0.021)) 
    gap_resistors = [ThermalResistor(name=Symbol("R_gap_$idx"), R=R_radial_val) for idx in 1:length(geom.cell_edges)]
    
    for (idx, edge) in enumerate(geom.cell_edges)
        c1, c2 = edge
        push!(eqs, connect(cells[c1].port_shell, gap_resistors[idx].port_a))
        push!(eqs, connect(gap_resistors[idx].port_b, cells[c2].port_shell))
    end
    
    # Represents tunable contact resistance between cell and cooling ribbon
    R_contact_val = 2.7 
    
    contact_resistors = [ThermalResistor(name=Symbol("R_contact_$idx"), R=R_contact_val) 
                         for idx in 1:length(geom.convection_edges)]
    
    for (idx, edge) in enumerate(geom.convection_edges)
        c_idx, t_idx = edge
        push!(eqs, connect(cells[c_idx].port_shell, contact_resistors[idx].port_a))
        push!(eqs, connect(contact_resistors[idx].port_b, tms_nodes[t_idx].heat_port))
    end
    
    for edge in geom.flow_edges
        t1, t2 = edge
        push!(eqs, connect(tms_nodes[t1].port_b, tms_nodes[t2].port_a))
    end
    
    # Applies fluid boundary conditions
    first_tms = tms_nodes[geom.flow_edges[1][1]]
    last_tms  = tms_nodes[geom.flow_edges[end][2]]
    
    @named fluid_inlet = FluidSource(m_flow_val = params.mass_flow_rate, T_val = params.inlet_temperature)
    @named fluid_outlet = FluidSink(p_val = 101325.0)
    
    push!(eqs, connect(fluid_inlet.port, first_tms.port_a))
    push!(eqs, connect(last_tms.port_b, fluid_outlet.port))
    
    # Assembles final system
    systems = Any[ambient_temp, casing_mass, casing_convection, fluid_inlet, fluid_outlet]
    append!(systems, cells)
    append!(systems, tms_nodes)
    append!(systems, axial_resistors)
    append!(systems, gap_resistors)
    append!(systems, contact_resistors)

    return ODESystem(eqs, t, [], []; systems=systems, name=name)
end