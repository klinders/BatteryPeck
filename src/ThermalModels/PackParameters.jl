# ==============================================================================
# PackParameters.jl
# Define thermophysical properties, boundary conditions, and geometric constants
# Contains:
# 1. FluidProperties: Contain constant thermophysical properties for single-phase coolant fluid
# 2. SolidProperties: Contain constant thermophysical properties for solid pack materials
# 3. TMSGeometry: Contain cross-sectional dimensions of cooling ribbon for hydraulic calculations
# 4. PackParameters: Store environmental boundaries, geometries, and materials in master dictionary
# 5. get_coolant_properties: Return baseline fluid properties evaluated at nominal 25°C
# 6. get_solid_properties: Return baseline thermophysical properties for standard pack materials
# 7. build_pack_parameters: Generate full parameter set with standard defaults via convenience constructor
# ==============================================================================

export FluidProperties, SolidProperties, TMSGeometry, PackParameters
export get_coolant_properties, get_solid_properties

"""
    FluidProperties

Contain constant thermophysical properties for single-phase coolant fluid.
"""
Base.@kwdef struct FluidProperties
    # Define fluid thermophysical properties
    density::Float64
    specific_heat::Float64
    thermal_conductivity::Float64
    dynamic_viscosity::Float64
end

"""
    SolidProperties

Contain constant thermophysical properties for solid pack materials.
"""
Base.@kwdef struct SolidProperties
    # Define solid thermophysical properties
    density::Float64
    specific_heat::Float64
    thermal_conductivity::Float64
end

"""
    TMSGeometry

Contain cross-sectional dimensions of cooling ribbon for hydraulic calculations.
"""
Base.@kwdef struct TMSGeometry
    # Define internal channel dimensions and pipe wall thickness
    channel_width::Float64
    channel_height::Float64
    number_of_channels::Int
    wall_thickness::Float64
end

"""
    PackParameters

Store environmental boundaries, geometries, and materials in master dictionary.
"""
Base.@kwdef struct PackParameters
    # Assign material properties and channel dimensions
    fluid::FluidProperties
    pipe_wall::SolidProperties
    potting_material::SolidProperties
    casing_material::SolidProperties
    tms_geometry::TMSGeometry
    
    # Set physical gap and casing thicknesses
    cell_gap_thickness::Float64
    axial_potting_thickness::Float64
    casing_thickness::Float64
    
    # Set environmental temperatures and flow boundaries
    ambient_temperature::Float64
    inlet_temperature::Float64
    mass_flow_rate::Float64
    ambient_convection_coefficient::Float64
end

"""
    get_coolant_properties(coolant_type::Symbol)

Return baseline fluid properties evaluated at nominal 25°C.

# Arguments
- `coolant_type::Symbol`: Identifier for desired coolant fluid

# Returns
- Instantiated FluidProperties struct
"""
function get_coolant_properties(coolant_type::Symbol)
    # Evaluate and return properties based on specified coolant type
    if coolant_type == :water_glycol
        # 50/50 volume ethylene glycol and water mixture at 25°C
        # Source: 2001 ASHRAE Fundamentals Handbook (SI), Chapter 21
        return FluidProperties(
            density = 1071.11,
            specific_heat = 3300.0,
            thermal_conductivity = 0.384,
            dynamic_viscosity = 0.00257
        )
    elseif coolant_type == :air
        # Standard atmospheric air at 300 K
        # Source: Fundamentals of Heat and Mass Transfer, Appendix A.4
        return FluidProperties(
            density = 1.1614,
            specific_heat = 1007.0,
            thermal_conductivity = 0.0263,
            dynamic_viscosity = 1.846e-5
        )
    elseif coolant_type == :water_nema2026
        # Pure water properties specifically hardcoded in Nema et al. (2026) Table 1
        return FluidProperties(
            density = 998.0,
            specific_heat = 4180.0,
            thermal_conductivity = 0.61,
            dynamic_viscosity = 0.002 
        )
    # Handle unsupported coolant types
    else
        error("Unknown coolant type. Choose :water_glycol or :air.")
    end
end

"""
    get_solid_properties(material_type::Symbol)

Return baseline thermophysical properties for standard pack materials.

# Arguments
- `material_type::Symbol`: Identifier for desired solid material

# Returns
- Instantiated SolidProperties struct
"""
function get_solid_properties(material_type::Symbol)
    # Evaluate and return properties based on specified solid material
    if material_type == :aluminium_6061
        # Alloy 6061 properties
        # Source: Bohler Uddeholm Aluminium 6061 Data Sheet
        return SolidProperties(
            density = 2700.0,
            specific_heat = 895.0,
            thermal_conductivity = 166.0
        )
    elseif material_type == :bergquist_tgf_1500
        # Thermally conductive silicone liquid gap filler
        # Source: Bergquist Gap Filler TGF 1500 Technical Data Sheet
        return SolidProperties(
            density = 2700.0,
            specific_heat = 1000.0,
            thermal_conductivity = 1.8
        )
    elseif material_type == :aluminium_nema2026
        # Pure aluminium properties from Nema et al. (2026) Table 1
        return SolidProperties(
            density = 2700.0,
            specific_heat = 900.0,
            thermal_conductivity = 238.0
        )
    # Handle unsupported material types
    else
        error("Unknown material type. Choose :aluminium_6061 or :bergquist_tgf_1500.")
    end
end

"""
    build_pack_parameters(; coolant=:water_nema2026, ambient_temp=298.15, inlet_temp=298.15, flow_rate=0.02994, cell_pitch=0.025, cell_diameter=0.021)

Generate full parameter set with standard defaults via convenience constructor.

Editable values:
`coolant`: Alters fluid properties dictionary.
`ambient_temp`: Sets environmental baseline temperature.
`inlet_temp`: Sets fluid entry temperature.
`flow_rate`: Adjusts global mass flow rate.
`cell_pitch`: Changes physical spacing between cell centres.
`cell_diameter`: Adjusts active cell diameter for gap calculations.

# Arguments
- `coolant::Symbol`: Identifier for desired coolant fluid
- `ambient_temp::Float64`: Environmental baseline temperature
- `inlet_temp::Float64`: Fluid entry temperature
- `flow_rate::Float64`: Global mass flow rate
- `cell_pitch::Float64`: Physical spacing between cell centres
- `cell_diameter::Float64`: Active cell diameter

# Returns
- Instantiated PackParameters master dictionary
"""
function build_pack_parameters(;
        coolant::Symbol = :water_nema2026,
        ambient_temp::Float64 = 298.15,
        inlet_temp::Float64 = 298.15,
        flow_rate::Float64 = 0.02994,
        cell_pitch::Float64 = 0.025,
        cell_diameter::Float64 = 0.021
    )
    
    # Define baseline cooling channel geometry matching Nema 2026 1C/5C validation
    tms_geom = TMSGeometry(
        channel_width = 0.002,
        channel_height = 0.050,
        number_of_channels = 1,
        wall_thickness = 0.001
    )
    
    # Construct and return master parameter set
    return PackParameters(
        fluid = get_coolant_properties(coolant),
        pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500),
        casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom,
        
        cell_gap_thickness = (cell_pitch - cell_diameter) / 2.0,
        axial_potting_thickness = 0.005, 
        casing_thickness = 0.01,        
        
        ambient_temperature = ambient_temp,
        inlet_temperature = inlet_temp,
        mass_flow_rate = flow_rate,
        ambient_convection_coefficient = 5.0 
    )
end