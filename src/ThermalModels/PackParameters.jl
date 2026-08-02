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
    get_coolant_properties(coolant_type)

Return baseline fluid properties evaluated at nominal 25°C.

Evaluate and return properties based on specified coolant type. Throw error if unsupported.

# Arguments
- `coolant_type`: Symbol identifier for desired coolant fluid

# Returns
- Instantiated FluidProperties struct
"""
function get_coolant_properties(coolant_type::Symbol)
    # Evaluate properties based on specified coolant type using literature reference data
    if coolant_type == :water_glycol
        return FluidProperties(
            density = 1071.11,
            specific_heat = 3300.0,
            thermal_conductivity = 0.384,
            dynamic_viscosity = 0.00257
        )
    elseif coolant_type == :air
        return FluidProperties(
            density = 1.1614,
            specific_heat = 1007.0,
            thermal_conductivity = 0.0263,
            dynamic_viscosity = 1.846e-5
        )
    elseif coolant_type == :water_nema2026
        return FluidProperties(
            density = 998.0,
            specific_heat = 4180.0,
            thermal_conductivity = 0.61,
            dynamic_viscosity = 0.002 
        )
    else
        error("Unknown coolant type")
    end
end

"""
    get_solid_properties(material_type)

Return baseline thermophysical properties for standard pack materials.

Evaluate and return properties based on specified solid material. Throw error if unsupported.

# Arguments
- `material_type`: Symbol identifier for desired solid material

# Returns
- Instantiated SolidProperties struct
"""
function get_solid_properties(material_type::Symbol)
    # Evaluate properties based on specified solid material using material data sheets
    if material_type == :aluminium_6061
        return SolidProperties(
            density = 2700.0,
            specific_heat = 895.0,
            thermal_conductivity = 166.0
        )
    elseif material_type == :bergquist_tgf_1500
        return SolidProperties(
            density = 2700.0,
            specific_heat = 1000.0,
            thermal_conductivity = 1.8
        )
    elseif material_type == :aluminium_nema2026
        return SolidProperties(
            density = 2700.0,
            specific_heat = 900.0,
            thermal_conductivity = 238.0
        )
    else
        error("Unknown material type")
    end
end

"""
    build_pack_parameters(; coolant, ambient_temp, inlet_temp, flow_rate, cell_pitch, cell_diameter, h_conv, casing_th)

Generate full parameter set with standard defaults via convenience constructor.

Construct and return master parameter set using provided parameters or defaults.

# Arguments
- `coolant`: Symbol identifier for desired coolant fluid
- `ambient_temp`: Float64 environmental baseline temperature
- `inlet_temp`: Float64 fluid entry temperature
- `flow_rate`: Float64 global mass flow rate
- `cell_pitch`: Float64 physical spacing between cell centres
- `cell_diameter`: Float64 active cell diameter
- `h_conv`: Float64 ambient convection coefficient
- `casing_th`: Float64 casing thickness

# Returns
- Instantiated PackParameters master dictionary
"""
function build_pack_parameters(;
        coolant::Symbol = :water_nema2026,
        ambient_temp::Float64 = 298.15,
        inlet_temp::Float64 = 298.15,
        flow_rate::Float64 = 0.02994,
        cell_pitch::Float64 = 0.025,
        cell_diameter::Float64 = 0.021,
        h_conv::Float64 = 5.0,
        casing_th::Float64 = 0.01
    )
    
    # Define baseline cooling channel geometry matching validation data
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
        casing_thickness = casing_th,        
        
        ambient_temperature = ambient_temp,
        inlet_temperature = inlet_temp,
        mass_flow_rate = flow_rate,
        ambient_convection_coefficient = h_conv 
    )
end