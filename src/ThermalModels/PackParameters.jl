# ==============================================================================
# PackParameters.jl
# Define thermophysical properties, boundary conditions, and geometric constants
# ==============================================================================

export FluidProperties, SolidProperties, TMSGeometry, PackParameters
export get_coolant_properties, get_solid_properties

"""
    FluidProperties

Contains constant thermophysical properties for single-phase coolant fluid.
"""
Base.@kwdef struct FluidProperties
    # Define fluid density
    density::Float64
    # Define specific heat capacity
    specific_heat::Float64
    # Define thermal conductivity
    thermal_conductivity::Float64
    # Define dynamic viscosity
    dynamic_viscosity::Float64
end

"""
    SolidProperties

Contains constant thermophysical properties for solid pack materials.
"""
Base.@kwdef struct SolidProperties
    # Define solid density
    density::Float64
    # Define specific heat capacity
    specific_heat::Float64
    # Define thermal conductivity
    thermal_conductivity::Float64
end

"""
    TMSGeometry

Contains cross-sectional dimensions of cooling ribbon for hydraulic calculations.
"""
Base.@kwdef struct TMSGeometry
    # Set internal channel width
    channel_width::Float64
    # Set internal channel height
    channel_height::Float64
    # Define quantity of parallel flow channels
    number_of_channels::Int
    # Set thickness of enclosing pipe wall
    wall_thickness::Float64
end

"""
    PackParameters

Master dictionary containing environmental boundaries, geometries, and materials.
"""
Base.@kwdef struct PackParameters
    # Assign fluid properties
    fluid::FluidProperties
    # Assign pipe wall properties
    pipe_wall::SolidProperties
    # Assign thermal potting properties
    potting_material::SolidProperties
    # Assign outer casing properties
    casing_material::SolidProperties
    # Assign cooling channel dimensions
    tms_geometry::TMSGeometry
    
    # Set physical gap between adjacent cells
    cell_gap_thickness::Float64
    # Set thickness of potting material at cell base
    axial_potting_thickness::Float64
    # Set thickness of outer pack casing
    casing_thickness::Float64
    
    # Set ambient environmental temperature
    ambient_temperature::Float64
    # Set coolant inlet temperature
    inlet_temperature::Float64
    # Set total coolant mass flow rate
    mass_flow_rate::Float64
    # Set convective heat transfer coefficient for outer casing
    ambient_convection_coefficient::Float64
end

"""
    get_coolant_properties(coolant_type::Symbol)

Returns baseline fluid properties evaluated at nominal 25°C.
"""
function get_coolant_properties(coolant_type::Symbol)
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

Returns baseline thermophysical properties for standard pack materials.
"""
function get_solid_properties(material_type::Symbol)
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
    build_pack_parameters(; kwargs...)

Convenience constructor for generating full parameter set with standard defaults.

Editable values:
`coolant`: Alters fluid properties dictionary.
`ambient_temp`: Sets environmental baseline temperature.
`inlet_temp`: Sets fluid entry temperature.
`flow_rate`: Adjusts global mass flow rate.
`cell_pitch`: Changes physical spacing between cell centres.
`cell_diameter`: Adjusts active cell diameter for gap calculations.
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