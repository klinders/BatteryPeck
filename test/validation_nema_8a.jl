# ==============================================================================
# validation_nema_8a.jl     (passive cooling)
# Validates 1D/2D LPTN thermal network against Nema et al. (2026) 3D CFD data Fig 8a
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq 
using DataInterpolations
using DelimitedFiles
using Plots
using Logging

using Revise
using BatteryToolkit

Revise.revise()

"""
    load_and_preprocess(filepath)

Loads CSV data and filters out backward-pass digitisation artifacts.
Ensures strictly increasing time array for interpolation stability.
"""
function load_and_preprocess(filepath)
    # Read raw data from specified filepath
    raw_data = readdlm(filepath, ',')
    
    # Remove header row if present
    if typeof(raw_data[1,1]) <: AbstractString
        raw_data = raw_data[2:end, :]
    end
    
    # Cast array to float format
    raw_data = Float64.(raw_data)
    
    # Extract time array
    t_raw = raw_data[:, 1]
    # Extract dependent variable array
    y_raw = raw_data[:, 2]
    
    # Find index of maximum time value
    max_idx = argmax(t_raw)
    
    # Truncate time array to remove trailing artifacts
    t_clean = t_raw[1:max_idx]
    # Truncate value array accordingly
    y_clean = y_raw[1:max_idx]
    
    # Initialise array of valid indices
    valid_idx = [1]
    
    # Loop through truncated array to enforce monotonic increase
    for i in 2:length(t_clean)
        # Compare current time to last valid time
        if t_clean[i] > t_clean[valid_idx[end]]
            # Append index if strictly increasing
            push!(valid_idx, i)
        end
    end
    
    # Apply monotonic filter to time array
    t_final = t_clean[valid_idx]
    # Apply monotonic filter to value array
    y_final = y_clean[valid_idx]
    
    # Shift time array to start at zero
    t_shifted = t_final .- minimum(t_final)
    
    # Return cleaned arrays
    return t_shifted, y_final
end

"""
    compute_metrics(t_sim, y_sim, t_val, y_val)

Calculates RMSE (Root Mean Square Error) and R² (Coefficient of Determination).

Equations:
RMSE = sqrt( sum((y_sim - y_val)^2) / N )
R² = 1 - ( sum((y_val - y_sim)^2) / sum((y_val - y_mean)^2) )
"""
function compute_metrics(t_sim, y_sim, t_val, y_val)
    # Create linear interpolation of simulation data
    sim_interp = LinearInterpolation(y_sim, t_sim)
    
    # Find indices where validation time is within simulation bounds
    valid_idx = findall(t -> t <= t_sim[end], t_val)
    
    # Trim validation time array
    t_val_trimmed = t_val[valid_idx]
    # Trim validation value array
    y_val_trimmed = y_val[valid_idx]
    
    # Evaluate simulation interpolation at validation time points
    y_sim_at_val = [sim_interp(t) for t in t_val_trimmed]
    
    # Root mean square error
    rmse = sqrt(sum((y_sim_at_val .- y_val_trimmed).^2) / length(y_val_trimmed))
    
    # Mean of validation values
    y_mean = sum(y_val_trimmed) / length(y_val_trimmed)
    
    # Total sum of squares
    ss_tot = sum((y_val_trimmed .- y_mean).^2)
    # Residual sum of squares
    ss_res = sum((y_val_trimmed .- y_sim_at_val).^2)
    
    # Coefficient of determination
    r2 = 1.0 - (ss_res / ss_tot)
    
    # Return calculated error metrics
    return rmse, r2
end

# Define path to validation data directory
data_dir = joinpath(@__DIR__, "..", "data", "Nema2026")

# Load and clean 1C heat generation data
t_hg_1c, y_hg_1c = load_and_preprocess(joinpath(data_dir, "7a_298K_HG.csv"))
# Create linear interpolator for 1C data
const hg_interp_1c = LinearInterpolation(y_hg_1c, t_hg_1c)

"""
    q_gen_1c(t)

Returns volumetric heat generation for 1C discharge rate.
Interpolates Nema et al. (2026) Figure 7a data.
"""
function q_gen_1c(t)
    # Return initial value if time is before data start
    if t < t_hg_1c[1] return y_hg_1c[1]
    # Return final value if time exceeds data end
    elseif t > t_hg_1c[end] return y_hg_1c[end]
    # Return interpolated value within data bounds
    else return hg_interp_1c(t) end
end

# Register 1C function for symbolic execution
@register_symbolic q_gen_1c(t)

# Load and clean 5C heat generation data
t_hg_5c, y_hg_5c = load_and_preprocess(joinpath(data_dir, "7b_298K_HG.csv"))
# Create linear interpolator for 5C data
const hg_interp_5c = LinearInterpolation(y_hg_5c, t_hg_5c)

"""
    q_gen_5c(t)

Returns volumetric heat generation for 5C discharge rate.
Interpolates Nema et al. (2026) Figure 7b data.
"""
function q_gen_5c(t)
    # Return initial value if time is before data start
    if t < t_hg_5c[1] return y_hg_5c[1]
    # Return final value if time exceeds data end
    elseif t > t_hg_5c[end] return y_hg_5c[end]
    # Return interpolated value within data bounds
    else return hg_interp_5c(t) end
end

# Register 5C function for symbolic execution
@register_symbolic q_gen_5c(t)

# Set near zero coolant velocity to simulate passive cooling
v_coolant = 0.0000001
# Define fluid density
rho_water = 998.0
# Define channel width
W_channel = 0.004 
# Define channel height
H_channel = 0.050 

# Calculate cross sectional area
A_cross = W_channel * H_channel
# Calculate baseline mass flow rate
m_dot_baseline = rho_water * A_cross * v_coolant 

# Initialise cooling channel geometry
tms_geom_nema = TMSGeometry(
    channel_width = W_channel, 
    channel_height = H_channel, 
    number_of_channels = 1, 
    wall_thickness = 0.001  
)

# Initialise pack parameters with materials and boundary conditions
val_params = PackParameters(
    fluid = get_coolant_properties(:water_nema2026),
    pipe_wall = get_solid_properties(:aluminium_nema2026),
    potting_material = get_solid_properties(:bergquist_tgf_1500),
    casing_material = get_solid_properties(:aluminium_nema2026),
    tms_geometry = tms_geom_nema,
    
    cell_gap_thickness = (0.025 - 0.021) / 2.0, 
    axial_potting_thickness = 0.005,
    casing_thickness = 0.003,
    
    ambient_temperature = 298.15,
    inlet_temperature = 298.15,
    mass_flow_rate = m_dot_baseline, 
    ambient_convection_coefficient = 16.0
)

# Build spatial geometry for battery pack
geom = build_pack_geometry(
    rows = 4, cols = 7, pattern = :hexagonal, 
    tms_routing = :single_row, tms_encasement = :full, cell_pitch = 0.025
)

# Count total number of cells
num_cells = length(geom.cell_coords)

# Declare independent time variable
@parameters t

"""
    run_and_validate(rate_name, q_func, t_hg_end, val_filename)

Builds LPTN system, runs solver, calculates metrics, and plots validation results.
"""
function run_and_validate(rate_name, q_func, t_hg_end, val_filename)
    # Print validation start message to console
    println("\nStarting validation for: $rate_name discharge")
    
    # Construct and simplifies ODE system while suppressing console warnings
    sys_simplified, cells = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        # Build base thermal system
        sys_base = build_pack_system(Symbol("sys_$rate_name"), geom, val_params)
        
        # Extract cell components from system
        cls = [getproperty(sys_base, Symbol("cell_$i")) for i in 1:num_cells]
        
        # Apply specific heat generation function to each cell
        forcing_eqs = [cls[i].Q_volumetric_in.u ~ q_func(t) for i in 1:num_cells]
        
        # Couple base system with forcing equations
        full_sys = ODESystem(forcing_eqs, t; systems=[sys_base], name=Symbol("full_$rate_name"))
        
        # Return structurally simplified system and cell array
        (structural_simplify(full_sys), cls) 
    end

    # Define ODE problem with jacobian and sparse matrix forms
    prob = ODEProblem(sys_simplified, [], (0.0, t_hg_end), jac=true, sparse=true)
    
    # Solve ODE problem and times execution
    @time sol = solve(prob, QNDF(), saveat=1.0, dtmax=1.0) 

    # Extract simulation time array
    t_sim = sol.t
    
    # Extract cell temperatures in Kelvin
    T_cells_K = [sol[cells[i].core_cap.T] for i in 1:num_cells]
    
    # Find maximum cell temperature at each time step
    T_max_sim_K = maximum(hcat(T_cells_K...), dims=2)[:, 1]
    
    # Convert maximum temperature to Celsius
    T_max_sim_C = T_max_sim_K .- 273.15

    # Load validation data from specified filename
    t_val, y_val_K = load_and_preprocess(joinpath(data_dir, val_filename))
    
    # Convert validation temperature to Celsius
    y_val_C = y_val_K .- 273.15
    
    # Compute error metrics between simulation and validation data
    rmse, r2 = compute_metrics(t_sim, T_max_sim_C, t_val, y_val_C)
    
    # Print calculated RMSE to console
    println("  -> RMSE: $(round(rmse, digits=3)) °C")
    # Print calculated R² to console
    println("  -> R²:   $(round(r2, digits=4))")

    # Initialise plot with simulation data
    plt = plot(
        t_sim, T_max_sim_C, 
        label="LPTN Simulation", 
        linewidth=2, color=:blue,
        xlabel="Time (s)", ylabel="Temperature (°C)",
        title="Nema et al. (2026) - $rate_name validation",
        legend=:bottomright, grid=true
    )
    
    # Overlay CFD validation data as scatter points
    scatter!(
        plt, t_val, y_val_C, 
        label="Nema 3D CFD Data", 
        markershape=:circle, color=:red, markersize=4, alpha=0.7
    )
    
    # Annotate plot with calculated error metrics
    annotate!(plt, [(t_hg_end*0.1, maximum(y_val_C)*0.95, text("RMSE: $(round(rmse, digits=2)) °C\nR²: $(round(r2, digits=3))", 10, :left))])
    
    # Return completed plot object
    return plt
end

# Execute validation sequence for 1C discharge
p_1c = run_and_validate("1C", q_gen_1c, t_hg_1c[end], "8a_1C.csv")
# Execute validation sequence for 5C discharge
p_5c = run_and_validate("5C", q_gen_5c, t_hg_5c[end], "8a_5C.csv")

# Display combined plot layout
display(plot(p_1c, p_5c, layout=(1, 2), size=(1000, 400)))