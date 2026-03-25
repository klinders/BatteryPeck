# ==============================================================================
# validation_nema_13f.jl
# Validates pack temperature uniformity (Delta T) against Nema et al. (2026)
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
    # Read raw data from specified file path
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
    
    # Calculate root mean square error
    rmse = sqrt(sum((y_sim_at_val .- y_val_trimmed).^2) / length(y_val_trimmed))
    
    # Calculate mean of validation values
    y_mean = sum(y_val_trimmed) / length(y_val_trimmed)
    
    # Calculate total sum of squares
    ss_tot = sum((y_val_trimmed .- y_mean).^2)
    # Calculate residual sum of squares
    ss_res = sum((y_val_trimmed .- y_sim_at_val).^2)
    
    # Calculate coefficient of determination
    r2 = 1.0 - (ss_res / ss_tot)
    
    # Return calculated error metrics
    return rmse, r2
end

# Define path to data directory
data_dir = joinpath(@__DIR__, "..", "data", "Nema2026")

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

# Sets physical environment and baseline boundaries
# Editable values:
# `v_coolant`: Alters fluid velocity. directly impacts downstream heating and Delta T magnitude.
v_coolant = 0.3 
# Define fluid density
rho_water = 998.0
# Define channel width
W_channel = 0.050 
# Define channel height
H_channel = 0.004 
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
    ambient_convection_coefficient = 17.0
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
    run_uniformity_validation(t_hg_end, val_filename)

Builds active LPTN system, runs solver, calculates maximum temperature gradient, and plots results.

Equations:
ΔT = T_max - T_min
Extracts instantaneous maximum and minimum core temperatures across all simulated cells.
"""
function run_uniformity_validation(t_hg_end, val_filename)
    # Print validation start message to console
    println("\nStarting uniformity validation: 5C discharge (Figure 13f)")
    
    # Construct and simplify ODE system while suppressing console warnings
    sys_simplified, cells = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        # Build base thermal system
        sys_base = build_pack_system(Symbol("sys_5C"), geom, val_params)
        # Extract cell components from system
        cls = [getproperty(sys_base, Symbol("cell_$i")) for i in 1:num_cells]
        # Apply specific heat generation function to each cell
        forcing_eqs = [cls[i].Q_volumetric_in.u ~ q_gen_5c(t) for i in 1:num_cells]
        # Couple base system with forcing equations
        full_sys = ODESystem(forcing_eqs, t; systems=[sys_base], name=Symbol("full_5C"))
        # Return structurally simplified system and cell array
        (structural_simplify(full_sys), cls) 
    end

    # Define ODE problem with jacobian and sparse matrix forms
    prob = ODEProblem(sys_simplified, [], (0.0, t_hg_end), jac=true, sparse=true)
    # Solve ODE problem and time execution
    @time sol = solve(prob, QNDF(), saveat=1.0, dtmax=1.0) 

    # Extract simulation time array
    t_sim = sol.t
    
    # Extracts temperatures for all 28 cells into matrix
    T_cells_K = [sol[cells[i].core_cap.T] for i in 1:num_cells]
    T_mat_K = hcat(T_cells_K...)
    
    # Calculates instantaneous peak and minimum temperatures
    T_max_sim_K = maximum(T_mat_K, dims=2)[:, 1]
    T_min_sim_K = minimum(T_mat_K, dims=2)[:, 1]
    
    # Calculates absolute temperature uniformity gradient
    Delta_T_sim = T_max_sim_K .- T_min_sim_K

    # Load validation data from specified filename
    t_val, y_val_delta_T = load_and_preprocess(joinpath(data_dir, val_filename))
    
    # Compute error metrics between simulation and validation data
    rmse, r2 = compute_metrics(t_sim, Delta_T_sim, t_val, y_val_delta_T)
    
    # Print calculated RMSE to console
    println("  -> RMSE: $(round(rmse, digits=3)) °C")
    # Print calculated R² to console
    println("  -> R²:   $(round(r2, digits=4))")

    # Initialise plot with simulation data
    plt = plot(
        t_sim, Delta_T_sim, 
        label="LPTN ΔT", 
        linewidth=2, color=:purple,
        xlabel="Time (s)", ylabel="Temperature Difference ΔT (°C)",
        title="Nema et al. (2026) - Fig 13f (Pack Uniformity)",
        legend=:bottomright, grid=true
    )
    # Overlay CFD validation data as scatter points
    scatter!(
        plt, t_val, y_val_delta_T, 
        label="Nema 3D CFD Data", 
        markershape=:circle, color=:orange, markersize=4, alpha=0.7
    )
    # Annotate plot with calculated error metrics
    annotate!(plt, [(t_hg_end*0.1, maximum(y_val_delta_T)*0.95, text("RMSE: $(round(rmse, digits=2)) °C\nR²: $(round(r2, digits=3))", 10, :left))])
    
    # Return completed plot object, solver solution, and cell array
    return plt, sol, cells
end

"""
    animate_pack_temperatures(sol, cells, geom, filename="pack_thermal_animation.gif")

Generates top-down spatial heat map animation of pack discharge.
Uses exact index slicing to prevent SciML interpolation jitter.
Applies dynamic colour limits to maximise spatial gradient visibility.
"""
function animate_pack_temperatures(sol, cells, geom, filename="pack_thermal_animation.gif")
    # Print animation start message to console
    println("\nGenerating spatial thermal animation...")
    
    # Extracts physical layout
    x_coords = [c[1] for c in geom.cell_coords]
    y_coords = [c[2] for c in geom.cell_coords]
    
    # Define horizontal padding
    x_pad = 0.02
    # Define vertical padding
    y_pad = 0.02
    # Calculate horizontal plot limits
    x_lims = (minimum(x_coords) - x_pad, maximum(x_coords) + x_pad)
    # Calculate vertical plot limits
    y_lims = (minimum(y_coords) - y_pad, maximum(y_coords) + y_pad)
    
    # Downsamples by extracting exact saved indices to prevent interpolation jitter
    total_steps = length(sol.t)
    # Calculate frame step size to maintain framerate
    step_size = max(1, div(total_steps, 150)) 
    # Create range of indices for animation frames
    frame_indices = 1:step_size:total_steps
    
    # Loop through selected frames to build animation
    anim = @animate for idx in frame_indices
        # Extract current time step
        t = sol.t[idx]
        
        # Extracts exact saved values directly from the array
        T_instant = [sol[cells[i].core_cap.T][idx] .- 273.15 for i in 1:length(cells)]
        
        # Dynamic colour limits to highlight the spatial gradient at this exact second
        frame_min = floor(minimum(T_instant))
        frame_max = ceil(maximum(T_instant))
        
        # Prevents perfectly uniform frames (like t=0) from crashing clims
        if frame_max == frame_min
            # Apply offset to prevent colour scale crash
            frame_max += 1.0
        end
        
        # Render spatial heat map for current frame
        scatter(
            x_coords, y_coords, 
            marker_z = T_instant,
            markershape = :circle, 
            markersize = 22, 
            color = :turbo, 
            clims = (frame_min, frame_max),
            xlims = x_lims,
            ylims = y_lims,
            framestyle = :box, 
            aspect_ratio = :equal,
            xlabel = "X Position (m)", 
            ylabel = "Y Position (m)",
            title = "LPTN Spatial Gradient | t = $(round(t, digits=1)) s",
            colorbar = true, 
            colorbar_title = "Temperature (°C)",
            legend = false, 
            margin = 8Plots.mm, 
            right_margin = 15Plots.mm, 
            size = (750, 450)
        )
    end
    
    # Compile frames into gif file
    gif(anim, filename, fps=15)
    # Print completion message to console
    println("Animation complete! Saved to local directory as: $filename")
end

"""
    plot_diagnostic_disco(sol, cells, t_sim)

Diagnoses high-frequency spatial and temporal numerical oscillations.
Generates a 3-panel plot showing temperature envelopes, extreme cell locations, and raw traces.
"""
function plot_diagnostic_disco(sol, cells, t_sim)
    # Print diagnostic generation message to console
    println("\nGenerating diagnostic oscillation plots...")
    
    # Extracts exact temperature matrix
    T_mat_K = hcat([sol[cells[i].core_cap.T] for i in 1:length(cells)]...)
    # Convert temperature matrix to celsius
    T_mat_C = T_mat_K .- 273.15
    
    # Calculates max and min for envelope
    T_max = maximum(T_mat_C, dims=2)[:, 1]
    T_min = minimum(T_mat_C, dims=2)[:, 1]
    
    # Identifies which specific cells are the hottest and coldest
    hot_idx = [argmax(T_mat_C[i, :]) for i in 1:length(t_sim)]
    cold_idx = [argmin(T_mat_C[i, :]) for i in 1:length(t_sim)]
    
    # Initialise temperature envelope plot
    p1 = plot(t_sim, T_max, label="Max Temp", color=:red, linewidth=2, 
              ylabel="Temp (°C)", title="1. Temperature Envelope", grid=true)
    # Overlay minimum temperature trace
    plot!(p1, t_sim, T_min, label="Min Temp", color=:blue, linewidth=2)
    
    # Initialise extreme cell location plot
    p2 = scatter(t_sim, hot_idx, label="Hottest Cell", color=:red, markersize=3, 
                 ylabel="Cell ID (1-28)", title="2. Location of Extremes", grid=true, legend=:right)
    # Overlay coldest cell locations
    scatter!(p2, t_sim, cold_idx, label="Coldest Cell", color=:blue, markersize=3)
    
    # Initialise individual cell trace plot
    p3 = plot(t_sim, T_mat_C[:, 1], label="Cell 1 (Inlet)", color=:cyan, linewidth=1.5,
              xlabel="Time (s)", ylabel="Temp (°C)", title="3. Individual Cell Traces", grid=true)
    # Overlay middle cell trace
    plot!(p3, t_sim, T_mat_C[:, 14], label="Cell 14 (Middle)", color=:green, linewidth=1.5)
    # Overlay outlet cell trace
    plot!(p3, t_sim, T_mat_C[:, 28], label="Cell 28 (Outlet)", color=:magenta, linewidth=1.5)
    
    # Render diagnostic plots in vertical layout
    display(plot(p1, p2, p3, layout=(3, 1), size=(800, 900)))
end

"""
    plot_transverse_cross_section(sol, cells)

Plots a cross-sectional temperature profile across the 4 rows of the battery pack
at the end of the 5C discharge to diagnose transverse thermal conduction.
"""
function plot_transverse_cross_section(sol, cells)
    # Print cross section generation message to console
    println("\nGenerating Transverse Cross-Section...")
    
    # We take a vertical slice down the middle of the pack (Column 4)
    # Row 1 (Inlet side), Row 2, Row 3, Row 4 (Outlet side)
    slice_indices = [4, 11, 18, 25] 
    
    # Extract final temperatures for selected slice
    T_end = [sol[cells[i].core_cap.T][end] .- 273.15 for i in slice_indices]
    # Define labels for transverse rows
    row_labels = ["Row 1\n(Cold Side)", "Row 2", "Row 3", "Row 4\n(Hot Side)"]
    
    # Render bar chart of transverse temperature profile
    p = bar(row_labels, T_end, 
            title="Transverse Temperature Profile (End of 5C)",
            ylabel="Temperature (°C)", 
            legend=false, 
            color=[:cyan, :skyblue, :orange, :red],
            ylims=(minimum(T_end) - 1.0, maximum(T_end) + 1.0),
            grid=true)
            
    # Add exact temperature labels on top of the bars
    for (i, t) in enumerate(T_end)
        # Render exact temperature value above bar
        annotate!(p, i, t + 0.1, text("$(round(t, digits=2)) °C", 10, :bottom))
    end
    
    # Render transverse profile plot
    display(p)
end

# Execute uniformity validation and extract results
p_uniformity, sol_out, cells_out = run_uniformity_validation(t_hg_5c[end], "13f_SSCC.csv")

# Display pack uniformity validation plot
display(p_uniformity)

# Passes the extracted variables safely into the animation engine
#animate_pack_temperatures(sol_out, cells_out, geom, "5C_serpentine_gradient.gif")

#plot_diagnostic_disco(sol_out, cells_out, sol_out.t)

#plot_transverse_cross_section(sol_out, cells_out)