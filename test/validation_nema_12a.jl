# ==============================================================================
# validation_nema_12a.jl
# Validates pack temperature vs. channel height against Nema et al. (2026) Fig 12a
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq 
using DataInterpolations
using DelimitedFiles
using Plots
using Plots.PlotMeasures
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

# Define path to validation data directory
data_dir = joinpath(@__DIR__, "..", "data", "Nema2026")

# Load and clean 5C heat generation data
t_hg_5c, y_hg_5c = load_and_preprocess(joinpath(data_dir, "7b_298K_HG.csv"))
# Create linear interpolator for 5C data
const hg_interp_5c = LinearInterpolation(y_hg_5c, t_hg_5c)

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

# Define channel width
W_channel = 0.002
# Define fluid density
rho_water = 998.0

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
    run_figure_12a_validation(t_hg_end)

Sweeps through channel heights, validates against CFD CSVs, 
and compiles a 5-panel MATLAB-styled plot.
"""
function run_figure_12a_validation(t_hg_end)
    # Define array of target heights
    heights_mm = [40, 50, 60, 70]
    # Map 70mm to 60mm CSV since they represent same plateaued thermal state
    csv_files = ["12a_40mm.csv", "12a_50mm.csv", "12a_60mm.csv", "12a_60mm.csv"]
    # Set active coolant velocity
    v_coolant = 0.3 
    
    # Define standard MATLAB colour palette
    matlab_colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]
    
    # Initialise main overview plot
    p_main = plot(
        xlabel="Time (s)", ylabel="Temperature (K)",
        title="All Heights Overview (v = 0.3 m/s)",
        legend=:topleft, grid=true, ylims=(295, 340)
    )
    
    # Initialise array for individual subplots
    sub_plots = []
    
    # Print validation start message to console
    println("\nStarting Channel Height Validation (Figure 12a)...")

    # Iterate through target heights
    for (i, h_mm) in enumerate(heights_mm)
        
        # Convert height to metres
        H_channel = h_mm / 1000.0
        # Calculate cross sectional area
        A_cross = W_channel * H_channel
        # Calculate mass flow rate for current iteration
        m_dot_current = rho_water * A_cross * v_coolant
        
        # Initialise cooling channel geometry
        tms_geom_current = TMSGeometry(
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
            tms_geometry = tms_geom_current,
            
            cell_gap_thickness = (0.025 - 0.021) / 2.0, 
            axial_potting_thickness = 0.005,
            casing_thickness = 0.01,
            
            ambient_temperature = 298.15,
            inlet_temperature = 298.15,
            mass_flow_rate = m_dot_current, 
            ambient_convection_coefficient = 5
        )
        
        # Construct and simplify ODE system while suppressing console warnings
        sys_simplified, cells = with_logger(ConsoleLogger(stderr, Logging.Error)) do
            # Build base thermal system
            sys_base = build_pack_system(Symbol("sys_h$(h_mm)"), geom, val_params)
            
            # Extract cell components from system
            cls = [getproperty(sys_base, Symbol("cell_$j")) for j in 1:num_cells]
            
            # Apply specific heat generation function to each cell
            forcing_eqs = [cls[j].Q_volumetric_in.u ~ q_gen_5c(t) for j in 1:num_cells]
            
            # Couple base system with forcing equations
            full_sys = ODESystem(forcing_eqs, t; systems=[sys_base], name=Symbol("full_h$(h_mm)"))
            
            # Return structurally simplified system and cell array
            (structural_simplify(full_sys), cls) 
        end

        # Define ODE problem with jacobian and sparse matrix forms
        prob = ODEProblem(sys_simplified, [], (0.0, t_hg_end), jac=true, sparse=true)
        
        # Solve ODE problem and extract results
        sol = solve(prob, QNDF(), saveat=1.0, dtmax=1.0) 

        # Extract cell temperatures in kelvin
        T_cells_K = [sol[cells[j].core_cap.T] for j in 1:num_cells]
        
        # Find maximum cell temperature at each time step
        T_max_sim_K = maximum(hcat(T_cells_K...), dims=2)[:, 1]

        # Load validation data from corresponding file
        t_val, y_val_K = load_and_preprocess(joinpath(data_dir, csv_files[i]))
        
        # Compute error metrics between simulation and validation data
        rmse, r2 = compute_metrics(sol.t, T_max_sim_K, t_val, y_val_K)
        
        # Print calculated metrics to console
        println("  -> h = $h_mm mm  |  RMSE: $(round(rmse, digits=3)) K  |  R²: $(round(r2, digits=4))")

        # Overlay simulation data on main plot
        plot!(p_main, sol.t, T_max_sim_K, label="Sim $h_mm mm", color=matlab_colors[i], linewidth=2)
        # Overlay CFD validation data on main plot
        plot!(p_main, t_val, y_val_K, label="CFD $h_mm mm", color=matlab_colors[i], linestyle=:dash, linewidth=2)

        # Create dedicated subplot for current height
        p_sub = plot(
            sol.t, T_max_sim_K, 
            label="LPTN", color=matlab_colors[i], linewidth=2, 
            title="Height = $h_mm mm", xlabel="Time (s)", ylabel="Temp (K)", 
            grid=true, ylims=(295, 340), legend=:topleft
        )
        
        # Overlay validation data on subplot
        plot!(p_sub, t_val, y_val_K, label="CFD Data", color=matlab_colors[i], linestyle=:dash, linewidth=2)
        
        # Annotate subplot with calculated error metrics
        annotate!(p_sub, [(t_hg_end * 0.55, 303, text("RMSE: $(round(rmse, digits=2)) K\nR²: $(round(r2, digits=3))", 10, :left))])
        
        # Append completed subplot to array
        push!(sub_plots, p_sub)
    end
    
    # Compile individual plots into final grid layout leaving bottom right spot empty
    plt_final = plot(p_main, sub_plots..., layout=(2, 3), size=(1400, 800), margin=5Plots.mm)
    
    # Return completed figure
    return plt_final
end

# Execute height sweep validation
p_sweep = run_figure_12a_validation(t_hg_5c[end])

# Display final layout in plot pane
display(p_sweep)