# ==============================================================================
# run_sensitivity.jl
# One-At-A-Time (OAT) Sensitivity Analysis for Battery Pack Thermal LPTN
# Contains:
# 1. calculate_metrics: Extract thermal peak, max gradient, and thermal dose from dataframe
# 2. run_sensitivity_scenario: Execute single explicit co-simulation with specific thermal parameters
# 3. run_sensitivity_matrix: Orchestrate baseline and sweep evaluations and generate tornado plots
# ==============================================================================

using Plots
using Revise
using BatteryToolkit
using CSV, DataFrames
using SciMLBase
using Printf
using Statistics

Revise.revise()

# Define global configuration dashboard and baseline parameter dict
CONFIG = (
    rows_series = 4,    
    cols_parallel = 7,  
    cell_capacity_ah = 5.0,
    cell_nominal_v = 3.7,
    ambient_temp = 298.15,

    baseline = Dict{Symbol, Float64}(
        :R_axial_val => 0.25,
        :R_radial_val => 0.20,
        :R_contact_val => 4.0,
        :h_conv => 5.0,
        :casing_th => 0.01,
        :cell_pitch => 0.025
    ),
    
    sweep_variance = 0.10
)

"""
    calculate_metrics(df, num_cells)

Extract thermal peak, max gradient, and thermal dose from dataframe.

Process output time-series data to isolate absolute peak temperature, maximum spatial temperature gradient across all cells, and integral of thermal dose above ambient.

# Arguments
- `df`: DataFrame containing transient simulation history
- `num_cells`: Integer count of total cells in pack

# Returns
- Tuple containing peak temperature, maximum delta T, and thermal dose
"""
function calculate_metrics(df::DataFrame, num_cells::Int)
    # Extract absolute peak core temperature from recorded history
    peak_T = maximum(df.Max_Temp_C)
    
    # Isolate individual cell temperature columns and calculate maximum spatial divergence at any given time step
    T_cols = [Symbol("T_Core_Cell_$i") for i in 1:num_cells]
    T_matrix = Matrix(df[!, T_cols])
    delta_Ts = maximum(T_matrix, dims=2) .- minimum(T_matrix, dims=2)
    max_delta_T = maximum(delta_Ts)
    
    # Integrate average pack temperature above ambient over total simulation duration to determine cumulative thermal dose
    dose = 0.0
    for i in 1:(nrow(df)-1)
        dt = df.Time_s[i+1] - df.Time_s[i]
        T_avg = mean(T_matrix[i, :])
        dose += max(0.0, T_avg - 25.0) * dt
    end
    
    return peak_T, max_delta_T, dose
end

"""
    run_sensitivity_scenario(params, cfg, label)

Execute single explicit co-simulation with specific thermal parameters.

Initialise geometry, parameters, and orchestrator using provided dictionary, run discharge schedule, and compute summary metrics.

# Arguments
- `params`: Dictionary containing physical and thermal parameters for current run
- `cfg`: NamedTuple containing global configuration dashboard
- `label`: String identifier for current run

# Returns
- Tuple containing calculated thermal metrics
"""
function run_sensitivity_scenario(params::Dict{Symbol, Float64}, cfg, label::String)
    # Establish base electrochemical parameters and instantiate safe initial chemical state avoiding boundary spikes
    elec_params = Chen2020()
    elec_params.Vmin = 2.5
    elec_params.Vmax = 4.2
    
    soc_init = 0.95 
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    # Inject specific sensitivity parameters into geometric and fluid property builders
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols_parallel, wall_thickness=0.001)
    nema_params = build_pack_parameters(
        coolant = :water_nema2026,
        ambient_temp = cfg.ambient_temp,
        inlet_temp = cfg.ambient_temp,
        flow_rate = 1e-6,
        cell_pitch = params[:cell_pitch],
        h_conv = params[:h_conv],
        casing_th = params[:casing_th]
    )
    
    # Construct complete explicit orchestrator using overriding thermal resistance keyword arguments
    geom = build_pack_geometry(rows=cfg.rows_series, cols=cfg.cols_parallel, cell_pitch=params[:cell_pitch])
    cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows_series, cfg.cols_parallel; 
        verbose=false, geom_sigma=0.0,
        R_axial_val=params[:R_axial_val], R_radial_val=params[:R_radial_val], R_contact_val=params[:R_contact_val]
    )
    
    # Formulate dummy reactive cooling strategy to strictly enforce unmitigated thermal drift
    reactive_strat = ReactiveTMS(ambient_temp=cfg.ambient_temp, T_high=35.0, T_low=32.0)
    anticipative_strat = AnticipativeTMS(ambient_temp=cfg.ambient_temp, T_high=35.0, T_low=32.0, cell_load_threshold=5.0, lookahead_window=300.0)
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=Dict(:smooth => :reactive, :rest => :reactive))
    
    # Design experimental schedule discharging at half nominal rate until reaching target threshold followed by rest phase
    pack_1C_amps = cfg.cell_capacity_ah * cfg.cols_parallel 
    pack_v_min = 2.5 * cfg.rows_series
    
    schedule = Step[
        TargetCurrentStep(-0.5 * pack_1C_amps, 0.20, pack_v_min, 3.0 * 3600.0, :smooth), 
        RestStep(3600.0)
    ]
    exp = Experiment(schedule)
    total_time = sum([s.period for s in schedule])
    
    # Execute full explicit co-simulation leveraging pareto optimised time steps
    df = simulate_pack!(cosim, exp, nema_params, tms; 
        total_time=total_time, save_csv=false, print_prefix="[$label] ",
        force_even_current=false, is_isothermal=false, geom_sigma=0.0,
        m_active=1e-6, m_passive=1e-6, heat_multiplier=1.0,
        
        opt_w_ratio=Dict{Symbol, Float64}(:smooth => 1.75, :rest => 1.0), 
        opt_alpha_elec=Dict{Symbol, Float64}(:smooth => 30.0, :rest => 10.0), 
        opt_dt_max_elec=Dict{Symbol, Float64}(:smooth => 1000.0, :rest => 3600.0),
        opt_alpha_therm=Dict{Symbol, Float64}(:smooth => 0.1, :rest => 0.1),
        opt_dt_max_therm=Dict{Symbol, Float64}(:smooth => 558.3, :rest => 3600.0),
        
        dense_logging=false, sparse_logging=true, ultra_sparse_logging=false, verbose=true
    )
    
    # Extract numerical metrics and clean memory to prevent overflow during sequential matrix execution
    metrics = calculate_metrics(df, cfg.rows_series * cfg.cols_parallel)
    cosim = nothing; GC.gc() 
    return metrics
end

"""
    run_sensitivity_matrix(cfg)

Orchestrate baseline and sweep evaluations and generate tornado plots.

Execute ground truth simulation, sequentially perturb each thermal parameter by predefined variance, log results to disk, and render horizontal bar charts.

# Arguments
- `cfg`: NamedTuple containing global configuration dashboard
"""
function run_sensitivity_matrix(cfg)
    # Acquire ground truth measurements running baseline dictionary
    println("\n>>> INITIATING THERMAL SENSITIVITY SWEEP MATRIX")
    base_peak, base_grad, base_dose = run_sensitivity_scenario(cfg.baseline, cfg, "BASELINE")
    
    # Initialise arrays to store parameter deviations for graphical rendering and dataframe storage
    param_keys = collect(keys(cfg.baseline))
    param_names = String.(param_keys)
    
    peak_low_diffs = Float64[]; peak_high_diffs = Float64[]
    grad_low_diffs = Float64[]; grad_high_diffs = Float64[]
    dose_low_diffs = Float64[]; dose_high_diffs = Float64[]
    
    # Loop over parameters to evaluate low and high variance scenarios
    for key in param_keys
        base_val = cfg.baseline[key]
        low_val = base_val * (1.0 - cfg.sweep_variance)
        high_val = base_val * (1.0 + cfg.sweep_variance)
        
        # Execute low boundary deviation
        low_params = copy(cfg.baseline)
        low_params[key] = low_val
        l_peak, l_grad, l_dose = run_sensitivity_scenario(low_params, cfg, "$(key) LOW")
        
        # Execute high boundary deviation
        high_params = copy(cfg.baseline)
        high_params[key] = high_val
        h_peak, h_grad, h_dose = run_sensitivity_scenario(high_params, cfg, "$(key) HIGH")
        
        # Calculate percentage divergence relative to established ground truth
        push!(peak_low_diffs, ((l_peak - base_peak) / base_peak) * 100.0)
        push!(peak_high_diffs, ((h_peak - base_peak) / base_peak) * 100.0)
        
        push!(grad_low_diffs, ((l_grad - base_grad) / base_grad) * 100.0)
        push!(grad_high_diffs, ((h_grad - base_grad) / base_grad) * 100.0)
        
        push!(dose_low_diffs, ((l_dose - base_dose) / base_dose) * 100.0)
        push!(dose_high_diffs, ((h_dose - base_dose) / base_dose) * 100.0)
    end
    
    # Compile discrete results into structured dataframe
    results_df = DataFrame(
        Parameter = param_names,
        Peak_T_Low_Pct = peak_low_diffs,
        Peak_T_High_Pct = peak_high_diffs,
        Max_Grad_Low_Pct = grad_low_diffs,
        Max_Grad_High_Pct = grad_high_diffs,
        Dose_Low_Pct = dose_low_diffs,
        Dose_High_Pct = dose_high_diffs
    )
    
    # Construct output directory to store sensitivity artifacts
    timestamp = round(Int, time())
    results_dir = joinpath(pwd(), "results", "sensitivity_run_$timestamp")
    mkpath(results_dir)
    
    csv_path = joinpath(results_dir, "sensitivity_results.csv")
    txt_path = joinpath(results_dir, "sensitivity_results.txt")
    img_path = joinpath(results_dir, "sensitivity_tornado_dashboard.png")
    
    CSV.write(csv_path, results_df)
    
    # Open file and standard output simultaneously to log formatted results
    open(txt_path, "w") do io
        for out in (stdout, io)
            println(out, "\n>>> SENSITIVITY ANALYSIS RESULTS (±$(round(Int, cfg.sweep_variance * 100))% Variance)")
            println(out, "Baseline Peak T: $(round(base_peak, digits=2)) °C")
            println(out, "Baseline Max ΔT: $(round(base_grad, digits=3)) °C")
            println(out, "Baseline Dose:   $(round(base_dose, digits=0)) °C-s\n")
            
            header = @sprintf("%-20s | %11s | %11s | %11s | %11s | %11s | %11s", "Parameter", "Peak T Low", "Peak T High", "Grad Low", "Grad High", "Dose Low", "Dose High")
            println(out, header)
            println(out, repeat("-", length(header)))
            
            for i in 1:nrow(results_df)
                line = @sprintf("%-20s | %10.2f%% | %10.2f%% | %10.2f%% | %10.2f%% | %10.2f%% | %10.2f%%", 
                    results_df.Parameter[i], 
                    results_df.Peak_T_Low_Pct[i], results_df.Peak_T_High_Pct[i],
                    results_df.Max_Grad_Low_Pct[i], results_df.Max_Grad_High_Pct[i],
                    results_df.Dose_Low_Pct[i], results_df.Dose_High_Pct[i]
                )
                println(out, line)
            end
        end
    end
    
    # Render horizontal bar charts correlating variance magnitude against physical parameters using numerical mapping
    y_pos = 1:length(param_names)
    y_ticks = (y_pos, param_names)
    
    p1 = bar(y_pos, peak_low_diffs, orientation=:h, yticks=y_ticks, label="-10% Variance", color=:blue, title="Peak Temperature", xlabel="% Change from Baseline", legend=:bottomright)
    bar!(p1, y_pos, peak_high_diffs, orientation=:h, label="+10% Variance", color=:red, ylims=(0.5, 6.5))
    
    p2 = bar(y_pos, grad_low_diffs, orientation=:h, yticks=y_ticks, label="", color=:blue, title="Max Temperature Gradient", xlabel="% Change from Baseline")
    bar!(p2, y_pos, grad_high_diffs, orientation=:h, label="", color=:red, ylims=(0.5, 6.5))
    
    p3 = bar(y_pos, dose_low_diffs, orientation=:h, yticks=y_ticks, label="", color=:blue, title="Thermal Dose", xlabel="% Change from Baseline")
    bar!(p3, y_pos, dose_high_diffs, orientation=:h, label="", color=:red, ylims=(0.5, 6.5))
    
    # Compile individual plots into master dashboard and export image
    dashboard = plot(p1, p2, p3, layout=(3,1), size=(1000, 1200), margin=8Plots.mm)
    display(dashboard)
    savefig(dashboard, img_path)
    
    println("\n>>> Dashboard exported to: $img_path")
    println(">>> Raw data exported to: $csv_path")
end

Base.invokelatest(run_sensitivity_matrix, CONFIG)