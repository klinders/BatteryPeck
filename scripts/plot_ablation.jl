# ==============================================================================
# plot_ablation.jl
# Generate SoH overlapping plots, dense load traces, and summary tables.
# Contains:
# 1. generate_dense_24h_csv: Automatically run a 1-day dense simulation if the CSV is missing
# 2. extract_ablation_metrics: Calculate final SoH, peak temperature, and max gradient
# 3. plot_ablation_dashboard: Render combined thesis figures and print terminal summary
# ==============================================================================

using Revise
using CSV
using DataFrames
using Plots
using Statistics
using Printf

# Import simulation libraries to generate the dense CSV if missing
using BatteryToolkit
using ModelingToolkit
using OrdinaryDiffEq

Revise.revise()

# ==============================================================================
# --- USER CONFIGURATION DASHBOARD ---
# ==============================================================================
CONFIG = (
    timestamp_folder = "ablation_run_1783249813", # Change to your actual folder name!
    c_rate = 0.75,                                 # Easily switch between 0.1 and 0.75
    
    # Internal pack constants for the generator
    rows = 4,
    cols = 7,
    cell_cap_ah = 5.0
)

# Standard MATLAB Colour Palette
MATLAB_COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]
# ==============================================================================

"""
    generate_dense_24h_csv(cfg, out_path)

Automatically run a 1-day dense simulation if the CSV is missing.
Constructs Scenario 5 (Coupled, Thermal) and simulates the first 24 hours to capture a high-fidelity trace.
"""
function generate_dense_24h_csv(cfg, out_path::String)
    println("\n[!] Dense CSV not found. Generating now for Scenario 5 at $(cfg.c_rate)C (High Fidelity)...")
    
    elec_params = Chen2020()
    elec_params.Vmin = 2.5
    elec_params.Vmax = 4.2
    
    soc_init = 0.65 
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    pack_1C_amps = cfg.cell_cap_ah * cfg.cols 
    pack_v_max = 4.2 * cfg.rows
    pack_v_min = 2.5 * cfg.rows

    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols, wall_thickness=0.001)
    nema_params = build_pack_parameters(
        coolant=:water_nema2026, 
        ambient_temp=298.15, 
        inlet_temp=298.15, 
        flow_rate=1e-6
    )
    
    geom = build_pack_geometry(rows=cfg.rows, cols=cfg.cols, cell_pitch=0.025)
    cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows, cfg.cols; verbose=false, geom_sigma=0.005)
    
    reactive_strat = ReactiveTMS(ambient_temp=298.15, T_high=35.0, T_low=32.0)
    anticipative_strat = AnticipativeTMS(ambient_temp=298.15, T_high=35.0, T_low=32.0, cell_load_threshold=5.0, lookahead_window=300.0)
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=Dict(:aging => :reactive, :rest => :reactive))
    
    schedule = Step[
        TargetCurrentStep(-cfg.c_rate * pack_1C_amps, 0.35, pack_v_min, 2.5 * 3600.0, :aging), 
        RestStep(1.0 * 3600.0),
        TargetCurrentStep(cfg.c_rate * pack_1C_amps, 0.65, pack_v_max, 2.5 * 3600.0, :aging),   
        RestStep(18.0 * 3600.0)
    ]
    exp = Experiment(schedule)
    
    df = simulate_pack!(cosim, exp, nema_params, tms; 
        total_time=24.0 * 3600.0, save_csv=false, print_prefix="[Dense single cycle] ",
        force_even_current=false, is_isothermal=false, geom_sigma=0.005,
        m_active=1e-6, m_passive=1e-6, heat_multiplier=1.0,
        
        # ----------------------------------------------------------------------
        # HIGH FIDELITY TIME-STEPPING (Quartered dt_max and alpha)
        # ----------------------------------------------------------------------
        opt_w_ratio=Dict{Symbol, Float64}(:aging => 2.0, :rest => 1.0), 
        opt_alpha_elec=Dict{Symbol, Float64}(:aging => 2.5, :rest => 2.5), 
        opt_dt_max_elec=Dict{Symbol, Float64}(:aging => 10.0, :rest => 100.0),
        opt_alpha_therm=Dict{Symbol, Float64}(:aging => 1, :rest => 1),
        opt_dt_max_therm=Dict{Symbol, Float64}(:aging => 10.0, :rest => 100.0),
        # ----------------------------------------------------------------------
        
        dense_logging=true, sparse_logging=false, ultra_sparse_logging=false, verbose=true
    )
    
    CSV.write(out_path, df)
    println("\n>>> Successfully generated and saved $out_path")
end

"""
    extract_ablation_metrics(df)

Calculate final SoH, peak temperature, and max gradient dynamically.
"""
function extract_ablation_metrics(df::DataFrame)
    soh_cols = filter(n -> startswith(n, "SoH_Cell_"), names(df))
    final_soh = mean(Float64.(Matrix(df[end:end, soh_cols])))
    peak_t = maximum(df.Max_Temp_C)
    
    t_cols = filter(n -> startswith(n, "T_Core_Cell_"), names(df))
    if length(t_cols) > 1
        T_mat = Float64.(Matrix(df[!, t_cols]))
        max_dt = maximum(maximum(T_mat, dims=2) .- minimum(T_mat, dims=2))
    else
        max_dt = 0.0
    end
    
    return final_soh, peak_t, max_dt, soh_cols
end

"""
    plot_ablation_dashboard(cfg)

Render combined thesis figures and print terminal summary.
"""
function plot_ablation_dashboard(cfg)
    base_dir = joinpath(pwd(), "results", cfg.timestamp_folder, "crate_$(cfg.c_rate)")
    if !isdir(base_dir)
        error("[!] Cannot find directory: $base_dir\nPlease check your timestamp and c_rate in CONFIG.")
    end

    # Use a unique name so the file doesn't overwrite if you switch C-rates
    dense_path = joinpath(base_dir, "dense_ablation_24h_$(cfg.c_rate)C.csv")
    if !isfile(dense_path)
        Base.invokelatest(generate_dense_24h_csv, cfg, dense_path)
    end

    println("\n=======================================================")
    println(">>> ABLATION STUDY SUMMARY ($(cfg.c_rate)C)")
    println("=======================================================")
    
    header = @sprintf("%-16s | %-15s | %-15s | %-15s", "Scenario", "Final SoH (%)", "Peak temp (°C)", "Max ΔT (°C)")
    println(header)
    println(repeat("-", length(header)))
    
    p_soh = plot(title="Ablation study: degradation ($(cfg.c_rate)C)", xlabel="Time [days]", ylabel="Pack average SoH [%]", legend=:bottomleft, lw=2)
    labels = ["S1 (iso, dec)", "S2 (therm, dec)", "S3 (iso, coup)", "S4 (therm, dec)", "S5 (therm, coup)", "S6 (double var)"]
    
    valid_plots = 0
    for i in 1:6
        scn_csv = joinpath(base_dir, "scn_$i", "ablation_1_year_master.csv")
        if !isfile(scn_csv)
            scn_csv = joinpath(base_dir, "scn_$i", "master_results.csv")
            if !isfile(scn_csv) continue end
        end
        
        df = CSV.read(scn_csv, DataFrame)
        final_soh, peak_t, max_dt, soh_cols = Base.invokelatest(extract_ablation_metrics, df)
        
        line = @sprintf("%-16s | %14.2f%% | %14.2f | %14.3f", labels[i], final_soh, peak_t, max_dt)
        println(line)
        
        avg_soh = mean(Float64.(Matrix(df[!, soh_cols])), dims=2)[:, 1]
        
        # Apply MATLAB colors cyclically
        plot!(p_soh, df.Time_s ./ (24*3600), avg_soh, label=labels[i], lw=2, color=MATLAB_COLORS[i])
        valid_plots += 1
    end
    
    if valid_plots > 0
        soh_out = "ablation_soh_comparison_$(cfg.c_rate)C.svg"
        savefig(p_soh, soh_out)
        println("\n>>> Exported SoH overlap plot to '$soh_out'")
    else
        println("\n[!] No scenario CSVs found in $base_dir to plot SoH.")
    end
    
    # ----------------------------------------------------------------------
    # Render Dense Dashboard (Voltage/Current Dual Axis + Temperature)
    # ----------------------------------------------------------------------
    df_dense = CSV.read(dense_path, DataFrame)
    t_hrs = df_dense.Time_s ./ 3600.0
    
    # Isolate the x-axis to the first 7.5 hours (2.5h + 1h + 2.5h + 1.5h settling buffer)
    cycle_xlims = (0.0, 7.5)
    
    # Dual Axis Plot: Voltage and Current
    p_vi = plot(t_hrs, df_dense.Pack_Voltage_V, title="Pack voltage and current ($(cfg.c_rate)C)", ylabel="Voltage [V]", xlabel="", label="Voltage", color=MATLAB_COLORS[1], lw=2, legend=:topleft, xlims=cycle_xlims)
    p_i = twinx(p_vi)
    plot!(p_i, t_hrs, df_dense.Pack_Current_A, ylabel="Current [A]", label="Current", color=MATLAB_COLORS[2], lw=2, legend=:bottomright, xlims=cycle_xlims)
    
    # Temperature Plot
    p_t = plot(t_hrs, df_dense.Max_Temp_C, title="Maximum cell temperature", ylabel="Temperature [°C]", xlabel="Time [hours]", label="Temperature", color=MATLAB_COLORS[7], lw=2, legend=false, xlims=cycle_xlims)
    
    dash_dense = plot(p_vi, p_t, layout=(2,1), size=(800, 600), margin=6Plots.mm)
    dense_out = "ablation_dense_cycle_$(cfg.c_rate)C.svg"
    savefig(dash_dense, dense_out)
    println(">>> Exported single cycle dense traces to '$dense_out'")
end

Base.invokelatest(plot_ablation_dashboard, CONFIG)