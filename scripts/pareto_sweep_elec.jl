# ==============================================================================
# pareto_sweep_elec.jl
# Configurable regime pareto tensor using explicit orchestrator logic.
# Contains:
# 1. calc_rmse: Calculate root mean square error
# 2. calc_r2: Calculate coefficient of determination
# 3. load_fcr_data: Load frequency containment reserve data or generate synthetic profile
# 4. run_elec_pareto_tensor: Execute parameter sweep to identify optimal electrical solver configuration
# ==============================================================================

using BatteryToolkit
using CSV, DataFrames, DataInterpolations
using Plots
using Statistics
using Printf

# Configure architecture and scaling, thermal management parameters, and pareto bounds
CONFIG = (
    rows_series = 1,
    cols_parallel = 1,
    cell_capacity_ah = 5.0,
    cell_nominal_v = 3.7,
    id4_energy_wh = 62000.0,
    id4_fcr_power_w = 11000.0,
    
    tms_m_active  = 0.3,
    tms_m_passive = 1e-6,
    tms_t_high    = 35.0,
    tms_t_low     = 32.0,
    tms_lookahead_s = 300.0,
    tms_anticipative_I_thresh = 5.0,
    tms_regime_map = Dict(:grid => :anticipative, :drive => :reactive, :aging => :reactive, :smooth => :reactive, :rest => :reactive),

    target_rmse = Dict(:smooth => 50.0, :drive => 50.0, :grid => 50.0, :aging => 150.0),
    target_r2   = Dict(:smooth => 0.95, :drive => 0.95, :grid => 0.95, :aging => 0.85),
    
    alpha_range = Dict(:smooth => (10.0, 30.0), :drive => (40.0, 80.0), :grid => (1.0, 10.0), :aging => (5.0, 25.0)),
    w_range     = Dict(:smooth => (1.0, 4.0), :drive => (2.0, 5.0), :grid => (1.0, 4.0), :aging => (1.0, 5.0)),
    dt_range    = Dict(:smooth => (500.0, 1500.0), :drive => (0.5, 2.0), :grid => (5.0, 15.0), :aging => (200.0, 600.0)),
    
    # Define grid points per axis
    pts = 5 
)

"""
    calc_rmse(y_true, y_pred)

Calculate root mean square error.

# Arguments
- `y_true`: Array of ground truth numerical values
- `y_pred`: Array of predicted numerical values

# Returns
- Error metric evaluated over specified domain
"""
function calc_rmse(y_true, y_pred) 
    return sqrt(mean((y_true .- y_pred).^2)) 
end

"""
    calc_r2(y_true, y_pred)

Calculate coefficient of determination.

# Arguments
- `y_true`: Array of ground truth numerical values
- `y_pred`: Array of predicted numerical values

# Returns
- Correlation fraction representing predictive accuracy
"""
function calc_r2(y_true, y_pred)
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return 1.0 - (ss_res / (ss_tot + 1e-12))
end

"""
    load_fcr_data(duration)

Load frequency containment reserve data or generate synthetic profile.

# Arguments
- `duration`: Target playback duration

# Returns
- Tuple containing temporal array and coupled frequency signal
"""
function load_fcr_data(duration)
    # Initialise arrays and attempt file read to parse RTE dataset
    t_data = Float64[]
    f_data = Float64[]
    try
        df = CSV.read(joinpath("data", "V2G", "RTE_Frequence_2024", "RTE_Frequence_2024_01.txt"), DataFrame; delim=';')
        freq_col = names(df)[occursin.(r"freq"i, names(df))][1]
        
        raw_f = df[!, freq_col]
        f_data = map(raw_f) do val
            if ismissing(val) return 50.0 
            elseif val isa Number return Float64(val)
            elseif val isa AbstractString return tryparse(Float64, replace(strip(val), "," => ".")) === nothing ? 50.0 : tryparse(Float64, replace(strip(val), "," => "."))
            else return 50.0 end
        end
        f_data = convert(Vector{Float64}, f_data)
        t_data = collect(0.0 : 10.0 : (length(f_data)-1)*10.0)
    catch
        # Generate simulated harmonic sine series upon extraction failure
        t_data = collect(0.0 : 10.0 : duration)
        f_data = 50.0 .+ 0.15 .* sin.(t_data ./ 300.0) .+ 0.05 .* randn(length(t_data))
    end
    return t_data, f_data
end

"""
    run_elec_pareto_tensor(cfg)

Execute parameter sweep to identify optimal electrical solver configuration.

# Arguments
- `cfg`: Defined pareto control map

# Returns
- Nothing
"""
function run_elec_pareto_tensor(cfg)
    # Define scaling attributes and assemble thermal package architecture
    pack_ah = cfg.cols_parallel * cfg.cell_capacity_ah
    pack_wh = pack_ah * (cfg.rows_series * cfg.cell_nominal_v)
    energy_scale = pack_wh / cfg.id4_energy_wh
    
    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2

    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols_parallel, wall_thickness=0.001)
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = cfg.tms_m_passive, ambient_convection_coefficient = 5
    )
    
    reactive_strat = ReactiveTMS(T_high = cfg.tms_t_high, T_low = cfg.tms_t_low)
    anticipative_strat = AnticipativeTMS(
        T_high = cfg.tms_t_high, T_low = cfg.tms_t_low, 
        cell_load_threshold = cfg.tms_anticipative_I_thresh, 
        lookahead_window = cfg.tms_lookahead_s
    )
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=cfg.tms_regime_map)
    
    geom = build_pack_geometry(rows=cfg.rows_series, cols=cfg.cols_parallel, cell_pitch=0.025)

    regimes = [:smooth, :drive, :grid, :aging]
    PACK_1C = pack_ah

    # Process operational regime configurations sequentially
    for regime in regimes
        println("\n=======================================================")
        println(" INITIATING PARETO SWEEP: REGIME [$regime]")
        println("=======================================================")
        
        soc_init = (regime == :smooth) ? 0.45 : 0.65
        elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
        elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊
        
        exp = nothing
        duration = 3600.0
        
        # Branch logical layout parsing based on provided regime block
        if regime == :smooth
            exp = Experiment([TargetCurrentStep(0.25 * PACK_1C, 0.65, 4.2 * cfg.rows_series, duration)])
            
        elseif regime == :drive
            wltp_path = joinpath("data", "V2G", "driving_power_wltp.csv")
            if isfile(wltp_path)
                f_wltp = CSV.File(wltp_path) |> Tables.matrix
                t_wltp = f_wltp[:,1]
                p_wltp = -f_wltp[:,2] .* energy_scale
                dt_wltp = diff(t_wltp)
                tend = findfirst(t_wltp .>= duration)
                if isnothing(tend) tend = length(dt_wltp) end
                exp = Experiment([DriveStep(Any[dt_wltp[1:tend], p_wltp[1:tend]], duration)])
            else
                exp = Experiment([CurrentStep(-0.5 * PACK_1C, duration)])
            end
            
        elseif regime == :grid
            t_data, f_data = load_fcr_data(duration)
            max_v2g_pwr = cfg.id4_fcr_power_w * energy_scale
            exp = Experiment([FCRStep(t_data, f_data, max_v2g_pwr, duration)])
            
        elseif regime == :aging
            duration = 24.0 * 3600.0
            exp = Experiment([
                TargetCurrentStep(-0.5 * PACK_1C, 0.20, 2.5 * cfg.rows_series, 2.5 * 3600.0), 
                RestStep(1.0 * 3600.0),                                     
                TargetCurrentStep(0.5 * PACK_1C, 0.65, 4.2 * cfg.rows_series, 2.5 * 3600.0),  
                RestStep(18.0 * 3600.0)                                     
            ])
        end
        
        cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows_series, cfg.cols_parallel; verbose=false)
        
        # Build baseline benchmark solution targeting highly restricted tolerance bounds
        println("[!] Generating Ground Truth Benchmark (dt = 0.1s)...")
        truth_start = time()
        df_truth = simulate_pack!(cosim, exp, nema_params, tms; 
            total_time=duration, is_isothermal=true, save_csv=false, dense_logging=true, verbose=false,
            m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive,
            opt_alpha_elec=Dict(regime => 0.0), opt_dt_max_elec=Dict(regime => 0.1)
        )
        truth_dur = time() - truth_start
        t_true = df_truth.Time_s
        v_true = df_truth.Pack_Voltage_V
        
        println("    -> Benchmark generated in $(round(truth_dur, digits=2)) seconds.")
        
        best_alpha, best_w, best_dt, best_speedup, best_rmse, best_r2 = 0.0, 0.0, 0.0, 0.0, Inf, 0.0
        lowest_rmse_ever, best_fail_a, best_fail_w, best_fail_dt, best_fail_speedup, best_fail_r2 = Inf, 0.0, 0.0, 0.0, 0.0, 0.0
        
        println("[!] Sweeping parameter tensor...")
        
        target_rmse = cfg.target_rmse[regime]
        target_r2 = cfg.target_r2[regime]
        
        alpha_grid = collect(range(cfg.alpha_range[regime][1], cfg.alpha_range[regime][2], length=cfg.pts))
        w_grid = collect(range(cfg.w_range[regime][1], cfg.w_range[regime][2], length=cfg.pts))
        dt_grid = collect(range(cfg.dt_range[regime][1], cfg.dt_range[regime][2], length=cfg.pts))
        
        total_runs = length(alpha_grid) * length(w_grid) * length(dt_grid)
        run_idx = 0
        
        # Evaluate multi dimensional sweep over simulation attributes
        for a in alpha_grid
            for w in w_grid
                for dt in dt_grid
                    run_idx += 1
                    current_best_display = best_rmse == Inf ? (lowest_rmse_ever == Inf ? "N/A" : round(lowest_rmse_ever, digits=2)) : round(best_rmse, digits=2)
                    print("\r    -> Progress: $(lpad(run_idx, 3, ' '))/$(total_runs) | a: $(round(a, digits=1)), w: $(round(w, digits=1)), dt: $(round(dt, digits=1)) | Best RMSE: $(current_best_display) mV        ")
                    flush(stdout)
                    
                    reset_simulator!(cosim)
                    sweep_start = time()
                    df_sweep = simulate_pack!(cosim, exp, nema_params, tms;
                        total_time=duration, is_isothermal=true, save_csv=false, dense_logging=false, verbose=false,
                        m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive,
                        opt_alpha_elec=Dict(regime => a), opt_w_ratio=Dict(regime => w), opt_dt_max_elec=Dict(regime => dt)
                    )
                    sweep_dur = time() - sweep_start
                    
                    t_sparse = df_sweep.Time_s
                    v_sparse = df_sweep.Pack_Voltage_V
                    if length(t_sparse) < 5 continue end
                    
                    # Compute resultant validation metrics and update parameter assignments
                    interp = LinearInterpolation(v_sparse, t_sparse)
                    t_eval = clamp.(t_true, t_sparse[1], t_sparse[end])
                    v_pred = interp.(t_eval)
                    
                    rmse_mv = calc_rmse(v_true, v_pred) * 1000.0
                    r2 = calc_r2(v_true, v_pred)
                    speedup = truth_dur / sweep_dur
                    
                    if rmse_mv < lowest_rmse_ever
                        lowest_rmse_ever, best_fail_a, best_fail_w, best_fail_dt, best_fail_speedup, best_fail_r2 = rmse_mv, a, w, dt, speedup, r2
                    end
                    
                    if rmse_mv <= target_rmse && r2 >= target_r2 && (speedup > best_speedup || best_speedup == 0.0)
                        best_alpha, best_w, best_dt, best_speedup, best_rmse, best_r2 = a, w, dt, speedup, rmse_mv, r2
                    end
                end
            end
        end
        
        println("\n") 
        
        # Handle conditional output based on whether acceptable configuration bounds were reached
        if best_speedup == 0.0
            println(">>> [!] NO CONFIG MET TARGETS. SHOWING CLOSEST MATCH FOR [$regime] <<<")
            println("    -> Alpha:        $(round(best_fail_a, digits=2))")
            println("    -> W_Ratio:      $(round(best_fail_w, digits=2))")
            println("    -> dt_max:       $(round(best_fail_dt, digits=2))")
            println("    -> Speedup:      $(round(best_fail_speedup, digits=1))x")
            println("    -> Error (RMSE): $(round(lowest_rmse_ever, digits=2)) mV")
            
            reset_simulator!(cosim)
            df_optimal = simulate_pack!(cosim, exp, nema_params, tms;
                total_time=duration, is_isothermal=true, save_csv=false, dense_logging=false, verbose=false,
                m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive,
                opt_alpha_elec=Dict(regime => best_fail_a), opt_w_ratio=Dict(regime => best_fail_w), opt_dt_max_elec=Dict(regime => best_fail_dt)
            )
        else
            println(">>> OPTIMAL SETTINGS FOR [$regime] <<<")
            println("    -> Alpha:        $(round(best_alpha, digits=2))")
            println("    -> W_Ratio:      $(round(best_w, digits=2))")
            println("    -> dt_max:       $(round(best_dt, digits=2))")
            println("    -> Speedup:      $(round(best_speedup, digits=1))x")
            println("    -> Error (RMSE): $(round(best_rmse, digits=2)) mV")
            
            reset_simulator!(cosim)
            df_optimal = simulate_pack!(cosim, exp, nema_params, tms;
                total_time=duration, is_isothermal=true, save_csv=false, dense_logging=false, verbose=false,
                m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive,
                opt_alpha_elec=Dict(regime => best_alpha), opt_w_ratio=Dict(regime => best_w), opt_dt_max_elec=Dict(regime => best_dt)
            )
        end
        
        # Produce comparative visual charts
        p = plot(t_true ./ 3600.0, v_true, label="Ground Truth (0.1s)", lw=2, color=:black)
        plot!(p, df_optimal.Time_s ./ 3600.0, df_optimal.Pack_Voltage_V, label="Optimized", lw=1.5, color=:red, linestyle=:dash)
        title!("Regime [$regime] - Voltage Comparison")
        xlabel!("Time [Hours]")
        ylabel!("Pack Voltage [V]")
        savefig(p, "pareto_voltage_$(regime).png")
        println("=======================================================")
    end
end

Base.invokelatest(run_elec_pareto_tensor, CONFIG)