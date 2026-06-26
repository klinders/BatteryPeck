# ==============================================================================
# pareto_sweep_elec.jl
#
# Configurable regime pareto tensor using explicit orchestrator logic.
# Contains:
# calc_rmse: Calculate root mean square error between two arrays
# calc_r2: Calculate coefficient of determination between two arrays
# load_fcr_data: Safely load or synthesize grid frequency data
# run_elec_pareto_tensor: Execute parameter sweep and identify optimal configuration
# ==============================================================================

using BatteryToolkit
using CSV, DataFrames, DataInterpolations
using Plots
using Statistics
using Printf

# --- CONFIGURATION ---
const TARGET_RMSE_MV = Dict(:smooth => 50.0, :drive => 50.0, :grid => 50.0, :aging => 150.0)
const TARGET_R2 = Dict(:smooth => 0.95, :drive => 0.95, :grid => 0.95, :aging => 0.85)

# Tightly bracketed around previous optimal runs, adjusted to fix ringing
const ALPHA_RANGE = Dict(
    :smooth => (10.0, 30.0), 
    :drive => (40.0, 80.0), 
    :grid => (1.0, 10.0), 
    :aging => (5.0, 25.0)
)
const ALPHA_POINTS = 5

const W_RANGE = Dict(
    :smooth => (1.0, 4.0), 
    :drive => (2.0, 5.0), 
    :grid => (1.0, 4.0), 
    :aging => (1.0, 5.0)
)
const W_POINTS = 5

const DT_RANGE = Dict(
    :smooth => (500.0, 1500.0), # Lowered from 400s to prevent CV phase ringing
    :drive => (0.5, 2.0),     # Shifted down to capture high frequency WLTP spikes
    :grid => (5.0, 15.0),     # Shifted down slightly to push R2 over 0.95
    :aging => (200.0, 600.0)  
)
const DT_POINTS = 5
# ---------------------

"""
    calc_rmse(y_true, y_pred)
"""
function calc_rmse(y_true, y_pred) 
    return sqrt(mean((y_true .- y_pred).^2)) 
end

"""
    calc_r2(y_true, y_pred)
"""
function calc_r2(y_true, y_pred)
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return 1.0 - (ss_res / (ss_tot + 1e-12))
end

"""
    load_fcr_data(duration)
"""
function load_fcr_data(duration)
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
        t_data = collect(0.0 : 10.0 : duration)
        f_data = 50.0 .+ 0.15 .* sin.(t_data ./ 300.0) .+ 0.05 .* randn(length(t_data))
    end
    return t_data, f_data
end

"""
    run_elec_pareto_tensor()
"""
function run_elec_pareto_tensor()
    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2

    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=1, wall_thickness=0.001)
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = 0.00, ambient_convection_coefficient = 5
    )
    geom = build_pack_geometry(rows=1, cols=1, cell_pitch=0.025)
    tms = ReactiveTMS(T_high = 45.0, T_low = 42.0)

    regimes = [:smooth, :drive, :grid, :aging]
    PACK_1C = 5.0

    for regime in regimes
        println("\n=======================================================")
        println(" INITIATING PARETO SWEEP: REGIME [$regime]")
        println("=======================================================")
        
        # Initialize isolated cell SoC dynamically based on the regime
        soc_init = (regime == :smooth) ? 0.45 : 0.65
        elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
        elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊
        
        exp = nothing
        duration = 3600.0
        
        if regime == :smooth
            exp = Experiment([TargetCurrentStep(0.25 * PACK_1C, 0.65, 4.2, duration)])
            
        elseif regime == :drive
            wltp_path = joinpath("data", "V2G", "driving_power_wltp.csv")
            if isfile(wltp_path)
                f_wltp = CSV.File(wltp_path) |> Tables.matrix
                max_csv_power = maximum(abs.(f_wltp[:,2]))
                target_max_power = 2.25 * PACK_1C * 3.7  # Scaled exactly to your 24h script
                wltp_scaler = target_max_power / max_csv_power
                
                t_wltp = f_wltp[:,1]
                p_wltp = -f_wltp[:,2] .* wltp_scaler
                dt_wltp = diff(t_wltp)
                tend = findfirst(t_wltp .>= duration)
                if isnothing(tend) tend = length(dt_wltp) end
                exp = Experiment([DriveStep(Any[dt_wltp[1:tend], p_wltp[1:tend]], duration)])
            else
                exp = Experiment([CurrentStep(-0.5 * PACK_1C, duration)])
            end
            
        elseif regime == :grid
            t_data, f_data = load_fcr_data(duration)
            max_v2g_pwr = 2.0 * PACK_1C * 3.7 # Scaled exactly to your 24h script
            exp = Experiment([FCRStep(t_data, f_data, max_v2g_pwr, duration)])
            
        elseif regime == :aging
            duration = 24.0 * 3600.0
            exp = Experiment([
                TargetCurrentStep(-0.5 * PACK_1C, 0.20, 2.5, 2.5 * 3600.0), # Discharge to 20%
                RestStep(1.0 * 3600.0),                                     # Rest
                TargetCurrentStep(0.5 * PACK_1C, 0.65, 4.2, 2.5 * 3600.0),  # Charge to 65%
                RestStep(18.0 * 3600.0)                                     # Rest
            ])
        end
        
        cosim = build_pack_simulator(geom, elec_params, nema_params, 1, 1; verbose=false)
        
        println("[!] Generating Ground Truth Benchmark (dt = 0.1s)...")
        truth_start = time()
        df_truth = simulate_pack!(cosim, exp, nema_params, tms; 
            total_time=duration, is_isothermal=true, save_csv=false, dense_logging=true, verbose=false,
            opt_alpha_elec=Dict(regime => 0.0), opt_dt_max_elec=Dict(regime => 0.1)
        )
        truth_dur = time() - truth_start
        t_true = df_truth.Time_s
        v_true = df_truth.Pack_Voltage_V
        
        println("    -> Benchmark generated in $(round(truth_dur, digits=2)) seconds.")
        
        best_alpha = 0.0
        best_w = 0.0
        best_dt = 0.0
        best_speedup = 0.0
        best_rmse = Inf
        best_r2 = 0.0
        
        lowest_rmse_ever = Inf
        best_fail_a = 0.0
        best_fail_w = 0.0
        best_fail_dt = 0.0
        best_fail_speedup = 0.0
        best_fail_r2 = 0.0
        
        println("[!] Sweeping parameter tensor...")
        
        target_rmse = TARGET_RMSE_MV[regime]
        target_r2 = TARGET_R2[regime]
        
        alpha_grid = collect(range(ALPHA_RANGE[regime][1], ALPHA_RANGE[regime][2], length=ALPHA_POINTS))
        w_grid = collect(range(W_RANGE[regime][1], W_RANGE[regime][2], length=W_POINTS))
        dt_grid = collect(range(DT_RANGE[regime][1], DT_RANGE[regime][2], length=DT_POINTS))
        
        total_runs = length(alpha_grid) * length(w_grid) * length(dt_grid)
        run_idx = 0
        
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
                        opt_alpha_elec=Dict(regime => a), opt_w_ratio=Dict(regime => w), opt_dt_max_elec=Dict(regime => dt)
                    )
                    sweep_dur = time() - sweep_start
                    
                    t_sparse = df_sweep.Time_s
                    v_sparse = df_sweep.Pack_Voltage_V
                    if length(t_sparse) < 5 continue end
                    
                    interp = LinearInterpolation(v_sparse, t_sparse)
                    t_eval = clamp.(t_true, t_sparse[1], t_sparse[end])
                    v_pred = interp.(t_eval)
                    
                    rmse_mv = calc_rmse(v_true, v_pred) * 1000.0
                    r2 = calc_r2(v_true, v_pred)
                    speedup = truth_dur / sweep_dur
                    
                    if rmse_mv < lowest_rmse_ever
                        lowest_rmse_ever = rmse_mv
                        best_fail_a = a
                        best_fail_w = w
                        best_fail_dt = dt
                        best_fail_speedup = speedup
                        best_fail_r2 = r2
                    end
                    
                    if rmse_mv <= target_rmse && r2 >= target_r2 && (speedup > best_speedup || best_speedup == 0.0)
                        best_alpha = a
                        best_w = w
                        best_dt = dt
                        best_speedup = speedup
                        best_rmse = rmse_mv
                        best_r2 = r2
                    end
                end
            end
        end
        
        println("\n") 
        
        if best_speedup == 0.0
            println(">>> [!] NO CONFIG MET TARGETS. SHOWING CLOSEST MATCH FOR [$regime] <<<")
            println("    -> Alpha:        $(round(best_fail_a, digits=2))")
            println("    -> W_Ratio:      $(round(best_fail_w, digits=2))")
            println("    -> dt_max:       $(round(best_fail_dt, digits=2))")
            println("    -> Speedup:      $(round(best_fail_speedup, digits=1))x")
            println("    -> Error (RMSE): $(round(lowest_rmse_ever, digits=2)) mV")
            println("    -> R² Score:     $(round(best_fail_r2, digits=4))")
            println("    -> Target RMSE:  < $(target_rmse) mV")
            println("    -> Target R²:    > $(target_r2)")
            
            # Save the failure plot so we can diagnose it
            println("[!] Generating diagnostic plot for [$regime] closest match...")
            reset_simulator!(cosim)
            df_optimal = simulate_pack!(cosim, exp, nema_params, tms;
                total_time=duration, is_isothermal=true, save_csv=false, dense_logging=false, verbose=false,
                opt_alpha_elec=Dict(regime => best_fail_a), opt_w_ratio=Dict(regime => best_fail_w), opt_dt_max_elec=Dict(regime => best_fail_dt)
            )
        else
            println(">>> OPTIMAL SETTINGS FOR [$regime] <<<")
            println("    -> Alpha:        $(round(best_alpha, digits=2))")
            println("    -> W_Ratio:      $(round(best_w, digits=2))")
            println("    -> dt_max:       $(round(best_dt, digits=2))")
            println("    -> Speedup:      $(round(best_speedup, digits=1))x")
            println("    -> Error (RMSE): $(round(best_rmse, digits=2)) mV")
            println("    -> R² Score:     $(round(best_r2, digits=4))")
            println("    -> Target RMSE:  < $(target_rmse) mV")
            println("    -> Target R²:    > $(target_r2)")
            
            println("[!] Generating plot for [$regime] optimal settings...")
            reset_simulator!(cosim)
            df_optimal = simulate_pack!(cosim, exp, nema_params, tms;
                total_time=duration, is_isothermal=true, save_csv=false, dense_logging=false, verbose=false,
                opt_alpha_elec=Dict(regime => best_alpha), opt_w_ratio=Dict(regime => best_w), opt_dt_max_elec=Dict(regime => best_dt)
            )
        end
        
        p = plot(t_true ./ 3600.0, v_true, label="Ground Truth (0.1s)", lw=2, color=:black)
        plot!(p, df_optimal.Time_s ./ 3600.0, df_optimal.Pack_Voltage_V, label="Optimized", lw=1.5, color=:red, linestyle=:dash)
        title!("Regime [$regime] - Voltage Comparison")
        xlabel!("Time [Hours]")
        ylabel!("Pack Voltage [V]")
        savefig(p, "pareto_voltage_$(regime).png")
        
        println("=======================================================")
    end
end

Base.invokelatest(run_elec_pareto_tensor)