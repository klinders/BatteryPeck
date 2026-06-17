# ==============================================================================
# pareto_sweep.jl
# Evaluate optimal explicit time-stepping parameters via Pareto frontier analysis
# Contains:
# 1. calc_thermal_dose: Calculate total thermal dose via trapezoidal integration
# 2. run_pareto_sweep: Execute parameter sweep and identify optimal solver configuration using Kneedle algorithm
# ==============================================================================

using Plots
using BatteryToolkit
using DataFrames
using Statistics
using CSV
using LinearAlgebra
using DataInterpolations

"""
    calc_thermal_dose(t_array, T_array)

Calculate total thermal dose via trapezoidal integration.

Iterates through time and temperature arrays to compute accumulated thermal exposure.

# Arguments
- `t_array`: Array of time steps
- `T_array`: Array of corresponding temperature values

# Returns
- Calculated thermal dose
"""
function calc_thermal_dose(t_array, T_array)
    dose = 0.0
    for i in 1:(length(t_array)-1)
        dt = t_array[i+1] - t_array[i]
        T_avg = (T_array[i] + T_array[i+1]) / 2.0
        dose += dt * T_avg
    end
    return dose
end

"""
    run_pareto_sweep()

Execute parameter sweep and identify optimal solver configuration using Kneedle algorithm.

Generates ground truth reference simulation and compares multiple time-stepping configurations against it. Evaluates speedup, peak error, and thermal dose error to find optimal balance. Generates and optionally exports visual dashboard.

# Returns
- Nothing
"""
function run_pareto_sweep()
    # Configure electrical and thermal parameters alongside physical pack geometry
    elec_params = Chen2020()
    elec_params.Vmin = 1.0; elec_params.Vmax = 5.0

    soc_init = 0.4
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=1, wall_thickness=0.001)

    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026),
        pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500),
        casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom,
        cell_gap_thickness = 0.002,
        axial_potting_thickness = 0.005,
        casing_thickness = 0.01,
        ambient_temperature = 298.15,
        inlet_temperature = 298.15,
        mass_flow_rate = 0.01,
        ambient_convection_coefficient = 5
    )

    rows, cols = 2, 2
    geom = build_pack_geometry(rows=rows, cols=cols, cell_pitch=0.025)
    pack_1C_amps = 5.0 * cols

    # Define dynamic stress cycle and compile experiment
    stress_cycle = [
        RestStep(300.0),
        CurrentStep(-0.5 * pack_1C_amps, 1000.0),  
        RestStep(100.0),
        CurrentStep(2.0 * pack_1C_amps, 600.0),  
        CurrentStep(-2.0 * pack_1C_amps, 600.0), 
        RestStep(100.0),     
        CurrentStep(0.5 * pack_1C_amps, 1000.0), 
        RestStep(300.0),                                
    ]
    exp = Experiment(stress_cycle)
    total_t = sum([s.period for s in stress_cycle])
    tms = ReactiveTMS(T_high = 30.0, T_low = 27.0)

    # Execute ground truth reference simulation with strict time step limits
    println("\n=======================================================")
    println(">>> RUNNING GROUND TRUTH REFERENCE (dt_max=0.1s)...")
    println("=======================================================")
    cosim_ref = build_pack_simulator(geom, elec_params, nema_params, rows, cols)
    start_ref = time()
    
    df_ref = simulate_pack!(cosim_ref, exp, nema_params, tms; total_time=total_t, dt_max=0.1, alpha=0.0, save_csv=false,
                            m_active=0.03, m_passive=0.01)
    time_ref = time() - start_ref
    
    dose_ref = calc_thermal_dose(df_ref.Time_s, df_ref.Max_Temp_C)
    println("\n[!] Ground Truth completed in $(round(time_ref, digits=2)) seconds.")
    println("    -> Reference Thermal Dose: $(round(dose_ref, digits=2)) °C-s")

    # Define parameter space for explicit solver sweep
    dt_max_list = [1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    alpha_list  = [0.10, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0]

    matrix_sweep = Iterators.product(dt_max_list, alpha_list)
    total_runs = length(dt_max_list) * length(alpha_list)

    results = []
    run_count = 1

    println("\n=======================================================")
    println(">>> INITIATING 7x7 PARETO SWEEP ($total_runs Configurations)")
    println("=======================================================")

    # Iterate through parameter matrix and execute test simulations
    for (d_max, a_val) in matrix_sweep
        println("\n[Run $run_count/$total_runs] Testing: dt_max=$(d_max)s | alpha=$(a_val)")
        
        cosim = build_pack_simulator(geom, elec_params, nema_params, rows, cols; verbose=false)
        
        start_wall = time()
        df_test = simulate_pack!(cosim, exp, nema_params, tms; 
                                 total_time=total_t, dt_max=d_max, alpha=a_val, save_csv=false,
                                 m_active=0.03, m_passive=0.01)
        wall_time = time() - start_wall
        
        # Interpolate test results to match reference time steps and calculate absolute error
        interp_T = LinearInterpolation(df_test.Max_Temp_C, df_test.Time_s, extrapolation_right=ExtrapolationType.Constant)
        T_test_sync = [interp_T(t) for t in df_ref.Time_s]
        
        errors = T_test_sync .- df_ref.Max_Temp_C
        max_e = maximum(abs.(errors))
        speedup = time_ref / wall_time
        
        dose_test = calc_thermal_dose(df_test.Time_s, df_test.Max_Temp_C)
        dose_err_pct = abs(dose_test - dose_ref) / dose_ref * 100.0
        
        push!(results, (
            dt_max=d_max, alpha=a_val, 
            wall_time=wall_time, speedup=speedup, 
            dose_err_pct=dose_err_pct, max_err=max_e,
            raw_time=copy(df_test.Time_s), 
            raw_temp=copy(df_test.Max_Temp_C),
            sync_temp=T_test_sync
        ))
        
        # Force garbage collection to dump dead integrators from memory
        df_test = nothing
        cosim = nothing
        GC.gc() 

        run_count += 1
    end

    # Consolidate sweep metrics into data frame and export to disk
    sweep_df = DataFrame(results)

    csv_out = select(sweep_df, Not([:raw_time, :raw_temp, :sync_temp]))
    CSV.write("pareto_matrix_results.csv", csv_out)
    println("\n[!] Sweep Complete. Metrics saved to 'pareto_matrix_results.csv'")

    # Filter runs for physical validity and apply Kneedle algorithm to identify optimal configuration
    best_run = nothing
    valid_runs = filter(row -> row.max_err <= 3.0 && row.dose_err_pct <= 4.0, sweep_df)
    
    if !isempty(valid_runs)
        # Normalise metrics for comparative distance calculation
        s_min, s_max = minimum(valid_runs.speedup), maximum(valid_runs.speedup)
        s_range = s_max > s_min ? (s_max - s_min) : 1.0
        
        e_min, e_max = minimum(valid_runs.dose_err_pct), maximum(valid_runs.dose_err_pct)
        e_range = e_max > e_min ? (e_max - e_min) : 1.0
        
        max_kneedle_dist = -Inf
        
        # Find configuration maximising distance from diagonal chord to identify Pareto optimal point
        for row in eachrow(valid_runs)
            x_norm = (row.speedup - s_min) / s_range
            y_norm = (row.dose_err_pct - e_min) / e_range
            
            dist = x_norm - y_norm 
            
            if dist > max_kneedle_dist
                max_kneedle_dist = dist
                best_run = row
            end
        end
        println("\n🏆 KNEEDLE ALGORITHM OPTIMAL CONFIGURATION FOUND:")
    else
        println("\n[!] No configurations met the 3.0°C / 4.0% criteria. Defaulting to lowest dose error.")
        sort!(sweep_df, :dose_err_pct)
        best_run = sweep_df[1, :]
    end
    
    println("   dt_max: $(best_run.dt_max)s | alpha: $(best_run.alpha)")
    println("   Execution Time: $(round(best_run.wall_time, digits=2))s (Speedup: $(round(best_run.speedup, digits=2))x)")
    println("   Peak Abs Error: $(round(best_run.max_err, digits=4))°C")
    println("   Thermal Dose Error: $(round(best_run.dose_err_pct, digits=4))%")

    # Assemble scatter plots and transient traces into comprehensive dashboard
    # Plot thermal dose error correlations
    p1 = Plots.scatter(sweep_df.speedup, sweep_df.dose_err_pct, 
        zcolor=sweep_df.dt_max, markershape=:circle, markersize=7,
        markerstrokewidth=0.5, cmap=:viridis, colorbar_title="dt_max [s]",
        xlabel="Speedup factor [x]", ylabel="Dose error [%]",
        title="Dose error vs. dt_max", label="", ylims=(0, 1.0))
    Plots.scatter!(p1, [best_run.speedup], [best_run.dose_err_pct], markershape=:star5, markersize=14, color=:red, label="Optimal")

    p2 = Plots.scatter(sweep_df.speedup, sweep_df.dose_err_pct, 
        zcolor=sweep_df.alpha, markershape=:circle, markersize=7,
        markerstrokewidth=0.5, cmap=:plasma, colorbar_title="alpha [ΔT/step]",
        xlabel="Speedup factor [x]", ylabel="Dose error [%]",
        title="Dose error vs. alpha", label="", ylims=(0, 1.0))
    Plots.scatter!(p2, [best_run.speedup], [best_run.dose_err_pct], markershape=:star5, markersize=14, color=:red, label="Optimal")

    # Plot peak absolute error correlations
    p3 = Plots.scatter(sweep_df.speedup, sweep_df.max_err, 
        zcolor=sweep_df.dt_max, markershape=:circle, markersize=7,
        markerstrokewidth=0.5, cmap=:viridis, colorbar_title="dt_max [s]",
        xlabel="Speedup factor [x]", ylabel="Peak error [°C]",
        title="Peak error vs. dt_max", label="", ylims=(0, 3.0))
    Plots.scatter!(p3, [best_run.speedup], [best_run.max_err], markershape=:star5, markersize=14, color=:red, label="Optimal")

    p4 = Plots.scatter(sweep_df.speedup, sweep_df.max_err, 
        zcolor=sweep_df.alpha, markershape=:circle, markersize=7,
        markerstrokewidth=0.5, cmap=:plasma, colorbar_title="alpha [ΔT/step]",
        xlabel="Speedup factor [x]", ylabel="Peak error [°C]",
        title="Peak error vs. alpha", label="", ylims=(0, 3.0))
    Plots.scatter!(p4, [best_run.speedup], [best_run.max_err], markershape=:star5, markersize=14, color=:red, label="Optimal")

    # Plot transient equivalence traces comparing optimal run against ground truth
    p5 = Plots.plot(df_ref.Time_s, df_ref.Max_Temp_C, 
        label="Ground truth", color=:black, linewidth=2,
        xlabel="Time [s]", ylabel="Max. temp. [°C]", 
        title="Equivalence trace")
    Plots.plot!(p5, best_run.raw_time, best_run.raw_temp, 
        label="Optimal", color=:red, linestyle=:dash, linewidth=2)

    p6 = Plots.plot(df_ref.Time_s, abs.(best_run.sync_temp .- df_ref.Max_Temp_C), 
        label="Abs. error", color=:purple, linewidth=1, fill=(0, 0.2, :purple),
        xlabel="Time [s]", ylabel="Abs. error [°C]",
        title="Transient error trace")

    # Combine individual plots into final grid layout
    final_dashboard = Plots.plot(p1, p2, p5, p3, p4, p6, layout=(2, 3), size=(1500, 800), 
                                 left_margin=15Plots.mm, bottom_margin=8Plots.mm)
    Plots.display(final_dashboard)

    # Set toggle to automatically export vector graphics file of final plot
    export_vector_image = true

    # Automatically save vector version if toggle is true
    if export_vector_image
        output_filename = "pareto_sweep_plot.svg" 
        Plots.savefig(final_dashboard, output_filename)
        println("Vector image successfully exported to: $output_filename")
    end
end

Base.invokelatest(run_pareto_sweep)