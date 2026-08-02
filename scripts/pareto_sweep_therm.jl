# ==============================================================================
# pareto_sweep_therm.jl
# Evaluate optimal explicit thermal time-stepping parameters via High-Speed Playback Spoofing
# Contains:
# 1. calc_thermal_dose: Calculate total thermal dose via trapezoidal integration
# 2. load_fcr_data: Load frequency containment reserve data or generate synthetic profile
# 3. run_ghost_loop: Execute spoofed thermal loop using interpolated heat generation tape
# 4. run_thermal_pareto_tensor: Execute parameter sweep to identify optimal thermal solver configuration
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal
using OrdinaryDiffEq
using BatteryToolkit
using Plots
using DataFrames
using Statistics
using CSV
using DataInterpolations
using Printf

# Configuration dashboard
CONFIG = (
    # Architecture and scaling
    rows_series = 2,
    cols_parallel = 2,
    cell_capacity_ah = 5.0,
    cell_nominal_v = 3.7,
    id4_energy_wh = 62000.0,
    id4_fcr_power_w = 11000.0,
    
    # Thermal playback settings
    heat_multiplier = 1.0, 
    
    # By default both active and passive are set to 1e-6 to force purely passive pareto sweep
    # Change tms_m_active to 0.3 to run active pareto sweep
    tms_m_active  = 1e-6,
    tms_m_passive = 1e-6,
    tms_t_high    = 35.0,
    tms_t_low     = 32.0,
    tms_lookahead_s = 300.0,
    tms_anticipative_I_thresh = 15.0,
    tms_regime_map = Dict(:grid => :reactive, :drive => :reactive, :aging => :reactive, :smooth => :reactive, :rest => :reactive),

    # Pareto targets and ranges
    target_max_err  = Dict(:smooth => 0.07, :drive => 0.1, :grid => 0.001, :aging => 0.4),
    target_dose_err = Dict(:smooth => 1.0, :drive => 1.0, :grid => 1.0, :aging => 1.0),
    
    alpha_range = Dict(:smooth => (0.1, 20.0), :drive => (0.1, 20.0), :grid => (0.1, 20.0), :aging => (0.5, 20.0)),
    dt_range    = Dict(:smooth => (100.0, 1200.0), :drive => (5.0, 500.0), :grid => (20.0, 500.0), :aging => (500.0, 3000.0)),
    
    # Points per grid axis
    pts = 25 
)

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
    # Iterate through arrays to accumulate average temperature over time steps
    dose = 0.0
    for i in 1:(length(t_array)-1)
        dt = t_array[i+1] - t_array[i]
        T_avg = (T_array[i] + T_array[i+1]) / 2.0
        dose += dt * T_avg
    end
    return dose
end

"""
    load_fcr_data(duration)

Load frequency containment reserve data or generate synthetic profile.

Attempts to read RTE frequency data from disk. Falls back to synthetic sine wave generation upon failure.

# Arguments
- `duration`: Total duration of required data profile

# Returns
- Tuple containing time array and frequency array
"""
function load_fcr_data(duration)
    t_data = Float64[]
    f_data = Float64[]
    
    # Attempt to read RTE frequency data from disk and parse values
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
        
    # Generate synthetic sine wave profile if file reading fails
    catch
        t_data = collect(0.0 : 10.0 : duration)
        f_data = 50.0 .+ 0.15 .* sin.(t_data ./ 300.0) .+ 0.05 .* randn(length(t_data))
    end
    return t_data, f_data
end

"""
    run_ghost_loop(cosim, total_time, dt_max_therm, alpha_therm, Q_interp, nema_params, tms_strategy, m_active, m_passive, exp_func, cols, rows)

Execute spoofed thermal loop using interpolated heat generation tape.

Initialises thermal integrator and steps through simulation time. Evaluates thermal management strategy and updates heat parameters dynamically.

# Arguments
- `cosim`: Initialised pack simulator struct
- `total_time`: Total simulation duration
- `dt_max_therm`: Maximum allowed thermal time step
- `alpha_therm`: Time step scaling factor based on temperature derivative
- `Q_interp`: Interpolation object for generated heat tape
- `nema_params`: Thermal parameter definitions
- `tms_strategy`: Thermal management operational logic
- `m_active`: Active cooling mass flow rate
- `m_passive`: Passive resting mass flow rate
- `exp_func`: Compiled experiment load function
- `cols`: Number of parallel columns
- `rows`: Number of series rows

# Returns
- Tuple containing history of time steps and maximum temperatures
"""
function run_ghost_loop(cosim, total_time, dt_max_therm, alpha_therm, Q_interp, nema_params, tms_strategy, m_active, m_passive, exp_func, cols, rows)
    # Initialise thermal integrator and history tracking arrays
    cosim.therm_integrator = init(remake(cosim.prob_therm), QNDF(autodiff=false); reltol=1e-2, abstol=1e-3, save_everystep=false, verbose=false)
    
    history_t = Float64[0.0]
    history_T = Float64[nema_params.ambient_temperature - 273.15]
    
    current_t = 0.0
    dt = 0.1
    prev_T_max = nema_params.ambient_temperature - 273.15
    T_max_C = prev_T_max
    
    get_load_future = compile_experiment(exp_func)
    
    # Step through simulation time and calculate adaptive time step
    while current_t < total_time
        dT_dt = abs(T_max_C - prev_T_max) / dt
        
        if alpha_therm == 0.0
            dt = dt_max_therm
        else
            proposed_dt = alpha_therm / (dT_dt + 1e-4)
            dt = clamp(proposed_dt, 0.1, dt_max_therm)
        end
        
        if current_t + dt > total_time dt = total_time - current_t end
        
        # Query load function and evaluate future load magnitude
        _, _, current_regime = get_load_future(current_t, 4.0 * rows, 0.5)
        window = BatteryToolkit.get_lookahead(tms_strategy, current_regime)
        
        if window > 0.0
            future_I, _, _ = get_load_future(current_t + window, 4.0 * rows, 0.5)
            cell_future_load = abs(future_I) / cols
        else
            cell_future_load = 0.0
        end
        
        # Evaluate thermal management strategy and assign target flow parameters
        safe_m_passive = max(m_passive, 1e-6)

        if m_active == 0.0 && safe_m_passive == 1e-6
            target_flow = safe_m_passive
            target_T = nema_params.inlet_temperature
        else
            current_flow = cosim.therm_integrator.ps[cosim.m_flow_sym]
            target_flow, dynamic_T = BatteryToolkit.evaluate_tms_state(tms_strategy, current_regime, T_max_C, cell_future_load, current_flow, safe_m_passive, m_active)
            target_T = (target_flow >= (m_active * 0.99) && m_active > 0.0) ? dynamic_T : nema_params.inlet_temperature
        end
        
        cosim.therm_integrator.ps[cosim.m_flow_sym] = target_flow
        cosim.therm_integrator.ps[cosim.T_inlet_sym] = target_T
        
        # Query interpolated heat generation and calculate rate of change
        Q_now = Q_interp(current_t)
        eval_t = (current_t + dt >= total_time) ? (total_time - 1e-6) : (current_t + dt)
        Q_next = Q_interp(eval_t)
        dQ_dt_val = (Q_next - Q_now) / dt
        
        # Update cell heat generation parameters and step integrator forward
        for i in 1:cosim.num_cells
            cosim.therm_integrator.ps[cosim.Q_syms[i]] = Q_now
            cosim.therm_integrator.ps[cosim.dQ_dt_syms[i]] = dQ_dt_val
        end
        cosim.therm_integrator.ps[cosim.t_start_therm_sym] = current_t
        
        t_target = current_t + dt
        SciMLBase.add_tstop!(cosim.therm_integrator, t_target)
        while cosim.therm_integrator.t < t_target
            SciMLBase.step!(cosim.therm_integrator)
        end
        
        # Extract core temperatures and record maximum to history arrays
        T_cores = [cosim.therm_integrator[getproperty(getproperty(cosim.sys_therm.therm_pack_base, Symbol("cell_$i")), :core_cap).T] for i in 1:cosim.num_cells]
        
        prev_T_max = T_max_C
        T_max_C = maximum(T_cores) - 273.15
        
        current_t = t_target
        push!(history_t, current_t)
        push!(history_T, T_max_C)
    end
    
    return history_t, history_T
end

"""
    run_thermal_pareto_tensor(cfg)

Execute parameter sweep to identify optimal thermal solver configuration.

Generates ground truth thermal benchmark using high resolution integration. Sweeps through grid of alpha and maximum time step parameters to evaluate error and speedup. Identifies and plots optimal configuration.

# Arguments
- `cfg`: Configuration tuple containing architecture and sweep parameters

# Returns
- Nothing
"""
function run_thermal_pareto_tensor(cfg)
    # Calculate pack scale and configure thermophysical properties
    pack_ah = cfg.cols_parallel * cfg.cell_capacity_ah
    pack_wh = pack_ah * (cfg.rows_series * cfg.cell_nominal_v)
    energy_scale = pack_wh / cfg.id4_energy_wh
    
    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2
    
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=1, wall_thickness=0.001)
    
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = cfg.tms_m_passive, ambient_convection_coefficient = 5.0
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

    # Iterate through defined simulation regimes
    for regime in regimes
        println("\n=======================================================")
        println(" INITIATING THERMAL PARETO SWEEP: REGIME [$regime]")
        println("=======================================================")
        
        # Assign initial state of charge and construct experiment profile based on regime
        soc_init = (regime == :smooth) ? 0.45 : 0.65
        elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
        elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊
        
        exp = nothing
        duration = 3600.0
        
        if regime == :smooth
            exp = Experiment([TargetCurrentStep(0.25 * PACK_1C, 0.65, 4.2 * cfg.rows_series, duration)])
            
        elseif regime == :drive
            duration = 1800.0 
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
        
        # Configure initial cell system to generate master heat tape
        println("[!] Generating Master Heat Tape (Multiplied by $(cfg.heat_multiplier)x)...")
        start_tape = time()
        sys_tape = BatteryToolkit.build_isolated_cell_template(elec_params)
        I_sym = BatteryToolkit.get_sym(sys_tape, "I_app")
        dI_dt_sym = BatteryToolkit.get_sym(sys_tape, "dI_dt_app")
        t_start_sym = BatteryToolkit.get_sym(sys_tape, "t_start")
        T_sym = BatteryToolkit.get_sym(sys_tape, "T_ext")
        
        prob_tape = ODEProblem(sys_tape, [I_sym=>0.0, dI_dt_sym=>0.0, t_start_sym=>0.0, T_sym=>nema_params.ambient_temperature], (0.0, duration); warn_initialize_determined=false)
        int_tape = init(prob_tape, QNDF(autodiff=true); reltol=1e-5, abstol=1e-5, save_everystep=false, verbose=false)
        
        t_tape = Float64[0.0]
        q_tape = Float64[0.0]
        tape_t = 0.0
        tape_dt = 1.0
        get_load = compile_experiment(exp)
        
        # Iterate through time steps to populate heat tape arrays
        while tape_t < duration
            if tape_t + tape_dt > duration tape_dt = duration - tape_t end
            
            v_term = int_tape[sys_tape.cell.v]
            soc = clamp(int_tape[sys_tape.cell.soc], 0.0, 1.0)
            
            target_I, _, _ = get_load(tape_t, v_term * cfg.rows_series, soc)
            next_I, _, _ = get_load(tape_t + tape_dt, v_term * cfg.rows_series, soc)
            
            int_tape.ps[I_sym] = target_I / cfg.cols_parallel
            int_tape.ps[dI_dt_sym] = ((next_I - target_I) / tape_dt) / cfg.cols_parallel
            int_tape.ps[t_start_sym] = tape_t
            
            SciMLBase.add_tstop!(int_tape, tape_t + tape_dt)
            while int_tape.t < tape_t + tape_dt
                SciMLBase.step!(int_tape)
            end
            
            tape_t += tape_dt
            push!(t_tape, tape_t)
            push!(q_tape, int_tape[sys_tape.cell.Q_total] * cfg.heat_multiplier)
        end
        
        # Create interpolation object and execute thermal ground truth benchmark
        Q_interp = LinearInterpolation(q_tape, t_tape)
        println("    -> Tape Generated in $(round(time() - start_tape, digits=2)) seconds.")
        
        println("[!] Generating Thermal Ground Truth Benchmark (dt_max_therm = 0.1s)...")
        cosim = redirect_stdout(() -> build_pack_simulator(geom, elec_params, nema_params, cfg.rows_series, cfg.cols_parallel; verbose=false), devnull)
        
        start_ref = time()
        t_ref, T_ref = run_ghost_loop(cosim, duration, 0.1, 0.0, Q_interp, nema_params, tms, cfg.tms_m_active, cfg.tms_m_passive, exp, cfg.cols_parallel, cfg.rows_series)
        truth_dur = time() - start_ref
        dose_ref = calc_thermal_dose(t_ref, T_ref)
        
        println("    -> Benchmark generated in $(round(truth_dur, digits=2)) seconds.")
        
        best_alpha, best_dt, best_speedup, best_max_err, best_dose_err = 0.0, 0.0, 0.0, Inf, Inf
        lowest_err_ever, best_fail_a, best_fail_dt, best_fail_speedup, best_fail_dose = Inf, 0.0, 0.0, 0.0, 0.0
        best_fail_T, best_fail_t, best_T, best_t = [], [], [], []
        
        println("[!] Sweeping thermal parameter tensor...")
        
        target_max_err = cfg.target_max_err[regime]
        target_dose_err = cfg.target_dose_err[regime]
        
        # Define grid parameters and iterate through solver configurations
        alpha_grid = collect(range(cfg.alpha_range[regime][1], cfg.alpha_range[regime][2], length=cfg.pts))
        dt_grid = collect(range(cfg.dt_range[regime][1], cfg.dt_range[regime][2], length=cfg.pts))
        
        total_runs = length(alpha_grid) * length(dt_grid)
        run_idx = 0
        
        for a in alpha_grid
            for dt in dt_grid
                run_idx += 1
                current_best_display = best_max_err == Inf ? (lowest_err_ever == Inf ? "N/A" : round(lowest_err_ever, digits=2)) : round(best_max_err, digits=2)
                print("\r    -> Progress: $(lpad(run_idx, 3, ' '))/$(total_runs) | alpha: $(round(a, digits=2)), dt_max: $(round(dt, digits=1)) | Best Peak Err: $(current_best_display) °C        ")
                flush(stdout)
                
                # Execute spoofed thermal loop and evaluate performance metrics
                sweep_start = time()
                t_test, T_test = run_ghost_loop(cosim, duration, dt, a, Q_interp, nema_params, tms, cfg.tms_m_active, cfg.tms_m_passive, exp, cfg.cols_parallel, cfg.rows_series)
                sweep_dur = time() - sweep_start
                
                interp_T = LinearInterpolation(T_test, t_test, extrapolation_right=ExtrapolationType.Constant)
                T_sync = [interp_T(t) for t in t_ref]
                
                max_err = maximum(abs.(T_sync .- T_ref))
                dose_test = calc_thermal_dose(t_test, T_test)
                dose_err = abs(dose_test - dose_ref) / dose_ref * 100.0
                speedup = truth_dur / sweep_dur
                
                # Track optimal configurations meeting defined targets
                if max_err < lowest_err_ever
                    lowest_err_ever, best_fail_a, best_fail_dt, best_fail_speedup, best_fail_dose = max_err, a, dt, speedup, dose_err
                    best_fail_T, best_fail_t = T_test, t_test
                end
                
                if max_err <= target_max_err && dose_err <= target_dose_err && (speedup > best_speedup || best_speedup == 0.0)
                    best_alpha, best_dt, best_speedup, best_max_err, best_dose_err = a, dt, speedup, max_err, dose_err
                    best_T, best_t = T_test, t_test
                end
            end
        end
        
        println("\n") 
        
        # Plot optimal result or closest match against ground truth
        if best_speedup == 0.0
            println(">>> [!] NO CONFIG MET TARGETS. SHOWING CLOSEST MATCH FOR [$regime] <<<")
            println("    -> Alpha Therm:  $(round(best_fail_a, digits=2))")
            println("    -> dt_max_therm: $(round(best_fail_dt, digits=2))")
            println("    -> Speedup:      $(round(best_fail_speedup, digits=1))x")
            println("    -> Peak Error:   $(round(lowest_err_ever, digits=3)) °C")
            p = plot(t_ref ./ 3600.0, T_ref, label="Ground Truth (0.1s)", lw=2, color=:black)
            plot!(p, best_fail_t ./ 3600.0, best_fail_T, label="Optimized ($(round(best_fail_speedup, digits=1))x)", lw=1.5, color=:red, linestyle=:dash)
        else
            println(">>> OPTIMAL SETTINGS FOR [$regime] <<<")
            println("    -> Alpha Therm:  $(round(best_alpha, digits=2))")
            println("    -> dt_max_therm: $(round(best_dt, digits=2))")
            println("    -> Speedup:      $(round(best_speedup, digits=1))x")
            println("    -> Peak Error:   $(round(best_max_err, digits=3)) °C")
            p = plot(t_ref ./ 3600.0, T_ref, label="Ground Truth (0.1s)", lw=2, color=:black)
            plot!(p, best_t ./ 3600.0, best_T, label="Optimized ($(round(best_speedup, digits=1))x)", lw=1.5, color=:red, linestyle=:dash)
        end
        
        title!("Regime [$regime] - Temperature Comparison")
        xlabel!("Time [Hours]")
        ylabel!("Max Core Temperature [°C]")
        savefig(p, "pareto_temperature_$(regime).png")
        println("=======================================================")
    end
end

Base.invokelatest(run_thermal_pareto_tensor, CONFIG)