# ==============================================================================
# pareto_sweep_therm.jl
# Evaluate optimal explicit thermal time-stepping parameters via High-Speed Playback Spoofing
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

# --- CONFIGURATION ---
const HEAT_MULTIPLIER = 2.0  # Reduced to a realistic multiplier now that TMS is disabled!

const TARGET_MAX_ERR_C = Dict(:smooth => 0.1, :drive => 2.0, :grid => 1.0, :aging => 1.5)
const TARGET_DOSE_ERR_PCT = Dict(:smooth => 1.0, :drive => 1.0, :grid => 1.0, :aging => 1.0)

const ALPHA_RANGE = Dict(
    :smooth => (0.1, 20.0), 
    :drive => (0.1, 20.0), 
    :grid => (0.1, 20.0), 
    :aging => (0.5, 20.0)
)
const ALPHA_POINTS = 25

const DT_RANGE = Dict(
    :smooth => (100.0, 800.0), 
    :drive => (5.0, 50.0), 
    :grid => (20.0, 150.0), 
    :aging => (500.0, 2000.0) 
)
const DT_POINTS = 25
# ---------------------

function calc_thermal_dose(t_array, T_array)
    dose = 0.0
    for i in 1:(length(t_array)-1)
        dt = t_array[i+1] - t_array[i]
        T_avg = (T_array[i] + T_array[i+1]) / 2.0
        dose += dt * T_avg
    end
    return dose
end

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

function run_ghost_loop(cosim, total_time, dt_max_therm, alpha_therm, Q_interp, nema_params, tms_strategy, exp_func, cols, rows)
    cosim.therm_integrator = init(remake(cosim.prob_therm), QNDF(autodiff=false); reltol=1e-2, abstol=1e-3, save_everystep=false, verbose=false)
    
    history_t = Float64[0.0]
    history_T = Float64[nema_params.ambient_temperature - 273.15]
    
    current_t = 0.0
    dt = 0.1
    prev_T_max = nema_params.ambient_temperature - 273.15
    T_max_C = prev_T_max
    
    get_load_future = compile_experiment(exp_func)
    window = 0.0
    try window = BatteryToolkit.get_lookahead(tms_strategy) catch end
    
    while current_t < total_time
        dT_dt = abs(T_max_C - prev_T_max) / dt
        
        if alpha_therm == 0.0
            dt = dt_max_therm
        else
            proposed_dt = alpha_therm / (dT_dt + 1e-4)
            dt = clamp(proposed_dt, 0.1, dt_max_therm)
        end
        
        if current_t + dt > total_time dt = total_time - current_t end
        
        # Completely disable TMS to match the unmitigated ablation baseline
        target_flow = 0.0
        target_T = nema_params.inlet_temperature
        
        cosim.therm_integrator.ps[cosim.m_flow_sym] = target_flow
        cosim.therm_integrator.ps[cosim.T_inlet_sym] = target_T
        
        Q_now = Q_interp(current_t)
        
        # Anti-aliasing boundary fix
        eval_t = (current_t + dt >= total_time) ? (total_time - 1e-6) : (current_t + dt)
        Q_next = Q_interp(eval_t)
        
        dQ_dt_val = (Q_next - Q_now) / dt
        
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
        
        T_cores = [cosim.therm_integrator[getproperty(getproperty(cosim.sys_therm.therm_pack_base, Symbol("cell_$i")), :core_cap).T] for i in 1:cosim.num_cells]
        
        prev_T_max = T_max_C
        T_max_C = maximum(T_cores) - 273.15
        
        current_t = t_target
        push!(history_t, current_t)
        push!(history_T, T_max_C)
    end
    
    return history_t, history_T
end

function run_thermal_pareto_tensor()
    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2
    
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=1, wall_thickness=0.001)
    
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = 0.00, ambient_convection_coefficient = 5.0
    )
    
    rows, cols = 2, 2
    geom = build_pack_geometry(rows=rows, cols=cols, cell_pitch=0.025)
    tms = ReactiveTMS(T_high = 40.0, T_low = 37.0)

    regimes = [:smooth, :drive, :grid, :aging]
    PACK_1C = 5.0 * cols

    for regime in regimes
        println("\n=======================================================")
        println(" INITIATING THERMAL PARETO SWEEP: REGIME [$regime]")
        println("=======================================================")
        
        soc_init = (regime == :smooth) ? 0.45 : 0.65
        elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
        elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊
        
        exp = nothing
        duration = 3600.0
        
        if regime == :smooth
            exp = Experiment([TargetCurrentStep(0.25 * PACK_1C, 0.65, 4.2 * rows, duration)])
            
        elseif regime == :drive
            duration = 1800.0 # 30 mins
            wltp_path = joinpath("data", "V2G", "driving_power_wltp.csv")
            if isfile(wltp_path)
                f_wltp = CSV.File(wltp_path) |> Tables.matrix
                max_csv_power = maximum(abs.(f_wltp[:,2]))
                target_max_power = 2.25 * PACK_1C * (3.7 * rows)
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
            max_v2g_pwr = 2.0 * PACK_1C * (3.7 * rows)
            exp = Experiment([FCRStep(t_data, f_data, max_v2g_pwr, duration)])
            
        elseif regime == :aging
            duration = 24.0 * 3600.0
            exp = Experiment([
                TargetCurrentStep(-0.5 * PACK_1C, 0.20, 2.5 * rows, 2.5 * 3600.0), 
                RestStep(1.0 * 3600.0),                                     
                TargetCurrentStep(0.5 * PACK_1C, 0.65, 4.2 * rows, 2.5 * 3600.0),  
                RestStep(18.0 * 3600.0)                                     
            ])
        end
        
        println("[!] Generating Master Heat Tape (Multiplied by $(HEAT_MULTIPLIER)x)...")
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
        
        while tape_t < duration
            if tape_t + tape_dt > duration tape_dt = duration - tape_t end
            
            v_term = int_tape[sys_tape.cell.v]
            soc = clamp(int_tape[sys_tape.cell.soc], 0.0, 1.0)
            
            target_I, _, _ = get_load(tape_t, v_term * rows, soc)
            next_I, _, _ = get_load(tape_t + tape_dt, v_term * rows, soc)
            
            int_tape.ps[I_sym] = target_I / cols
            int_tape.ps[dI_dt_sym] = ((next_I - target_I) / tape_dt) / cols
            int_tape.ps[t_start_sym] = tape_t
            
            SciMLBase.add_tstop!(int_tape, tape_t + tape_dt)
            while int_tape.t < tape_t + tape_dt
                SciMLBase.step!(int_tape)
            end
            
            tape_t += tape_dt
            push!(t_tape, tape_t)
            push!(q_tape, int_tape[sys_tape.cell.Q_total] * HEAT_MULTIPLIER)
        end
        
        Q_interp = LinearInterpolation(q_tape, t_tape)
        println("    -> Tape Generated in $(round(time() - start_tape, digits=2)) seconds.")
        
        println("[!] Generating Thermal Ground Truth Benchmark (dt_max_therm = 0.1s)...")
        cosim = redirect_stdout(() -> build_pack_simulator(geom, elec_params, nema_params, rows, cols; verbose=false), devnull)
        
        start_ref = time()
        t_ref, T_ref = run_ghost_loop(cosim, duration, 0.1, 0.0, Q_interp, nema_params, tms, exp, cols, rows)
        truth_dur = time() - start_ref
        dose_ref = calc_thermal_dose(t_ref, T_ref)
        
        println("    -> Benchmark generated in $(round(truth_dur, digits=2)) seconds.")
        
        best_alpha = 0.0
        best_dt = 0.0
        best_speedup = 0.0
        best_max_err = Inf
        best_dose_err = Inf
        
        lowest_err_ever = Inf
        best_fail_a = 0.0
        best_fail_dt = 0.0
        best_fail_speedup = 0.0
        best_fail_dose = 0.0
        best_fail_T = []
        best_fail_t = []
        
        best_T = []
        best_t = []
        
        println("[!] Sweeping thermal parameter tensor...")
        
        target_max_err = TARGET_MAX_ERR_C[regime]
        target_dose_err = TARGET_DOSE_ERR_PCT[regime]
        
        alpha_grid = collect(range(ALPHA_RANGE[regime][1], ALPHA_RANGE[regime][2], length=ALPHA_POINTS))
        dt_grid = collect(range(DT_RANGE[regime][1], DT_RANGE[regime][2], length=DT_POINTS))
        
        total_runs = length(alpha_grid) * length(dt_grid)
        run_idx = 0
        
        for a in alpha_grid
            for dt in dt_grid
                run_idx += 1
                
                current_best_display = best_max_err == Inf ? (lowest_err_ever == Inf ? "N/A" : round(lowest_err_ever, digits=2)) : round(best_max_err, digits=2)
                print("\r    -> Progress: $(lpad(run_idx, 3, ' '))/$(total_runs) | alpha: $(round(a, digits=2)), dt_max: $(round(dt, digits=1)) | Best Peak Err: $(current_best_display) °C        ")
                flush(stdout)
                
                sweep_start = time()
                t_test, T_test = run_ghost_loop(cosim, duration, dt, a, Q_interp, nema_params, tms, exp, cols, rows)
                sweep_dur = time() - sweep_start
                
                interp_T = LinearInterpolation(T_test, t_test, extrapolation_right=ExtrapolationType.Constant)
                T_sync = [interp_T(t) for t in t_ref]
                
                max_err = maximum(abs.(T_sync .- T_ref))
                dose_test = calc_thermal_dose(t_test, T_test)
                dose_err = abs(dose_test - dose_ref) / dose_ref * 100.0
                speedup = truth_dur / sweep_dur
                
                if max_err < lowest_err_ever
                    lowest_err_ever = max_err
                    best_fail_a = a
                    best_fail_dt = dt
                    best_fail_speedup = speedup
                    best_fail_dose = dose_err
                    best_fail_T = T_test
                    best_fail_t = t_test
                end
                
                if max_err <= target_max_err && dose_err <= target_dose_err && (speedup > best_speedup || best_speedup == 0.0)
                    best_alpha = a
                    best_dt = dt
                    best_speedup = speedup
                    best_max_err = max_err
                    best_dose_err = dose_err
                    best_T = T_test
                    best_t = t_test
                end
            end
        end
        
        println("\n") 
        
        if best_speedup == 0.0
            println(">>> [!] NO CONFIG MET TARGETS. SHOWING CLOSEST MATCH FOR [$regime] <<<")
            println("    -> Alpha Therm:  $(round(best_fail_a, digits=2))")
            println("    -> dt_max_therm: $(round(best_fail_dt, digits=2))")
            println("    -> Speedup:      $(round(best_fail_speedup, digits=1))x")
            println("    -> Peak Error:   $(round(lowest_err_ever, digits=3)) °C")
            println("    -> Dose Error:   $(round(best_fail_dose, digits=3)) %")
            println("    -> Target Peak:  < $(target_max_err) °C")
            println("    -> Target Dose:  < $(target_dose_err) %")
            
            println("[!] Generating diagnostic plot for [$regime] closest match...")
            p = plot(t_ref ./ 3600.0, T_ref, label="Ground Truth (0.1s)", lw=2, color=:black)
            plot!(p, best_fail_t ./ 3600.0, best_fail_T, label="Optimized ($(round(best_fail_speedup, digits=1))x)", lw=1.5, color=:red, linestyle=:dash)
        else
            println(">>> OPTIMAL SETTINGS FOR [$regime] <<<")
            println("    -> Alpha Therm:  $(round(best_alpha, digits=2))")
            println("    -> dt_max_therm: $(round(best_dt, digits=2))")
            println("    -> Speedup:      $(round(best_speedup, digits=1))x")
            println("    -> Peak Error:   $(round(best_max_err, digits=3)) °C")
            println("    -> Dose Error:   $(round(best_dose_err, digits=3)) %")
            println("    -> Target Peak:  < $(target_max_err) °C")
            println("    -> Target Dose:  < $(target_dose_err) %")
            
            println("[!] Generating plot for [$regime] optimal settings...")
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

Base.invokelatest(run_thermal_pareto_tensor)