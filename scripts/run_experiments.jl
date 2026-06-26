# ==============================================================================
# run_experiments.jl
# Master script for 6-Year Multi-Fidelity Ablation Experiments
# Contains:
# run_single_scenario: Execute individual ablation study configuration and log results
# run_master_matrix: Route execution to single scenario or sequential full matrix
# ==============================================================================

using Plots
using Revise
using BatteryToolkit
using CSV, DataFrames
using Serialization
using SciMLBase
using Printf

Revise.revise()

# ==============================================================================
# --- SIMULATION TIMELINE CONFIGURATION ---
# ==============================================================================
SIM_DAYS_AGING = 5 * 365  # 5 Years of Accelerated Degradation
SIM_DAYS_V2G   = 1 * 365  # 1 Year of V2G Base Cycling

# Select scenario (Scenario 5 = The Truth Model: Thermal + Coupled)
TARGET_SCENARIO = 5

# Target starting State of Health (1.0 = Pristine)
TARGET_STARTING_SOH = 1.0

# Physical geometry tolerance using normal distribution
BASE_GEOM_SIGMA = 0.005     # 99.7% of cells will sit within 3sigma of 1.5%
BAD_BATCH_SIGMA = 0.01      # 99.7% of cells will sit within 3sigma of 3%

PACK_ROWS_SERIES = 4    
PACK_COLS_PARALLEL = 7  
# ==============================================================================

# --- CHECKPOINTING SYSTEM ---
function save_checkpoint(cosim::ExplicitPackSimulator, filepath::String)
    data = Dict(
        "therm_u" => cosim.therm_integrator.u,
        "therm_p" => cosim.therm_integrator.p,
        "therm_t" => cosim.therm_integrator.t,
        "cells_u" => [int.u for int in cosim.cell_integrators],
        "cells_p" => [int.p for int in cosim.cell_integrators],
        "cells_t" => [int.t for int in cosim.cell_integrators]
    )
    serialize(filepath, data)
    println("[!] Checkpoint successfully saved to: $filepath")
end

function load_checkpoint!(cosim::ExplicitPackSimulator, filepath::String)
    data = deserialize(filepath)
    SciMLBase.reinit!(cosim.therm_integrator, data["therm_u"]; t0=data["therm_t"], reset_dt=true)
    cosim.therm_integrator.p = data["therm_p"]

    for i in 1:cosim.num_cells
        SciMLBase.reinit!(cosim.cell_integrators[i], data["cells_u"][i]; t0=data["cells_t"][i], reset_dt=true)
        cosim.cell_integrators[i].p = data["cells_p"][i]
    end
    println("[!] Checkpoint successfully loaded from: $filepath")
end

# --- STATE DIAGNOSTIC TABLE ---
function print_pack_state_table(cosim::ExplicitPackSimulator, title::String)
    sys = cosim.sys_cell_template
    println("\n=======================================================")
    println(" $title")
    println("=======================================================")
    @printf("%-6s | %-12s | %-12s | %-10s | %-10s | %-10s\n", "Cell", "Voltage [V]", "Current [A]", "SoH [%]", "SoC [%]", "Temp [°C]")
    println("-"^75)
    
    for i in 1:cosim.num_cells
        int_c = cosim.cell_integrators[i]
        v_term = int_c[sys.cell.v]
        curr = int_c.ps[cosim.I_sym]
        soh = int_c[sys.cell.SoH]
        soc = int_c[sys.cell.soc] * 100.0
        
        core_T_var = getproperty(getproperty(cosim.sys_therm.therm_pack_base, Symbol("cell_$i")), :core_cap).T
        temp = cosim.therm_integrator[core_T_var] - 273.15
        
        flag = (v_term <= 2.501 || v_term >= 4.199) ? "<<< CRITICAL" : ""
        
        @printf("Cell %-2d| %-12.4f | %-12.4f | %-10.2f | %-10.2f | %-10.2f %s\n", i, v_term, curr, soh, soc, temp, flag)
    end
    println("=======================================================\n")
end

# --- PHASE 1.5: WORKSHOP ACTIVE BALANCER ---
function balance_pack!(cosim::ExplicitPackSimulator, target_soc::Float64, max_duration_hrs::Float64)
    println("\n>>> [PHASE 1.5] INITIATING WORKSHOP ACTIVE BALANCING...")
    
    num_cells = cosim.num_cells
    dt = 60.0 # 1-minute steps for smooth balancing
    duration = max_duration_hrs * 3600.0
    t_bal = 0.0
    start_t_abs = cosim.therm_integrator.t
    
    K_prop = 2.0           # Gentle real-world proportional gain
    max_bal_current = 0.2  # 200mA max balance current
    
    sys = cosim.sys_cell_template
    socs = zeros(num_cells)
    
    # Store live current states to prevent the "Ghost Parameter" bug
    I_prev = zeros(num_cells)
    for i in 1:num_cells
        I_prev[i] = cosim.cell_integrators[i].ps[cosim.I_sym]
    end
    
    any_crashed = false

    # FIX: Loop terminates slightly before duration to prevent dt=0.0 infinite loop
    while t_bal < duration - 1e-6
        # Read Current SoCs
        for i in 1:num_cells
            socs[i] = cosim.cell_integrators[i][sys.cell.soc]
        end
        
        # Live-updating terminal "Progress Bar" every hour
        if mod(t_bal, 3600.0) == 0.0 || t_bal == 0.0
            min_soc = minimum(socs) * 100.0
            max_soc = maximum(socs) * 100.0
            spread = max_soc - min_soc
            @printf("\r    -> Time: %5.1f h / %5.1f h | Target: %.1f%% | Spread: %5.3f%% | Min: %5.2f%% | Max: %5.2f%%        ", 
                    t_bal/3600.0, max_duration_hrs, target_soc*100.0, spread, min_soc, max_soc)
            flush(stdout)
            
            if spread < 0.05 && abs(max_soc - target_soc*100.0) < 0.05
                print("\n    -> Pack fully synchronized early!\n")
                break
            end
        end
        
        if t_bal + dt > duration
            dt = duration - t_bal
        end
        
        any_crashed = false
        
        # Step Electrical Solver (Using EXACT PROVEN test_balancer.jl ramp logic!)
        for i in 1:num_cells
            int_c = cosim.cell_integrators[i]
            
            if int_c.sol.retcode != SciMLBase.ReturnCode.Success && int_c.sol.retcode != SciMLBase.ReturnCode.Default
                continue
            end
            
            err = target_soc - socs[i]
            
            # K_prop polarity fix + continuous ramp
            I_target = clamp(K_prop * err, -max_bal_current, max_bal_current)
            I_curr = I_prev[i]
            dI_dt = (I_target - I_curr) / dt
            
            int_c.ps[cosim.I_sym] = I_curr
            int_c.ps[cosim.dI_dt_sym] = dI_dt
            int_c.ps[cosim.t_start_sym] = start_t_abs + t_bal
            
            t_target = start_t_abs + t_bal + dt
            SciMLBase.add_tstop!(int_c, t_target)
            while int_c.t < t_target
                SciMLBase.step!(int_c)
                if int_c.sol.retcode != SciMLBase.ReturnCode.Success && int_c.sol.retcode != SciMLBase.ReturnCode.Default
                    any_crashed = true
                    break
                end
            end
            
            I_prev[i] = I_target
        end
        
        if any_crashed
            break
        end
        
        # Step Thermal Solver (Virtually zero heat generation during balancing)
        for i in 1:num_cells
            cosim.therm_integrator.ps[cosim.Q_syms[i]] = 0.0
            cosim.therm_integrator.ps[cosim.dQ_dt_syms[i]] = 0.0
        end
        cosim.therm_integrator.ps[cosim.t_start_therm_sym] = start_t_abs + t_bal
        
        t_target_therm = start_t_abs + t_bal + dt
        SciMLBase.add_tstop!(cosim.therm_integrator, t_target_therm)
        while cosim.therm_integrator.t < t_target_therm
            SciMLBase.step!(cosim.therm_integrator)
        end
        
        t_bal += dt
    end
    
    if any_crashed
        print("\n\n[X] BALANCER HALTED: A cell hit a physical limit (Terminated). Gracefully exiting...\n")
        print_pack_state_table(cosim, "BALANCER CRASH AUTOPSY TABLE (ALL 28 CELLS)")
    else
        print("\n>>> WORKSHOP BALANCING COMPLETE.\n")
        print_pack_state_table(cosim, "POST-BALANCING VERIFICATION TABLE (ALL 28 CELLS)")
    end
end
# ----------------------------

function load_fcr_data(duration)
    t_data = Float64[]
    f_data = Float64[]
    try
        df = CSV.read(joinpath("data", "V2G", "RTE_Frequence_2024", "RTE_Frequence_2024_02.txt"), DataFrame; delim=';')
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

function run_single_scenario(scn, base_sigma, bad_batch_sigma, target_rows, target_cols)
    start_time_real = time()
    println("\n=======================================================")
    println(">>> INITIATING 6-YEAR ABLATION STUDY: SCENARIO $scn")
    println("=======================================================")

    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2
    
    lost_soc = 1.0 - TARGET_STARTING_SOH
    soc_init = 0.65 - lost_soc 
    
    lost_ah = elec_params.Q₀ * lost_soc
    moles_Li_trapped = lost_ah * 3600.0 / 96485.0
    
    Area = elec_params.Hcc * elec_params.Wcc * elec_params.n_el
    Volume_n = Area * elec_params.e.Lₙ
    moles_sei = moles_Li_trapped / elec_params.n.side_reactions[1].z
    delta_sei_thickness = (moles_sei * elec_params.n.side_reactions[1].V̄) / (Volume_n * elec_params.n.aₖ)
    
    elec_params.n.side_reactions[1].Lf₀ += delta_sei_thickness
    println("[!] Time Machine Engaged: Cell aged to $(TARGET_STARTING_SOH*100)%. Added $(round(delta_sei_thickness*1e6, digits=3)) μm of SEI.")

    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=target_cols, wall_thickness=0.001)

    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = 0.01, ambient_convection_coefficient = 5.0 
    )
    tms = ReactiveTMS(T_high = 40.0, T_low = 37.0)

    rows, cols = target_rows, target_cols
    is_iso = false; force_even = false; current_sigma = 0.0

    if scn == 1
        println("Configuration: Single Cell | Isothermal (Absolute Baseline)")
        rows, cols = 1, 1; is_iso = true; force_even = true
    elseif scn == 2
        println("Configuration: Single Cell | Thermal (Self-Heating Baseline)")
        rows, cols = 1, 1; is_iso = false; force_even = true
    elseif scn == 3
        println("Configuration: $(rows)s$(cols)p Pack | Isothermal | Coupled | $(base_sigma*100)% Defects")
        is_iso = true; force_even = false; current_sigma = base_sigma
    elseif scn == 4
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | DECOUPLED | $(base_sigma*100)% Defects")
        is_iso = false; force_even = true; current_sigma = base_sigma
    elseif scn == 5
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | Coupled | $(base_sigma*100)% Defects (The Truth Model)")
        is_iso = false; force_even = false; current_sigma = base_sigma
    elseif scn == 6
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | Coupled | $(bad_batch_sigma*100)% Defects (The Bad Batch)")
        is_iso = false; force_even = false; current_sigma = bad_batch_sigma
    else
        error("Invalid Scenario ID.")
    end

    pack_1C_amps = 5.0 * cols 
    pack_v_max = 4.2 * rows
    pack_v_min = 2.5 * rows

    println("[!] Compiling Split Phase Schedules...")
    
    # --- PHASE 1: 5 YEARS AGING SCHEDULE ---
    aging_schedule = Step[]
    for d in 1:SIM_DAYS_AGING
        push!(aging_schedule, TargetCurrentStep(-0.5 * pack_1C_amps, 0.35, pack_v_min, 2.5 * 3600.0)) 
        push!(aging_schedule, RestStep(1.0 * 3600.0))                                                 
        push!(aging_schedule, TargetCurrentStep(0.5 * pack_1C_amps, 0.65, pack_v_max, 2.5 * 3600.0))  
        push!(aging_schedule, RestStep(18.0 * 3600.0))                                                
    end
    exp_aging = Experiment(aging_schedule)
    total_t_aging = sum([s.period for s in aging_schedule])

    # --- PHASE 2: 1 YEAR V2G SCHEDULE ---
    t_fcr, f_fcr = load_fcr_data(13.5 * 3600.0)
    
    # RESTORED FULL POWER (2.25C WLTP, 2.0C FCR)
    max_v2g_power_w = 2.0 * pack_1C_amps * (3.7 * rows)
    wltp_path = joinpath("data", "V2G", "driving_power_wltp.csv")
    if isfile(wltp_path)
        f_wltp = CSV.File(wltp_path) |> Tables.matrix
        target_max_power = 2.25 * pack_1C_amps * (3.7 * rows)
        wltp_scaler = target_max_power / maximum(abs.(f_wltp[:,2]))
        t_wltp = f_wltp[:,1]; p_wltp = -f_wltp[:,2] .* wltp_scaler; dt_wltp = diff(t_wltp)
        tend = findfirst(t_wltp .>= 1800.0)
        if isnothing(tend) tend = length(dt_wltp) end
        wltp_step = DriveStep(Any[dt_wltp[1:tend], p_wltp[1:tend]], 1800.0)
    else
        wltp_step = CurrentStep(-0.5 * pack_1C_amps, 1800.0)
    end

    v2g_schedule = Step[]
    for d in 1:SIM_DAYS_V2G
        push!(v2g_schedule, wltp_step)                                                              
        push!(v2g_schedule, TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 8.5 * 3600.0)) 
        push!(v2g_schedule, wltp_step)                                                              
        push!(v2g_schedule, FCRStep(t_fcr, f_fcr, max_v2g_power_w, 13.5 * 3600.0))                  
        push!(v2g_schedule, TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 1.0 * 3600.0)) 
    end
    exp_v2g = Experiment(v2g_schedule)
    total_t_v2g = sum([s.period for s in v2g_schedule])

    # ==========================================================
    # EXECUTE MULTI-PHASE SIMULATION
    # ==========================================================
    geom = build_pack_geometry(rows=rows, cols=cols, cell_pitch=0.025)
    cosim = build_pack_simulator(geom, elec_params, nema_params, rows, cols; verbose=true, geom_sigma=current_sigma)
    
    aging_cp = joinpath(pwd(), "aging_5_years_scn_$(scn).jls")
    
    if !isfile(aging_cp)
        println("\n>>> [PHASE 1] RUNNING 5-YEAR PRE-DEGRADATION...")
        df_aging = simulate_pack!(cosim, exp_aging, nema_params, tms; 
            total_time=total_t_aging, save_csv=false, 
            force_even_current=force_even, is_isothermal=is_iso, geom_sigma=current_sigma,
            m_active=0.01, m_passive=0.01, heat_multiplier=15.0,
            opt_w_ratio=Dict(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
            opt_alpha_elec=Dict(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
            opt_dt_max_elec=Dict(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
            opt_alpha_therm=Dict(:smooth => 0.1, :drive => 2.59, :grid => 5.9, :aging => 0.5, :rest => 0.1),
            opt_dt_max_therm=Dict(:smooth => 158.33, :drive => 6.88, :grid => 57.92, :aging => 500.0, :rest => 3600.0),
            dense_logging=false
        )
        save_checkpoint(cosim, aging_cp)
    else
        println("\n>>> [PHASE 1] RAW CHECKPOINT FOUND! Loading 5-Year state into memory...")
        load_checkpoint!(cosim, aging_cp)
    end
    
    # Phase 1.5: Workshop Balancing (Will output beautiful single-line progress!)
    balance_pack!(cosim, 0.65, 48.0) 

    println("\n>>> [PHASE 2] RUNNING 1-YEAR V2G BASE CYCLES...")
    df_v2g = simulate_pack!(cosim, exp_v2g, nema_params, tms; 
        total_time=total_t_v2g, save_csv=true, 
        force_even_current=force_even, is_isothermal=is_iso, geom_sigma=current_sigma,
        m_active=0.01, m_passive=0.01, heat_multiplier=15.0,
        opt_w_ratio=Dict(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
        opt_alpha_elec=Dict(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
        opt_dt_max_elec=Dict(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
        opt_alpha_therm=Dict(:smooth => 0.1, :drive => 2.59, :grid => 5.9, :aging => 0.5, :rest => 0.1),
        opt_dt_max_therm=Dict(:smooth => 158.33, :drive => 6.88, :grid => 57.92, :aging => 500.0, :rest => 3600.0),
        dense_logging=false
    )

    # ==========================================================
    # PHASE 2 DIAGNOSTIC CHECK
    # ==========================================================
    any_v2g_crash = false
    for i in 1:cosim.num_cells
        if cosim.cell_integrators[i].sol.retcode != SciMLBase.ReturnCode.Success && cosim.cell_integrators[i].sol.retcode != SciMLBase.ReturnCode.Default
            any_v2g_crash = true
            break
        end
    end

    if any_v2g_crash
        println("\n[X] V2G PHASE HALTED: A cell hit a physical limit (Terminated).")
        print_pack_state_table(cosim, "V2G CRASH AUTOPSY TABLE (ALL 28 CELLS)")
    else
        println("\n>>> V2G PHASE SUCCESSFULLY COMPLETED 1 FULL YEAR.")
        print_pack_state_table(cosim, "END OF 1-YEAR V2G: FINAL STATE TABLE (ALL 28 CELLS)")
    end
    # ==========================================================

    # Correct Time Machine SoH Offset
    soh_offset = lost_soc * 100.0
    for i in 1:(rows*cols)
        col = Symbol("SoH_Cell_$i")
        if hasproperty(df_v2g, col) df_v2g[!, col] .-= soh_offset end
    end
    
    run_dirs = filter(isdir, readdir(pwd()))
    results_dirs = filter(d -> startswith(d, "results_run_"), run_dirs)
    if !isempty(results_dirs)
        latest_dir = sort(results_dirs)[end]
        CSV.write(joinpath(latest_dir, "master_results.csv"), df_v2g)
    end

    duration_minutes = round((time() - start_time_real) / 60.0, digits=1)
    println("\n[!] Scenario $scn completed.")
    println(">>> Computation Duration: $duration_minutes minutes")
    
    cosim = nothing; GC.gc() 
end

function run_master_matrix(target_scn, base_sigma, bad_batch_sigma, rows, cols)
    if target_scn == 7
        start_matrix = time()
        println("\n=======================================================")
        println("INITIATING FULL ABLATION MATRIX (SCENARIOS 1-6)")
        println("=======================================================")
        for s in 1:6 run_single_scenario(s, base_sigma, bad_batch_sigma, rows, cols) end
        total_dur = round((time() - start_matrix) / 3600.0, digits=2)
        println("\n ALL SCENARIOS COMPLETED SUCCESSFULLY! (Total Time: $total_dur Hours)")
    else
        run_single_scenario(target_scn, base_sigma, bad_batch_sigma, rows, cols)
    end
end

run_master_matrix(TARGET_SCENARIO, BASE_GEOM_SIGMA, BAD_BATCH_SIGMA, PACK_ROWS_SERIES, PACK_COLS_PARALLEL)