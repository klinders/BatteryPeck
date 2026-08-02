# ==============================================================================
# run_ablation.jl
# Master script for 1-Year Multi-Fidelity Unmitigated Ablation Experiments.
# Runs the baseline aging profile with TMS completely disabled to verify drift.
# Contains:
# 1. save_checkpoint: Serialize current system state to disk
# 2. load_checkpoint!: Restore previous system state from disk
# 3. run_scenario: Execute multi fidelity ablation simulation for specific configuration
# 4. run_ablation_matrix: Orchestrate full ablation matrix across multiple discharge rates
# ==============================================================================

using Plots
using Revise
using BatteryToolkit
using CSV, DataFrames
using Serialization
using SciMLBase
using Printf

Revise.revise()

# Define global configuration dashboard parameters
CONFIG = (
    target_scenario = 7,        
    target_starting_soh = 1.0,  
    target_days = 180,          
    ablation_c_rates = [0.1, 0.75],
    
    base_geom_sigma = 0.005,    
    bad_batch_sigma = 0.01,     
    
    rows_series = 4,    
    cols_parallel = 7,  
    cell_capacity_ah = 5.0,
    cell_nominal_v = 3.7,

    tms_m_active  = 1e-6,    
    tms_m_passive = 1e-6,   
    tms_t_high    = 35.0,   
    tms_t_low     = 32.0,   
    tms_lookahead_s = 300.0, 
    tms_anticipative_I_thresh = 5.0, 
    tms_regime_map = Dict(:aging => :reactive, :rest => :reactive)
)

"""
    save_checkpoint(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)

Serialize current system state to disk.

# Arguments
- `cosim::ExplicitPackSimulator`: Active pack simulator instance
- `filepath::String`: Destination path for saved state
- `silent::Bool`: Flag to suppress terminal output
"""
function save_checkpoint(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)
    # Extract active states and parameters and serialize to defined path
    data = Dict(
        "therm_u" => cosim.therm_integrator.u, "therm_p" => cosim.therm_integrator.p, "therm_t" => cosim.therm_integrator.t,
        "cells_u" => [int.u for int in cosim.cell_integrators], "cells_p" => [int.p for int in cosim.cell_integrators], "cells_t" => [int.t for int in cosim.cell_integrators]
    )
    serialize(filepath, data)
end

"""
    load_checkpoint!(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)

Restore previous system state from disk.

# Arguments
- `cosim::ExplicitPackSimulator`: Target pack simulator instance
- `filepath::String`: Source path for saved state
- `silent::Bool`: Flag to suppress terminal output
"""
function load_checkpoint!(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)
    # Deserialize state file and reinitialise all integrators
    data = deserialize(filepath)
    SciMLBase.reinit!(cosim.therm_integrator, data["therm_u"]; t0=data["therm_t"], reset_dt=true)
    cosim.therm_integrator.p = data["therm_p"]
    for i in 1:cosim.num_cells
        SciMLBase.reinit!(cosim.cell_integrators[i], data["cells_u"][i]; t0=data["cells_t"][i], reset_dt=true)
        cosim.cell_integrators[i].p = data["cells_p"][i]
    end
end

"""
    run_scenario(scn, crate_dir, c_rate, cfg)

Execute multi fidelity ablation simulation for specific configuration.

# Arguments
- `scn`: Target scenario identifier
- `crate_dir`: Destination directory for output data
- `c_rate`: Applied discharge rate
- `cfg`: Global configuration parameters
"""
function run_scenario(scn, crate_dir, c_rate, cfg)
    start_time_real = time()
    println("\n=======================================================")
    println(">>> INITIATING 1-YEAR ABLATION: SCENARIO $scn | $(c_rate)C RATE")
    println("=======================================================")

    # Initialise base parameters and apply artificial aging penalty
    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2
    
    lost_soc = 1.0 - cfg.target_starting_soh
    soc_init = 0.65 - lost_soc 
    lost_ah = elec_params.Q₀ * lost_soc
    moles_Li_trapped = lost_ah * 3600.0 / 96485.0
    Area = elec_params.Hcc * elec_params.Wcc * elec_params.n_el
    Volume_n = Area * elec_params.e.Lₙ
    moles_sei = moles_Li_trapped / elec_params.n.side_reactions[1].z
    delta_sei_thickness = (moles_sei * elec_params.n.side_reactions[1].V̄) / (Volume_n * elec_params.n.aₖ)
    elec_params.n.side_reactions[1].Lf₀ += delta_sei_thickness
    println("[!] Time Machine Engaged: Cell aged to $(cfg.target_starting_soh*100)%. Added $(round(delta_sei_thickness*1e6, digits=3)) μm of SEI.")
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    # Configure pack layout and assign variance logic based on requested scenario
    rows, cols = cfg.rows_series, cfg.cols_parallel
    is_iso = false; force_even = false; current_sigma = 0.0

    if scn == 1
        println("Configuration: Single Cell | Isothermal (Absolute Baseline)")
        rows, cols = 1, 1; is_iso = true; force_even = true
    elseif scn == 2
        println("Configuration: Single Cell | Thermal (Self-Heating Baseline)")
        rows, cols = 1, 1; is_iso = false; force_even = true
    elseif scn == 3
        println("Configuration: $(rows)s$(cols)p Pack | Isothermal | Coupled | $(cfg.base_geom_sigma*100)% Defects")
        is_iso = true; force_even = false; current_sigma = cfg.base_geom_sigma
    elseif scn == 4
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | DECOUPLED | $(cfg.base_geom_sigma*100)% Defects")
        is_iso = false; force_even = true; current_sigma = cfg.base_geom_sigma
    elseif scn == 5
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | Coupled | $(cfg.base_geom_sigma*100)% Defects (The Truth Model)")
        is_iso = false; force_even = false; current_sigma = cfg.base_geom_sigma
    elseif scn == 6
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | Coupled | $(cfg.bad_batch_sigma*100)% Defects (The Bad Batch)")
        is_iso = false; force_even = false; current_sigma = cfg.bad_batch_sigma
    else
        error("Invalid Scenario ID.")
    end

    # Build underlying geometry and instantiate explicit orchestrator
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cols, wall_thickness=0.001)
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = cfg.tms_m_passive, ambient_convection_coefficient = 5.0 
    )
    
    reactive_strat = ReactiveTMS(ambient_temp=298.15, T_high = cfg.tms_t_high, T_low = cfg.tms_t_low)
    anticipative_strat = AnticipativeTMS(ambient_temp=298.15, T_high = cfg.tms_t_high, T_low = cfg.tms_t_low, cell_load_threshold = cfg.tms_anticipative_I_thresh, lookahead_window = cfg.tms_lookahead_s)
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=cfg.tms_regime_map)
    
    geom = build_pack_geometry(rows=rows, cols=cols, cell_pitch=0.025)
    cosim = build_pack_simulator(geom, elec_params, nema_params, rows, cols; verbose=false, geom_sigma=current_sigma)
    
    # Define experiment load profile limits and construct step sequence
    pack_1C_amps = cfg.cell_capacity_ah * cols 
    pack_v_max = 4.2 * rows
    pack_v_min = 2.5 * rows

    charge_amps = c_rate * pack_1C_amps
    discharge_amps = -c_rate * pack_1C_amps
    
    aging_daily_schedule = Step[
        TargetCurrentStep(discharge_amps, 0.35, pack_v_min, 2.5 * 3600.0, :aging),
        RestStep(1.0 * 3600.0),
        TargetCurrentStep(charge_amps, 0.65, pack_v_max, 2.5 * 3600.0, :aging),
        RestStep(18.0 * 3600.0)
    ]
    exp_aging_day = Experiment(aging_daily_schedule)
    total_t_aging_day = sum([s.period for s in aging_daily_schedule])

    # Ensure output directories exist and preallocate checkpoint targets
    scn_dir = joinpath(crate_dir, "scn_$scn")
    mkpath(scn_dir)
    
    master_csv_path = joinpath(scn_dir, "ablation_1_year_master.csv")
    daily_save_path = joinpath(scn_dir, "ablation_working_checkpoint.jls")
    summary_csv_path = joinpath(scn_dir, "final_metrics_scn_$scn.csv")
    
    if isfile(master_csv_path) rm(master_csv_path) end
    save_checkpoint(cosim, daily_save_path; silent=true)

    day = 1
    current_derate_exp = 1.0
    
    # Iterate through target duration handling potential physical bounds crashes safely
    while day <= cfg.target_days
        daily_crashes = 0
        day_success = false
        
        while !day_success
            prefix = "[Scn $scn | $(c_rate)C] Day $(lpad(day, 3))/$(cfg.target_days) | Derate ^$(round(current_derate_exp, digits=1)) "
            
            try
                df_day = simulate_pack!(cosim, exp_aging_day, nema_params, tms; 
                    total_time=total_t_aging_day, save_csv=false, 
                    derate_exponent=current_derate_exp, print_prefix=prefix,
                    force_even_current=force_even, is_isothermal=is_iso, geom_sigma=current_sigma,
                    m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive, heat_multiplier=1.0,
                    
                    opt_w_ratio=Dict{Symbol, Float64}(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
                    opt_alpha_elec=Dict{Symbol, Float64}(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
                    opt_dt_max_elec=Dict{Symbol, Float64}(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
                    
                    opt_alpha_therm=Dict{Symbol, Float64}(:smooth => 0.1, :drive => 12.0, :grid => 12.0, :aging => 0.5, :rest => 0.1),
                    opt_dt_max_therm=Dict{Symbol, Float64}(:smooth => 558.3, :drive => 150.0, :grid => 150.0, :aging => 2062.5, :rest => 3600.0),
                    
                    dense_logging=false, sparse_logging=false, ultra_sparse_logging=true, verbose=true
                )
                
                # Check for explicit terminal codes and propagate errors if limits reached
                crashed_cell = 0
                for i in 1:cosim.num_cells
                    rc = cosim.cell_integrators[i].sol.retcode
                    if rc != SciMLBase.ReturnCode.Success && rc != SciMLBase.ReturnCode.Default
                        crashed_cell = i
                        break
                    end
                end
                
                if crashed_cell > 0 error("Cell $crashed_cell crashed.") end
                
                # Adjust output state of health logging relative to manual initial degradation penalty
                soh_offset = lost_soc * 100.0
                for i in 1:(rows*cols)
                    col = Symbol("SoH_Cell_$i")
                    if hasproperty(df_day, col) df_day[!, col] .-= soh_offset end
                end

                CSV.write(master_csv_path, df_day, append=isfile(master_csv_path))
                save_checkpoint(cosim, daily_save_path; silent=true)
                
                day_success = true
                day += 1
                
            catch e
                if e isa InterruptException
                    rethrow(e)
                end
                
                daily_crashes += 1
                print("\n[!] Crash detected on Day $day (Crash $daily_crashes/4). Exception: ", e)
                
                # Terminate gracefully upon repeated unsolvable physics bounds
                if daily_crashes >= 10
                    print("\n")
                    println("FATAL: Maximum retries (10) exceeded on Day $day.")
                    println("       Battery physical limits have collapsed. Halting Scenario $scn.")
                    day = cfg.target_days + 1 
                    break
                end
                
                print("\nReloading checkpoint...")
                load_checkpoint!(cosim, daily_save_path; silent=true)
                
                if daily_crashes % 2 == 0
                    current_derate_exp += 0.1
                    print(" Increasing Derate Exponent to $(round(current_derate_exp, digits=1)).\n")
                else
                    print(" Retrying with current exponent...\n")
                end
            end
        end
    end

    # Produce output summary files and manually trigger garbage collector
    if isfile(master_csv_path)
        df_master = CSV.read(master_csv_path, DataFrame)
        if nrow(df_master) > 0
            df_final = df_master[end:end, :] 
            CSV.write(summary_csv_path, df_final)
        end
    end
    
    duration_minutes = round((time() - start_time_real) / 60.0, digits=1)
    println("\n[!] Scenario $scn ($(c_rate)C) completed.")
    println(">>> Computation Duration: $duration_minutes minutes")
    println(">>> Master data and final summary exported to: $scn_dir")
    
    cosim = nothing; GC.gc() 
end

"""
    run_ablation_matrix(cfg)

Orchestrate full ablation matrix across multiple discharge rates.

# Arguments
- `cfg`: Global configuration parameters
"""
function run_ablation_matrix(cfg)
    # Establish parent directory and route simulation sequence based on target scenario identifier
    timestamp = round(Int, time())
    master_dir = joinpath(pwd(), "results", "ablation_run_$timestamp")
    mkpath(master_dir)
    println(">>> Created centralized ablation directory: $master_dir")
    
    for c_rate in cfg.ablation_c_rates
        crate_dir = joinpath(master_dir, "crate_$(c_rate)")
        mkpath(crate_dir)
        
        if cfg.target_scenario == 7
            start_matrix = time()
            println("\n=======================================================")
            println("INITIATING FULL ABLATION MATRIX ($(c_rate)C | SCENARIOS 1-6)")
            println("=======================================================")
            for s in 1:6 
                run_scenario(s, crate_dir, c_rate, cfg) 
            end
            total_dur = round((time() - start_matrix) / 3600.0, digits=2)
            println("\n ALL SCENARIOS FOR $(c_rate)C COMPLETED! (Time: $total_dur Hours)")
        else
            run_scenario(cfg.target_scenario, crate_dir, c_rate, cfg)
        end
    end
end

Base.invokelatest(run_ablation_matrix, CONFIG)