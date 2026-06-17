# ==============================================================================
# run_experiments.jl
# Master script for ablation experiments
# Contains:
# 1. run_single_scenario: Execute individual ablation study configuration and log results
# 2. run_master_matrix: Route execution to single scenario or sequential full matrix
# ==============================================================================

using Plots

using Revise
using BatteryToolkit

Revise.revise()

# Set to 7 to run all scenarios sequentially
TARGET_SCENARIO = 7

# Physical geometry tolerance using normal distribution
BASE_GEOM_SIGMA = 0.005     # 99.7% of cells will sit within 3sigma of 1.5%
BAD_BATCH_SIGMA = 0.01      # 99.7% of cells will sit within 3sigma of 3%

PACK_ROWS_SERIES = 4    
PACK_COLS_PARALLEL = 7  

"""
    run_single_scenario(scn, base_sigma, bad_batch_sigma, target_rows, target_cols)

Execute individual ablation study configuration and log results.

# Arguments
- `scn`: Scenario identifier integer
- `base_sigma`: Standard geometric deviation for nominal defects
- `bad_batch_sigma`: Elevated geometric deviation for bad batch defects
- `target_rows`: Number of cells in series
- `target_cols`: Number of cells in parallel

# Returns
- Nothing
"""
function run_single_scenario(scn, base_sigma, bad_batch_sigma, target_rows, target_cols)
    # Print initialisation banner and configure electrical parameters
    start_time_real = time()
    println("\n=======================================================")
    println(">>> INITIATING ABLATION STUDY: SCENARIO $scn")
    println("=======================================================")

    elec_params = Chen2020()
    elec_params.Vmin = 1.0; elec_params.Vmax = 5.0
    soc_init = 0.85 # Updated SoC
    
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    # Define cooling channel geometry and master pack parameters
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=target_cols, wall_thickness=0.001)

    # Mass flow rate for active and passive is regulated down below when passing "df"
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026),
        pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500),
        casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom,
        cell_gap_thickness = 0.002,
        axial_potting_thickness = 0.005,
        casing_thickness = 0.01,
        ambient_temperature = 303.15, # 30°C
        inlet_temperature = 303.15,   # Coolant starts at ambient

        # NOTE: This parameter is structurally required but functionally bypassed.
        # Active and Passive flow rates are explicitly controlled in simulate_pack!() and "df"
        mass_flow_rate = 0.0,

        ambient_convection_coefficient = 5
    )

    tms = ReactiveTMS(T_high = 45.0, T_low = 42.0)

    # Route matrix configuration based on specified scenario identifier
    rows, cols = target_rows, target_cols
    is_iso = false
    force_even = false
    
    current_sigma = 0.0

    if scn == 1
        println("Configuration: Single Cell | Isothermal (Absolute Baseline)")
        rows, cols = 1, 1; is_iso = true; force_even = true
    elseif scn == 2
        println("Configuration: Single Cell | Thermal (Self-Heating Baseline)")
        rows, cols = 1, 1; is_iso = false; force_even = true
    elseif scn == 3
        println("Configuration: $(rows)s$(cols)p Pack | Isothermal | Coupled | $(base_sigma*100)% Defects (The Electrical Effect)")
        is_iso = true; force_even = false; current_sigma = base_sigma
    elseif scn == 4
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | DECOUPLED | $(base_sigma*100)% Defects (The Thermal Effect)")
        is_iso = false; force_even = true; current_sigma = base_sigma
    elseif scn == 5
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | Coupled | $(base_sigma*100)% Defects (The Truth Model)")
        is_iso = false; force_even = false; current_sigma = base_sigma
    elseif scn == 6
        println("Configuration: $(rows)s$(cols)p Pack | Thermal | Coupled | $(bad_batch_sigma*100)% Defects (The Bad Batch)")
        is_iso = false; force_even = false; current_sigma = bad_batch_sigma
    else
        error("Invalid Scenario ID. Choose 1 through 7.")
    end

    # Calculate pack voltage limits and construct base experiment cycle
    pack_1C_amps = 5.0 * cols 
    
    pack_v_max = 4.2 * rows
    pack_v_min = 2.5 * rows

    base_cycle = [
        RestStep(3600.0*3),
        # TargetCurrentStep(Current, Target SoC, Target Voltage, Max Duration)
        TargetCurrentStep(0.25 * pack_1C_amps, 0.1, pack_v_min, 3600*3),
        TargetCurrentStep(-1.25 * pack_1C_amps, 0.9, pack_v_max, 3600.0*0.2),     
        TargetCurrentStep(1.25 * pack_1C_amps, 0.1, pack_v_min, 3600.0*0.2),  
        TargetCurrentStep(-0.25 * pack_1C_amps, 0.9, pack_v_max, 3600*3)
    ]
    
    # Extend base cycle to reach target simulation duration
    # 91 days
    target_simulation_time = 3 * 30 * 24 * 3600 + 24 * 3600
    cycle_duration = sum(s.period for s in base_cycle)
    num_repeats = floor(Int, target_simulation_time / cycle_duration)
    stress_cycle = repeat(base_cycle, num_repeats)
    
    # Safely scale load for single cell baselines or use full pack cycle
    if rows == 1
        scaled_steps = map(stress_cycle) do s
            if s isa RestStep 
                return RestStep(s.period)
            elseif s isa TargetCurrentStep
                return TargetCurrentStep(s.value / cols, s.target_soc, s.target_v / rows, s.period)
            else 
                return CurrentStep((hasproperty(s, :I) ? s.I : s.value) / cols, s.period) 
            end
        end
        exp = Experiment(scaled_steps)
    else
        exp = Experiment(stress_cycle)
    end
    
    total_t = sum([s.period for s in stress_cycle])

    # Execute pack simulation and save output
    geom = build_pack_geometry(rows=rows, cols=cols, cell_pitch=0.025)
    
    cosim = build_pack_simulator(geom, elec_params, nema_params, rows, cols; 
                                 verbose=true, geom_sigma=current_sigma)
    
    df = simulate_pack!(cosim, exp, nema_params, tms; 
                        total_time=total_t, dt_max=25.0, alpha=1.0, 
                        save_csv=true, force_even_current=force_even, 
                        is_isothermal=is_iso, geom_sigma=current_sigma,
                        m_active=0.0, m_passive=0.0)
    
    # Log completion time and force garbage collection to manage memory
    end_time_real = time()
    duration_minutes = round((end_time_real - start_time_real) / 60.0, digits=1)
    
    println("\n[!] Scenario $scn completed and saved to CSV.")
    println(">>> Duration: $duration_minutes minutes")
    
    # Brutal Garbage Collection so Scenarios don't overflow your computer's RAM
    cosim = nothing
    df = nothing
    GC.gc() 
end

"""
    run_master_matrix(target_scn, base_sigma, bad_batch_sigma, rows, cols)

Route execution to single scenario or sequential full matrix.

# Arguments
- `target_scn`: Identifier for specific scenario or trigger for full matrix
- `base_sigma`: Standard geometric deviation for nominal defects
- `bad_batch_sigma`: Elevated geometric deviation for bad batch defects
- `rows`: Number of cells in series
- `cols`: Number of cells in parallel

# Returns
- Nothing
"""
function run_master_matrix(target_scn, base_sigma, bad_batch_sigma, rows, cols)
    # Iterate through all scenarios if target matches full matrix trigger
    if target_scn == 7
        start_matrix = time()
        println("\n=======================================================")
        println("INITIATING FULL ABLATION MATRIX (SCENARIOS 1-6)")
        println("=======================================================")
        for s in 1:6
            run_single_scenario(s, base_sigma, bad_batch_sigma, rows, cols)
        end
        end_matrix = time()
        total_dur = round((end_matrix - start_matrix) / 3600.0, digits=2)
        println("\n ALL SCENARIOS COMPLETED SUCCESSFULLY! (Total Time: $total_dur Hours)")
    # Execute specific scenario if individual target provided
    else
        run_single_scenario(target_scn, base_sigma, bad_batch_sigma, rows, cols)
    end
end

Base.invokelatest(run_master_matrix, TARGET_SCENARIO, BASE_GEOM_SIGMA, BAD_BATCH_SIGMA, PACK_ROWS_SERIES, PACK_COLS_PARALLEL)