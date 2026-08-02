# ==============================================================================
# run_v2g_degradation.jl
# Dedicated script for the full 6-Year Lifecycle (5-Yr Aging + 1-Yr V2G).
# Features autonomous daily checkpointing, adaptive BMS peak shaving, 
# bulletproof on-the-fly CSV appending, and a 4-crash circuit breaker.
# Contains:
# 1. load_fcr_data: Retrieve grid frequency profile
# 2. save_checkpoint: Serialize current system state to disk
# 3. load_checkpoint!: Restore previous system state from disk
# 4. run_full_v2g_degradation: Execute complete lifecycle degradation protocol
# ==============================================================================

using Plots
using Revise
using BatteryToolkit
using CSV, DataFrames
using Serialization
using SciMLBase
using Printf

Revise.revise()

# Define structural checkpoint parameters and experiment boundaries
CONFIG = (
    start_mode = :fresh, 
    
    target_years_aging = 5,
    target_days_v2g = 365,
    
    input_5yr_checkpoint = "aging_5_years_scn_5.jls", 
    resume_folder = "results/lifecycle_run_1782566228",
    resume_file   = "v2g_working_checkpoint.jls",

    enable_fcr = false,

    rows_series = 4,
    cols_parallel = 7,
    cell_capacity_ah = 5.0,
    cell_nominal_v = 3.7,

    id4_energy_wh = 62000.0,
    id4_fcr_power_w = 11000.0,

    tms_m_active  = 0.03,    
    tms_m_passive = 1e-6,   
    tms_t_high    = 35.0,   
    tms_t_low     = 32.0,   
    tms_lookahead_s = 300.0, 
    tms_anticipative_I_thresh = 5.0, 
    
    tms_regime_map = Dict(:grid => :reactive, :drive => :reactive, :aging => :reactive, :smooth => :reactive, :rest => :reactive)
)

"""
    load_fcr_data(duration)

Retrieve grid frequency profile.

# Arguments
- `duration`: Defined length to cover array span

# Returns
- Tuple containing time array and signal profile array
"""
function load_fcr_data(duration)
    # Read grid control data handling absent measurements safely
    t_data = Float64[]
    f_data = Float64[]
    
    try
        df = CSV.read(joinpath("data", "V2G", "RTE_Frequence_2024", "RTE_Frequence_2024_02.txt"), DataFrame; delim=';')
        freq_col = names(df)[occursin.(r"freq"i, names(df))][1]
        raw_f = df[!, freq_col]
        
        f_data = map(raw_f) do val
            if ismissing(val)
                return 50.0 
            elseif val isa Number
                return Float64(val)
            elseif val isa AbstractString
                parsed = tryparse(Float64, replace(strip(val), "," => "."))
                return parsed === nothing ? 50.0 : parsed
            else
                return 50.0
            end
        end
        
        f_data = convert(Vector{Float64}, f_data)
        t_data = collect(0.0 : 10.0 : (length(f_data)-1)*10.0)
    catch e
        println("[!] Could not parse 2024 frequency data. Falling back to synthetic 50Hz noise.")
        t_data = collect(0.0 : 10.0 : duration)
        f_data = 50.0 .+ 0.15 .* sin.(t_data ./ 300.0) .+ 0.05 .* randn(length(t_data))
    end
    
    return t_data, f_data
end

"""
    save_checkpoint(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)

Serialize current system state to disk.

# Arguments
- `cosim::ExplicitPackSimulator`: Pack simulator instance evaluating equations
- `filepath::String`: Target file output directory
- `silent::Bool`: Toggles terminal logging flag
"""
function save_checkpoint(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)
    # Map active internal integrators and write to specified disk payload
    data = Dict(
        "therm_u" => cosim.therm_integrator.u,
        "therm_p" => cosim.therm_integrator.p,
        "therm_t" => cosim.therm_integrator.t,
        "cells_u" => [int.u for int in cosim.cell_integrators],
        "cells_p" => [int.p for int in cosim.cell_integrators],
        "cells_t" => [int.t for int in cosim.cell_integrators]
    )
    serialize(filepath, data)
    if !silent println(">>> Checkpoint perfectly saved at: $filepath") end
end

"""
    load_checkpoint!(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)

Restore previous system state from disk.

# Arguments
- `cosim::ExplicitPackSimulator`: Active pack simulator instance replacing variables
- `filepath::String`: Saved file source path
- `silent::Bool`: Toggles terminal logging flag
"""
function load_checkpoint!(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)
    # Parse existing payload and strictly override sequential integrators
    data = deserialize(filepath)
    SciMLBase.reinit!(cosim.therm_integrator, data["therm_u"]; t0=data["therm_t"], reset_dt=true)
    cosim.therm_integrator.p = data["therm_p"]
    
    for i in 1:cosim.num_cells
        SciMLBase.reinit!(cosim.cell_integrators[i], data["cells_u"][i]; t0=data["cells_t"][i], reset_dt=true)
        cosim.cell_integrators[i].p = data["cells_p"][i]
    end
    if !silent println(">>> Checkpoint successfully restored from: $filepath") end
end

"""
    run_full_v2g_degradation(cfg)

Execute complete lifecycle degradation protocol.

# Arguments
- `cfg`: Initial configuration structure parameters
"""
function run_full_v2g_degradation(cfg)
    println("\n=======================================================")
    println(">>> INITIATING LIFECYCLE DEGRADATION PROTOCOL")
    println("=======================================================")
    
    # Establish default parameter limits and align thermal pack dimensions
    elec_params = Chen2020()
    elec_params.Vmin = 2.5
    elec_params.Vmax = 4.2
    
    soc_init = 0.65 
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    pack_ah = cfg.cols_parallel * cfg.cell_capacity_ah
    pack_wh = pack_ah * (cfg.rows_series * cfg.cell_nominal_v)
    energy_scale = pack_wh / cfg.id4_energy_wh
    
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols_parallel, wall_thickness=0.001)
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = cfg.tms_m_passive, ambient_convection_coefficient = 5.0 
    )
    
    reactive_strat = ReactiveTMS(ambient_temp=298.15, T_high=cfg.tms_t_high, T_low=cfg.tms_t_low)
    anticipative_strat = AnticipativeTMS(ambient_temp=298.15, T_high=cfg.tms_t_high, T_low=cfg.tms_t_low, cell_load_threshold=cfg.tms_anticipative_I_thresh, lookahead_window=cfg.tms_lookahead_s)
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=cfg.tms_regime_map)
    
    geom = build_pack_geometry(rows=cfg.rows_series, cols=cfg.cols_parallel, cell_pitch=0.025)
    
    pack_1C_amps = cfg.cell_capacity_ah * cfg.cols_parallel 
    pack_v_max = 4.2 * cfg.rows_series
    pack_v_min = 2.5 * cfg.rows_series

    # Format baseline aging schedule representing default vehicle life
    charge_amps = 0.1 * pack_1C_amps
    discharge_amps = -0.1 * pack_1C_amps
    
    aging_daily_schedule = Step[
        TargetCurrentStep(discharge_amps, 0.35, pack_v_min, 2.5 * 3600.0, :aging),
        RestStep(1.0 * 3600.0),
        TargetCurrentStep(charge_amps, 0.65, pack_v_max, 2.5 * 3600.0, :aging),
        RestStep(18.0 * 3600.0)
    ]
    exp_aging_day = Experiment(aging_daily_schedule)
    total_t_aging_day = sum([s.period for s in aging_daily_schedule])

    # Assign vehicle to grid profile elements depending on defined flags
    t_fcr, f_fcr = load_fcr_data(13.5 * 3600.0)
    max_v2g_power_w = cfg.id4_fcr_power_w * energy_scale 
    
    wltp_path = joinpath("data", "V2G", "driving_power_wltp.csv")
    if isfile(wltp_path)
        f_wltp = CSV.File(wltp_path) |> Tables.matrix
        t_wltp = f_wltp[:,1]
        p_wltp = -f_wltp[:,2] .* energy_scale
        dt_wltp = diff(t_wltp)
        tend = findfirst(t_wltp .>= 1800.0)
        if isnothing(tend) tend = length(dt_wltp) end
        wltp_step = DriveStep(Any[dt_wltp[1:tend], p_wltp[1:tend]], 1800.0)
    else
        wltp_step = CurrentStep(-0.5 * pack_1C_amps, 1800.0)
    end

    fcr_or_rest = cfg.enable_fcr ? FCRStep(t_fcr, f_fcr, max_v2g_power_w, 13.5 * 3600.0) : RestStep(13.5 * 3600.0)

    v2g_daily_schedule = Step[
        wltp_step,                                                              
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 8.5 * 3600.0, :smooth), 
        wltp_step,                                                              
        fcr_or_rest,                  
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 1.0 * 3600.0, :smooth)  
    ]
    exp_v2g_day = Experiment(v2g_daily_schedule)
    total_t_v2g_day = sum([s.period for s in v2g_daily_schedule])

    # Construct main testing loop isolating continuous 5 year pre degradation phase
    master_dir = joinpath(pwd(), "results", "lifecycle_run_$(round(Int, time()))")
    
    if cfg.start_mode == :fresh || cfg.start_mode == :after_balancing
        mkpath(master_dir)
        println(">>> Phase 1: Commencing 5-Year Pre-Degradation")
        cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows_series, cfg.cols_parallel; verbose=true, geom_sigma=0.005)
        
        if cfg.start_mode == :after_balancing
            load_checkpoint!(cosim, "balanced_checkpoint.jls")
            println(">>> Starting perfectly balanced from initialization checkpoint.")
        end

        aging_csv_path = joinpath(master_dir, "aging_5_years_master.csv")
        daily_aging_save_path = joinpath(master_dir, "aging_working_checkpoint.jls")
        save_checkpoint(cosim, daily_aging_save_path; silent=true)
        
        target_days_aging = cfg.target_years_aging * 365
        day = 1
        current_derate_exp = 1.0
        
        while day <= target_days_aging
            daily_crashes = 0
            day_success = false
            
            while !day_success
                prefix = "[Pre-Aging] Year $(ceil(Int, day/365)) | Day $(lpad(day, 4))/$(target_days_aging) | Derate ^$(round(current_derate_exp, digits=1)) "
                
                try
                    df_day = simulate_pack!(cosim, exp_aging_day, nema_params, tms; 
                        total_time=total_t_aging_day, save_csv=false, 
                        derate_exponent=current_derate_exp, print_prefix=prefix,
                        force_even_current=false, is_isothermal=false, geom_sigma=0.005,
                        m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive, heat_multiplier=1.0,
                        
                        opt_w_ratio=Dict{Symbol, Float64}(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
                        opt_alpha_elec=Dict{Symbol, Float64}(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
                        opt_dt_max_elec=Dict{Symbol, Float64}(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
                        
                        opt_alpha_therm=Dict{Symbol, Float64}(:smooth => 0.1, :drive => 12.0, :grid => 12.0, :aging => 0.5, :rest => 0.1),
                        opt_dt_max_therm=Dict{Symbol, Float64}(:smooth => 558.3, :drive => 150.0, :grid => 150.0, :aging => 2062.5, :rest => 3600.0),
                        
                        dense_logging=false, sparse_logging=false, ultra_sparse_logging=true, verbose=true
                    )
                    
                    crashed_cell = 0
                    for i in 1:cosim.num_cells
                        rc = cosim.cell_integrators[i].sol.retcode
                        if rc != SciMLBase.ReturnCode.Success && rc != SciMLBase.ReturnCode.Default
                            crashed_cell = i; break
                        end
                    end
                    if crashed_cell > 0 error("Cell $crashed_cell crashed.") end
                    
                    CSV.write(aging_csv_path, df_day, append=isfile(aging_csv_path))
                    save_checkpoint(cosim, daily_aging_save_path; silent=true)
                    
                    day_success = true
                    day += 1
                    
                catch e
                    if e isa InterruptException
                        rethrow(e)
                    end
                    
                    daily_crashes += 1
                    print("\n[!] Crash detected on Day $day (Crash $daily_crashes/4). Exception: ", e)
                    
                    if daily_crashes >= 4
                        print("\n")
                        error("FATAL: Maximum retries (4) exceeded on Day $day. Battery physical limits have completely collapsed. Halting simulation.")
                    end
                    
                    print("\nReloading checkpoint...")
                    load_checkpoint!(cosim, daily_aging_save_path; silent=true)
                    
                    if daily_crashes % 2 == 0
                        current_derate_exp += 0.1
                        print(" Increasing Derate Exponent to $(round(current_derate_exp, digits=1)).\n")
                    else
                        print(" Retrying with current exponent...\n")
                    end
                end
            end
        end
        
        save_checkpoint(cosim, "aging_5_years_scn_5.jls"; silent=false)
        println(">>> 5-Year Pre-Degradation Successfully Concluded!")
    end

    # Handle final discrete vehicle to grid implementation phase tracking independent output
    if cfg.start_mode == :resume_v2g || cfg.start_mode == :fresh || cfg.start_mode == :after_balancing
        
        if cfg.start_mode == :resume_v2g
            master_dir = cfg.resume_folder
        end
        
        master_csv_path = joinpath(master_dir, "v2g_master_data.csv")
        daily_v2g_save_path = joinpath(master_dir, "v2g_working_checkpoint.jls")
        
        println("\n>>> Phase 2: Loading State for 1-Year V2G Degradation Run")
        
        if cfg.start_mode == :resume_v2g
            cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows_series, cfg.cols_parallel; verbose=true, geom_sigma=0.005)
            
            if isfile(joinpath(cfg.resume_folder, cfg.resume_file))
                load_checkpoint!(cosim, joinpath(cfg.resume_folder, cfg.resume_file); silent=false)
                current_t = cosim.therm_integrator.t
                start_day = floor(Int, current_t / total_t_v2g_day) + 1
                println(">>> Resuming V2G Phase from Day $start_day (t=$(round(current_t, digits=1))s)")
            else
                load_checkpoint!(cosim, cfg.input_5yr_checkpoint; silent=false)
                start_day = 1
                println(">>> Starting V2G Phase exactly from Year 5 Handoff.")
            end
        else
            start_day = 1
        end
        
        day = start_day
        current_derate_exp = 1.0
        
        while day <= cfg.target_days_v2g
            daily_crashes = 0
            day_success = false
            
            while !day_success
                prefix = "[V2G Year] Day $(lpad(day, 3))/$(cfg.target_days_v2g) | Derate ^$(round(current_derate_exp, digits=1)) "
                
                try
                    df_day = simulate_pack!(cosim, exp_v2g_day, nema_params, tms; 
                        total_time=total_t_v2g_day, save_csv=false, 
                        derate_exponent=current_derate_exp, print_prefix=prefix,
                        force_even_current=false, is_isothermal=false, geom_sigma=0.005,
                        m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive, heat_multiplier=1.0,
                        
                        opt_w_ratio=Dict{Symbol, Float64}(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
                        opt_alpha_elec=Dict{Symbol, Float64}(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
                        opt_dt_max_elec=Dict{Symbol, Float64}(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
                        
                        opt_alpha_therm=Dict{Symbol, Float64}(:smooth => 0.1, :drive => 12.0, :grid => 12.0, :aging => 0.5, :rest => 0.1),
                        opt_dt_max_therm=Dict{Symbol, Float64}(:smooth => 558.3, :drive => 150.0, :grid => 150.0, :aging => 2062.5, :rest => 3600.0),
                        
                        dense_logging=false, sparse_logging=false, ultra_sparse_logging=true, verbose=true
                    )
                    
                    crashed_cell = 0
                    for i in 1:cosim.num_cells
                        rc = cosim.cell_integrators[i].sol.retcode
                        if rc != SciMLBase.ReturnCode.Success && rc != SciMLBase.ReturnCode.Default
                            crashed_cell = i; break
                        end
                    end
                    if crashed_cell > 0 error("Cell $crashed_cell crashed.") end
                    
                    CSV.write(master_csv_path, df_day, append=isfile(master_csv_path))
                    save_checkpoint(cosim, daily_v2g_save_path; silent=true)
                    
                    if day % 10 == 0
                        println("\n    [✓] Day $day Checkpoint safely written to disk.")
                    end
                    
                    day_success = true
                    day += 1
                    
                catch e
                    if e isa InterruptException
                        rethrow(e)
                    end
                    
                    daily_crashes += 1
                    print("\n[!] Crash detected on Day $day (Crash $daily_crashes/4). Exception: ", e)
                    
                    if daily_crashes >= 4
                        print("\n")
                        error("FATAL: Maximum retries (4) exceeded on Day $day. Battery physical limits have completely collapsed. Halting simulation.")
                    end
                    
                    print("\nReloading checkpoint...")
                    load_checkpoint!(cosim, daily_v2g_save_path; silent=true)
                    
                    if daily_crashes % 2 == 0
                        current_derate_exp += 0.1
                        print(" Increasing Derate Exponent to $(round(current_derate_exp, digits=1)).\n")
                    else
                        print(" Retrying with current exponent...\n")
                    end
                end
            end
        end
        
        println("\n=======================================================")
        println(">>> 6-YEAR DEGRADATION SIMULATION SUCCESSFULLY COMPLETED!")
        println(">>> Master Data Exported to: $master_csv_path")
        println("=======================================================")
    end
end

Base.invokelatest(run_full_v2g_degradation, CONFIG)