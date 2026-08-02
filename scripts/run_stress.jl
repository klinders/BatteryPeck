# ==============================================================================
# run_stress.jl
# 1-Week V2G Stress Test isolating Reactive and Anticipative TMS logic
# Runs at 30C Ambient with 10x scaled FCR profiles.
# Contains:
# 1. load_fcr_data: Retrieve grid frequency profile
# 2. calculate_pumping_energy: Calculate total mechanical energy consumed by cooling pump
# 3. run_stress_scenario: Execute individual stress test configuration
# 4. run_stress_matrix: Orchestrate complete stress testing matrix
# ==============================================================================

using Plots
using Revise
using BatteryToolkit
using CSV, DataFrames
using Serialization
using SciMLBase
using Printf

Revise.revise()

# Define global configuration parameters and sequential cases
CONFIG = (
    target_days = 7,             
    stress_multiplier = 10.0,    
    ambient_temp = 273.15 + 30.0,
    
    tms_cases = [:reactive, :anticipative], 
    
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
    
    tms_anticipative_I_thresh = 2.0, 
)

"""
    load_fcr_data(duration)

Retrieve grid frequency profile.

# Arguments
- `duration`: Target playback duration

# Returns
- Tuple containing temporal array and coupled frequency signal
"""
function load_fcr_data(duration)
    # Attempt to read RTE frequency data from disk and map missing values
    t_data = Float64[]
    f_data = Float64[]
    
    try
        df = CSV.read(joinpath("data", "V2G", "RTE_Frequence_2024", "RTE_Frequence_2024_02.txt"), DataFrame; delim=';')
        freq_col = names(df)[occursin.(r"freq"i, names(df))][1]
        raw_f = df[!, freq_col]
        
        f_data = map(raw_f) do val
            if ismissing(val) return 50.0 
            elseif val isa Number return Float64(val)
            elseif val isa AbstractString
                parsed = tryparse(Float64, replace(strip(val), "," => "."))
                return parsed === nothing ? 50.0 : parsed
            else return 50.0 end
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
    calculate_pumping_energy(df, cfg, nema_params)

Calculate total mechanical energy consumed by cooling pump.

# Arguments
- `df`: Output simulation data frame containing velocity array
- `cfg`: Named tuple containing global pack scaling options
- `nema_params`: Master dictionary containing geometrical channel limits

# Returns
- Calculated pumping energy in watt-hours
"""
function calculate_pumping_energy(df, cfg, nema_params)
    # Extract fluid constants and hydraulic channel geometry to compute pressure drop parameters
    rho = nema_params.fluid.density
    mu = nema_params.fluid.dynamic_viscosity
    
    W = nema_params.tms_geometry.channel_width
    H = nema_params.tms_geometry.channel_height
    N = nema_params.tms_geometry.number_of_channels
    
    A_flow = W * H * N
    P_wet = 2 * (W + H) * N
    D_h = 4 * A_flow / P_wet
    a_hyd = D_h / 2.0
    R_c = 0.0125
    
    L = cfg.rows_series * 0.025 
    
    total_energy_j = 0.0
    times = df.Time_s
    vels = df.Velocity_ms
    
    # Integrate pumping power across sequential time steps
    for i in 1:(length(times)-1)
        v = max(vels[i], 1e-6)
        dt = times[i+1] - times[i]
        
        Re = (rho * v * D_h) / mu
        De = Re * sqrt(a_hyd / R_c)
        Re_crit = 2100.0 * (1.0 + 12.0 * (R_c / a_hyd)^-0.5)
        
        f_curved_ratio = De < 30.0 ? 1.0 : (De < 300.0 ? 0.419 * De^0.275 : 0.1125 * sqrt(De))
        
        if Re < Re_crit
            f_major = (1.0 / Re) * f_curved_ratio
        else
            f_major = 4.0 * (sqrt(a_hyd / R_c) * (0.00725 + 0.076 * (Re * (a_hyd / R_c)^2)^-0.25))
        end
        
        dp = (f_major * (L / D_h) * (rho * v^2 / 2.0)) + ((42.0 * L) * (rho * v^2 / 2.0))
        
        power_w = dp * A_flow * v
        total_energy_j += power_w * dt
    end
    
    return total_energy_j / 3600.0 
end

"""
    run_stress_scenario(mode::Symbol, cfg, master_dir::String)

Execute individual stress test configuration.

# Arguments
- `mode::Symbol`: Defines target thermal management logic route
- `cfg`: Global configuration dashboard parameters
- `master_dir::String`: Assigned output directory string

# Returns
- Nothing
"""
function run_stress_scenario(mode::Symbol, cfg, master_dir::String)
    println("\n=======================================================")
    println(">>> INITIATING 1-WEEK STRESS TEST | MODE: $(uppercase(string(mode)))")
    println("=======================================================")

    start_time_real = time()

    # Define base parameters and initiate nominal state of charge
    elec_params = Chen2020()
    elec_params.Vmin = 2.5
    elec_params.Vmax = 4.2
    
    soc_init = 0.65 
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    pack_ah = cfg.cols_parallel * cfg.cell_capacity_ah
    pack_wh = pack_ah * (cfg.rows_series * cfg.cell_nominal_v)
    energy_scale = pack_wh / cfg.id4_energy_wh
    
    # Configure thermal geometries and instantiate logic components based on mode target
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols_parallel, wall_thickness=0.001)
    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = cfg.ambient_temp, inlet_temperature = cfg.ambient_temp, mass_flow_rate = cfg.tms_m_passive, ambient_convection_coefficient = 5.0 
    )
    
    m_act = cfg.tms_m_active
    m_pass = cfg.tms_m_passive
    
    if mode == :reactive
        regime_map = Dict(:grid => :reactive, :drive => :reactive, :aging => :reactive, :smooth => :reactive, :rest => :reactive)
    elseif mode == :anticipative
        regime_map = Dict(:grid => :anticipative, :drive => :reactive, :aging => :reactive, :smooth => :reactive, :rest => :reactive)
    else
        error("Unknown TMS mode.")
    end

    reactive_strat = ReactiveTMS(ambient_temp=cfg.ambient_temp, T_high=cfg.tms_t_high, T_low=cfg.tms_t_low)
    anticipative_strat = AnticipativeTMS(ambient_temp=cfg.ambient_temp, T_high=cfg.tms_t_high, T_low=cfg.tms_t_low, cell_load_threshold=cfg.tms_anticipative_I_thresh, lookahead_window=cfg.tms_lookahead_s)
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=regime_map)
    
    geom = build_pack_geometry(rows=cfg.rows_series, cols=cfg.cols_parallel, cell_pitch=0.025)
    cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows_series, cfg.cols_parallel; verbose=false, geom_sigma=0.005)
    
    pack_1C_amps = cfg.cell_capacity_ah * cfg.cols_parallel 
    pack_v_max = 4.2 * cfg.rows_series

    # Compose experimental driving and grid service schedule combining parsed real world data
    t_fcr, f_fcr = load_fcr_data(13.5 * 3600.0)
    max_v2g_power_w = cfg.id4_fcr_power_w * energy_scale * cfg.stress_multiplier
    
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

    fcr_step = FCRStep(t_fcr, f_fcr, max_v2g_power_w, 13.5 * 3600.0)

    v2g_daily_schedule = Step[
        wltp_step,                                                              
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 8.5 * 3600.0, :smooth), 
        wltp_step,                                                              
        fcr_step,                  
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 1.0 * 3600.0, :smooth)  
    ]
    exp_stress_day = Experiment(v2g_daily_schedule)
    total_t_stress_day = sum([s.period for s in v2g_daily_schedule])

    # Initiate sequentially appending daily runs handling explicit integration steps
    mode_dir = joinpath(master_dir, "mode_$(mode)")
    mkpath(mode_dir)
    master_csv_path = joinpath(mode_dir, "stress_master.csv")
    
    day = 1
    master_df = DataFrame()
    
    while day <= cfg.target_days
        prefix = "[$(uppercase(string(mode)))] Day $day/$(cfg.target_days) "
        try
            df_day = simulate_pack!(cosim, exp_stress_day, nema_params, tms; 
                total_time=total_t_stress_day, save_csv=false, 
                derate_exponent=1.0, print_prefix=prefix,
                force_even_current=false, is_isothermal=false, geom_sigma=0.005,
                m_active=m_act, m_passive=m_pass, heat_multiplier=1.0,
                
                opt_w_ratio=Dict{Symbol, Float64}(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
                opt_alpha_elec=Dict{Symbol, Float64}(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
                opt_dt_max_elec=Dict{Symbol, Float64}(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
                
                opt_alpha_therm=Dict{Symbol, Float64}(:smooth => 0.1, :drive => 12.0, :grid => 12.0, :aging => 0.5, :rest => 0.1),
                opt_dt_max_therm=Dict{Symbol, Float64}(:smooth => 558.3, :drive => 150.0, :grid => 150.0, :aging => 2062.5, :rest => 3600.0),
                
                dense_logging=false, sparse_logging=true, ultra_sparse_logging=false, verbose=true
            )
            
            append!(master_df, df_day)
            CSV.write(master_csv_path, df_day, append=isfile(master_csv_path))
            day += 1
        catch e
            println("\n[!] FATAL CRASH: Pack failed under stress on Day $day. Exception: ", e)
            break
        end
    end
    
    # Calculate pumping energy usage over duration and force garbage collection cleanup
    duration_minutes = round((time() - start_time_real) / 60.0, digits=1)
    energy_wh = calculate_pumping_energy(master_df, cfg, nema_params)
    
    println("\n[!] Scenario $(uppercase(string(mode))) completed in $duration_minutes minutes.")
    println(">>> Total TMS Pumping Energy Consumed: $(round(energy_wh, digits=4)) Wh")
    
    cosim = nothing; GC.gc() 
end

"""
    run_stress_matrix(cfg)

Orchestrate complete stress testing matrix.

# Arguments
- `cfg`: Defined global configuration settings array
"""
function run_stress_matrix(cfg)
    # Create directory sequence and iterate through defined evaluation cases
    timestamp = round(Int, time())
    master_dir = joinpath(pwd(), "results", "stress_run_$timestamp")
    mkpath(master_dir)
    println(">>> Created centralized stress test directory: $master_dir")
    
    for tms_case in cfg.tms_cases
        run_stress_scenario(tms_case, cfg, master_dir)
    end
    
    println("\n>>> FULL STRESS MATRIX COMPLETED SUCCESSFULLY!")
end

Base.invokelatest(run_stress_matrix, CONFIG)