# ==============================================================================
# test_24h_cycle.jl
# Integration test for Multi-Fidelity Handoffs (5 Days Aging + 1 Day V2G)
# Contains:
# 1. load_fcr_data: Retrieve grid frequency profile
# 2. run_test_cycle: Execute multi fidelity handoff integration test
# ==============================================================================

using BatteryToolkit
using CSV, DataFrames
using Plots

# Define parameter ranges and evaluation map
CONFIG = (
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
    # Target raw data table and return default harmonic oscillation array upon missing constraints
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
                clean_val = replace(strip(val), "," => ".")
                parsed = tryparse(Float64, clean_val)
                return isnothing(parsed) ? 50.0 : parsed
            else
                return 50.0
            end
        end
        
        f_data = convert(Vector{Float64}, f_data)
        t_data = collect(0.0 : 10.0 : (length(f_data)-1)*10.0)
        
    catch e
        println("[!] Warning: Could not properly load RTE FCR data. Generating synthetic grid frequency...")
        t_data = collect(0.0 : 10.0 : duration)
        f_data = 50.0 .+ 0.15 .* sin.(t_data ./ 300.0) .+ 0.05 .* randn(length(t_data))
    end
    return t_data, f_data
end

"""
    run_test_cycle(cfg)

Execute multi fidelity handoff integration test.

# Arguments
- `cfg`: Mapped array containing global targets
"""
function run_test_cycle(cfg)
    println(">>> Compiling $(cfg.rows_series)s$(cfg.cols_parallel)p Pack for 6-Day Integration Test...")
    
    # Establish operational constants and dictate geometric arrangement constraints
    pack_ah = cfg.cols_parallel * cfg.cell_capacity_ah
    pack_wh = pack_ah * (cfg.rows_series * cfg.cell_nominal_v)
    energy_scale = pack_wh / cfg.id4_energy_wh

    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2

    soc_init = 0.65
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols_parallel, wall_thickness=0.001)

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
    cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows_series, cfg.cols_parallel; verbose=false)

    pack_1C_amps = pack_ah
    pack_v_min = 2.5 * cfg.rows_series
    pack_v_max = 4.2 * cfg.rows_series

    # Merge dynamic target profiles evaluating specific experimental duration bounds
    println(">>> Constructing 6-Day Load Profile...")
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
        println("[!] Could not find WLTP file, using fallback step.")
        wltp_step = CurrentStep(-0.5 * pack_1C_amps, 1800.0) 
    end

    schedule = Step[]

    for d in 1:5
        push!(schedule, TargetCurrentStep(-0.5 * pack_1C_amps, 0.20, pack_v_min, 2.5 * 3600.0)) 
        push!(schedule, RestStep(1.0 * 3600.0))                                                 
        push!(schedule, TargetCurrentStep(0.5 * pack_1C_amps, 0.65, pack_v_max, 2.5 * 3600.0))  
        push!(schedule, RestStep(18.0 * 3600.0))                                                
    end

    push!(schedule, wltp_step)                                                                  
    push!(schedule, TargetCurrentStep(0.25 * pack_1C_amps, 0.65, pack_v_max, 8.5 * 3600.0))     
    push!(schedule, wltp_step)                                                                  
    push!(schedule, FCRStep(t_fcr, f_fcr, max_v2g_power_w, 13.5 * 3600.0))                      
    push!(schedule, TargetCurrentStep(0.25 * pack_1C_amps, 0.65, pack_v_max, 1.0 * 3600.0))     

    exp_protocol = Experiment(schedule)
    total_duration = 6.0 * 24.0 * 3600.0

    # Start multi day simulation handling internal limits securely
    println(">>> Executing 6-Day Integration Run...")
    df = simulate_pack!(cosim, exp_protocol, nema_params, tms; 
        total_time=total_duration, save_csv=false, dense_logging=false, verbose=true, heat_multiplier=1.0,
        m_active=cfg.tms_m_active, m_passive=cfg.tms_m_passive,
        opt_w_ratio=Dict(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
        opt_alpha_elec=Dict(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
        opt_dt_max_elec=Dict(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
        opt_alpha_therm=Dict(:smooth => 0.1, :drive => 2.59, :grid => 5.9, :aging => 0.5, :rest => 0.1),
        opt_dt_max_therm=Dict(:smooth => 158.33, :drive => 6.88, :grid => 57.92, :aging => 500.0, :rest => 3600.0)
    )

    # Plot continuous dynamic properties recording bounds accurately
    println("\n>>> Generating 4-Tier Validation Dashboard...")
    t_hrs = df.Time_s ./ 3600.0
    
    T_cols = [Symbol("T_Core_Cell_$i") for i in 1:(cfg.rows_series*cfg.cols_parallel)]
    T_matrix = Matrix(df[!, T_cols])
    T_max_array = maximum(T_matrix, dims=2)[:, 1]
    T_min_array = minimum(T_matrix, dims=2)[:, 1]

    p1 = plot(t_hrs, df.Pack_Voltage_V, title="Pack Voltage (6-Day Multi-Fidelity Handoff)", ylabel="Voltage [V]", color=:blue, lw=1.5, legend=false)
    hline!(p1, [pack_v_max], color=:red, linestyle=:dash, lw=2)
    hline!(p1, [pack_v_min], color=:red, linestyle=:dash, lw=2)

    p2 = plot(t_hrs, df.Pack_Current_A, title="Pack Current", ylabel="Current [A]", color=:orange, lw=1.5, legend=false)
    hline!(p2, [0.0], color=:black, lw=1)

    p3 = plot(t_hrs, df.Pack_SoC, title="True Chemical Pack SoC (%)", ylabel="SoC [%]", color=:green, lw=2, legend=false)
    hline!(p3, [65.0], color=:black, linestyle=:dot, lw=1.5)
    hline!(p3, [20.0], color=:black, linestyle=:dot, lw=1.5)

    p4 = plot(t_hrs, T_max_array, title="Cell Temperatures (1x Heat Multiplier)", xlabel="Time [Hours]", ylabel="Temperature [°C]", color=:red, lw=1.5, label="Hottest Cell")
    plot!(p4, t_hrs, T_min_array, color=:blue, lw=1.5, label="Coldest Cell", legend=:topleft)

    display(plot(p1, p2, p3, p4, layout=(4,1), size=(1000, 1200), margin=5 * Plots.Measures.mm))
    println(">>> Test Complete! Check the plots to verify regime shifting and thermal bounds.")
end

Base.invokelatest(run_test_cycle, CONFIG)