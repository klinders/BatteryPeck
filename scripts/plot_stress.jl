# ==============================================================================
# plot_stress.jl
# Standalone script to generate and plot the 1x vs 10x FCR Stress Test
# Proves that even at 10x scale in 30°C ambient, FCR does not require active cooling.
# Contains:
# 1. load_fcr_data: Retrieve grid frequency profile
# 2. generate_dense_24h_stress_csv: Autonomously run 1-day dense simulation for 30°C stress test
# 3. plot_stress_dashboard: Render 3-panel dense comparison dashboard for 1x vs 10x FCR stress test
# ==============================================================================

using Revise
using CSV
using DataFrames
using Plots
using Statistics
using Printf

using BatteryToolkit
using ModelingToolkit
using OrdinaryDiffEq

Revise.revise()

# Create output folder autonomously to store generated traces, enforce 30°C summer conditions, set internal pack constants, and define standard MATLAB colour palette
CONFIG = (
    output_folder = "results/stress_test_dashboard",
    
    ambient_temp = 273.15 + 30.0, 
    
    rows = 4,
    cols = 7,
    cell_cap_ah = 5.0,
    cell_nominal_v = 3.7,
    id4_energy_wh = 62000.0,
    id4_fcr_power_w = 11000.0
)

MATLAB_COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]

"""
    load_fcr_data(duration)

Retrieve grid frequency profile.

# Arguments
- `duration`: Target playback duration

# Returns
- Tuple containing temporal array and coupled frequency signal
"""
function load_fcr_data(duration)
    # Initialise arrays and attempt file read to parse RTE dataset
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
    generate_dense_24h_stress_csv(cfg, out_path::String, stress_multiplier::Float64)

Autonomously run 1-day dense simulation for 30°C stress test.

# Arguments
- `cfg`: Global configuration dashboard parameters
- `out_path::String`: Target file output directory
- `stress_multiplier::Float64`: Multiplier for V2G power

# Returns
- Nothing
"""
function generate_dense_24h_stress_csv(cfg, out_path::String, stress_multiplier::Float64)
    println("\n[!] Trace not found. Generating now for $(Int(stress_multiplier))x FCR (High fidelity)...")
    
    # Establish base electrochemical parameters and instantiate initial state of charge
    elec_params = Chen2020()
    elec_params.Vmin = 2.5
    elec_params.Vmax = 4.2
    
    soc_init = 0.65 
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    # Define experimental load profile limits and construct step sequence
    pack_1C_amps = cfg.cell_cap_ah * cfg.cols 
    pack_v_max = 4.2 * cfg.rows
    pack_v_min = 2.5 * cfg.rows
    
    pack_ah = cfg.cols * cfg.cell_cap_ah
    pack_wh = pack_ah * (cfg.rows * cfg.cell_nominal_v)
    energy_scale = pack_wh / cfg.id4_energy_wh

    # Enforce 30°C boundaries and build underlying geometry and explicit orchestrator
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols, wall_thickness=0.001)
    nema_params = build_pack_parameters(
        coolant=:water_nema2026, 
        ambient_temp=cfg.ambient_temp, 
        inlet_temp=cfg.ambient_temp, 
        flow_rate=1e-6
    )
    
    geom = build_pack_geometry(rows=cfg.rows, cols=cfg.cols, cell_pitch=0.025)
    
    cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows, cfg.cols; verbose=false, geom_sigma=0.005)
    
    # Properly instantiated dummy AnticipativeTMS to satisfy strict typing
    reactive_strat = ReactiveTMS(ambient_temp=cfg.ambient_temp, T_high=35.0, T_low=32.0)
    anticipative_strat = AnticipativeTMS(ambient_temp=cfg.ambient_temp, T_high=35.0, T_low=32.0, cell_load_threshold=5.0, lookahead_window=300.0)
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=Dict(:grid => :reactive, :drive => :reactive, :aging => :reactive, :smooth => :reactive, :rest => :reactive))
    
    # Apply stress multiplier to V2G power and parse WLTP driving cycle
    t_fcr, f_fcr = load_fcr_data(13.5 * 3600.0)
    
    max_v2g_power_w = cfg.id4_fcr_power_w * energy_scale * stress_multiplier
    
    wltp_path = joinpath("data", "V2G", "driving_power_wltp.csv")
    if isfile(wltp_path)
        f_wltp = Matrix(CSV.read(wltp_path, DataFrame))
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

    # Construct daily experimental schedule and simulate pack
    v2g_daily_schedule = Step[
        wltp_step,                                                              
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 8.5 * 3600.0, :smooth), 
        wltp_step,                                                              
        fcr_step,                  
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 1.0 * 3600.0, :smooth)  
    ]
    exp = Experiment(v2g_daily_schedule)
    
    # High fidelity time stepping to capture stress profile smoothly
    df = simulate_pack!(cosim, exp, nema_params, tms; 
        total_time=24.0 * 3600.0, save_csv=false, print_prefix="[$(Int(stress_multiplier))x FCR] ",
        force_even_current=false, is_isothermal=false, geom_sigma=0.005,
        m_active=0.03, m_passive=1e-6, heat_multiplier=1.0,
        
        opt_w_ratio=Dict{Symbol, Float64}(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
        opt_alpha_elec=Dict{Symbol, Float64}(:smooth => 7.5, :drive => 17.5, :grid => 1.9375, :aging => 5.0, :rest => 2.5), 
        opt_dt_max_elec=Dict{Symbol, Float64}(:smooth => 250.0, :drive => 0.3125, :grid => 3.125, :aging => 100.0, :rest => 900.0),
        opt_alpha_therm=Dict{Symbol, Float64}(:smooth => 0.025, :drive => 3.0, :grid => 3.0, :aging => 0.125, :rest => 0.025),
        opt_dt_max_therm=Dict{Symbol, Float64}(:smooth => 139.575, :drive => 37.5, :grid => 37.5, :aging => 515.625, :rest => 900.0),
        
        dense_logging=true, sparse_logging=false, ultra_sparse_logging=false, verbose=true
    )
    
    CSV.write(out_path, df)
    println("\n>>> Successfully generated and saved $out_path")
end

"""
    plot_stress_dashboard(cfg)

Render 3-panel dense comparison dashboard for 1x vs 10x FCR stress test.

# Arguments
- `cfg`: Global configuration dashboard parameters

# Returns
- Nothing
"""
function plot_stress_dashboard(cfg)
    # Create output directory and invoke dataset generation if missing
    base_dir = joinpath(pwd(), cfg.output_folder)
    if !isdir(base_dir) mkpath(base_dir) end

    dense_1x_path = joinpath(base_dir, "dense_stress_24h_1x.csv")
    dense_10x_path = joinpath(base_dir, "dense_stress_24h_10x.csv")
    
    if !isfile(dense_1x_path) Base.invokelatest(generate_dense_24h_stress_csv, cfg, dense_1x_path, 1.0) end
    if !isfile(dense_10x_path) Base.invokelatest(generate_dense_24h_stress_csv, cfg, dense_10x_path, 10.0) end
    
    println("\n=======================================================")
    println(">>> STRESS TEST DASHBOARD (1x vs 10x FCR at 30°C)")
    println("=======================================================")
    
    # Load simulated traces and format terminal table header
    df_1x = CSV.read(dense_1x_path, DataFrame)
    df_10x = CSV.read(dense_10x_path, DataFrame)
    
    t_hrs_1x = df_1x.Time_s ./ 3600.0
    t_hrs_10x = df_10x.Time_s ./ 3600.0
    
    header = @sprintf("%-15s | %-20s", "Scenario", "Peak temp (°C)")
    println(header)
    println(repeat("-", length(header)))
    @printf("%-15s | %14.2f\n", "1x FCR", maximum(df_1x.Max_Temp_C))
    @printf("%-15s | %14.2f\n", "10x FCR", maximum(df_10x.Max_Temp_C))
    
    # Render pack voltage plot drawing 10x trace first followed by 1x dashed overlay
    p_v = plot(title="Pack voltage", ylabel="Voltage [V]", xlabel="", legend=:bottomright)
    plot!(p_v, t_hrs_10x, df_10x.Pack_Voltage_V, label="10x FCR", color=MATLAB_COLORS[2], lw=2)
    plot!(p_v, t_hrs_1x, df_1x.Pack_Voltage_V, label="1x FCR", color=MATLAB_COLORS[1], lw=2.5, ls=:dash)
    
    # Render pack current plot drawing 10x trace first followed by 1x dashed overlay
    p_i = plot(title="Pack current", ylabel="Current [A]", xlabel="", legend=:bottomright)
    plot!(p_i, t_hrs_10x, df_10x.Pack_Current_A, label="10x FCR", color=MATLAB_COLORS[2], lw=2)
    plot!(p_i, t_hrs_1x, df_1x.Pack_Current_A, label="1x FCR", color=MATLAB_COLORS[1], lw=2.5, ls=:dash)
    
    # Render maximum cell temperature plot with trigger limit reference
    p_t = plot(title="Maximum cell temperature", ylabel="Temperature [°C]", xlabel="Time [hours]", legend=:topleft)
    plot!(p_t, t_hrs_10x, df_10x.Max_Temp_C, label="10x FCR", color=MATLAB_COLORS[2], lw=2)
    plot!(p_t, t_hrs_1x, df_1x.Max_Temp_C, label="1x FCR", color=MATLAB_COLORS[1], lw=2.5, ls=:dash)
    
    hline!(p_t, [35.0], label="Trigger (35°C)", lw=2, ls=:dot, color=:black)
    
    # Combine plots into stacked layout and export to disk
    dash_dense = plot(p_v, p_i, p_t, layout=(3,1), size=(800, 900), margin=6Plots.mm)
    out_name = joinpath(base_dir, "stress_dense_comparison_24h.svg")
    savefig(dash_dense, out_name)
    println("\n>>> Exported 24h stress test comparison traces to '$out_name'")
end

Base.invokelatest(plot_stress_dashboard, CONFIG)