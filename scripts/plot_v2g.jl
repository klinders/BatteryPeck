# ==============================================================================
# plot_v2g.jl
# Generate 6-year SoH lifecycle, voltage drift plots, and print EFC metrics.
# Contains:
# 1. load_fcr_data: Helper function to retrieve grid frequency profile
# 2. generate_dense_24h_csv: Automatically run a 1-day dense simulation if the CSV is missing
# 3. extract_v2g_metrics: Calculate final SoH, equivalent full cycles, and peak voltage drift
# 4. plot_v2g_dashboard: Render lifecycle figures and print terminal summary
# ==============================================================================

using Revise
using CSV
using DataFrames
using Plots
using Statistics
using Printf

# Import simulation libraries to generate the dense CSV if missing
using BatteryToolkit
using ModelingToolkit
using OrdinaryDiffEq

Revise.revise()

# ==============================================================================
# --- USER CONFIGURATION DASHBOARD ---
# ==============================================================================
CONFIG = (
    # Specify the folders where the respective master CSVs are located.
    folder_fcr  = "lifecycle_run_1783201650", 
    folder_rest = "lifecycle_run_1783235664",
    
    # Internal pack constants for the generator
    rows = 4,
    cols = 7,
    cell_cap_ah = 5.0,
    cell_nominal_v = 3.7,
    id4_energy_wh = 62000.0,
    id4_fcr_power_w = 11000.0
)

# Standard MATLAB Colour Palette (1=Blue, 2=Orange)
MATLAB_COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]
# ==============================================================================

"""
    load_fcr_data(duration)

Helper function to retrieve grid frequency profile.
"""
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
    generate_dense_24h_csv(cfg, out_path, is_fcr)

Automatically run a 1-day dense simulation if the CSV is missing.
"""
function generate_dense_24h_csv(cfg, out_path::String, is_fcr::Bool)
    mode_str = is_fcr ? "FCR enabled" : "REST enabled"
    println("\n[!] Dense CSV not found. Generating now for V2G ($mode_str) (High fidelity)...")
    
    elec_params = Chen2020()
    elec_params.Vmin = 2.5
    elec_params.Vmax = 4.2
    
    soc_init = 0.65 
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    pack_1C_amps = cfg.cell_cap_ah * cfg.cols 
    pack_v_max = 4.2 * cfg.rows
    pack_v_min = 2.5 * cfg.rows
    
    pack_ah = cfg.cols * cfg.cell_cap_ah
    pack_wh = pack_ah * (cfg.rows * cfg.cell_nominal_v)
    energy_scale = pack_wh / cfg.id4_energy_wh

    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=cfg.cols, wall_thickness=0.001)
    nema_params = build_pack_parameters(
        coolant=:water_nema2026, 
        ambient_temp=298.15, 
        inlet_temp=298.15, 
        flow_rate=1e-6
    )
    
    geom = build_pack_geometry(rows=cfg.rows, cols=cfg.cols, cell_pitch=0.025)
    cosim = build_pack_simulator(geom, elec_params, nema_params, cfg.rows, cfg.cols; verbose=false, geom_sigma=0.005)
    
    reactive_strat = ReactiveTMS(ambient_temp=298.15, T_high=35.0, T_low=32.0)
    anticipative_strat = AnticipativeTMS(ambient_temp=298.15, T_high=35.0, T_low=32.0, cell_load_threshold=5.0, lookahead_window=300.0)
    tms = HybridTMS(reactive=reactive_strat, anticipative=anticipative_strat, regime_map=Dict(:grid => :reactive, :drive => :reactive, :aging => :reactive, :smooth => :reactive, :rest => :reactive))
    
    t_fcr, f_fcr = load_fcr_data(13.5 * 3600.0)
    max_v2g_power_w = cfg.id4_fcr_power_w * energy_scale 
    
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

    fcr_or_rest = is_fcr ? FCRStep(t_fcr, f_fcr, max_v2g_power_w, 13.5 * 3600.0) : RestStep(13.5 * 3600.0)

    v2g_daily_schedule = Step[
        wltp_step,                                                              
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 8.5 * 3600.0, :smooth), 
        wltp_step,                                                              
        fcr_or_rest,                  
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 1.0 * 3600.0, :smooth)  
    ]
    exp = Experiment(v2g_daily_schedule)
    
    df = simulate_pack!(cosim, exp, nema_params, tms; 
        total_time=24.0 * 3600.0, save_csv=false, print_prefix="[Dense single cycle] ",
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
    extract_v2g_metrics(df, num_cells, pack_ah)

Calculate final SoH, equivalent full cycles, and peak voltage drift.
"""
function extract_v2g_metrics(df::DataFrame, num_cells::Int, pack_ah::Float64)
    soh_cols = [Symbol("SoH_Cell_$i") for i in 1:num_cells]
    final_soh = mean(Float64.(Matrix(df[end:end, soh_cols])))
    
    dt_arr = diff(df.Time_s)
    I_arr = abs.(df.Pack_Current_A[1:end-1])
    ah_throughput = sum(I_arr .* dt_arr) / 3600.0
    efc = ah_throughput / (pack_ah * 2.0)
    
    V_mat = Float64.(Matrix(df[!, [Symbol("V_Cell_$i") for i in 1:num_cells]]))
    max_drift_mv = maximum(maximum(V_mat, dims=2) .- minimum(V_mat, dims=2)) * 1000.0
    
    return final_soh, efc, max_drift_mv
end

"""
    extract_daily_voltage_drift(df, num_cells)

Extract the maximum voltage spread (in mV) recorded per 24-hour block.
"""
function extract_daily_voltage_drift(df::DataFrame, num_cells::Int)
    V_cols = [Symbol("V_Cell_$i") for i in 1:num_cells]
    
    # Calculate continuous spread in mV
    V_mat = Float64.(Matrix(df[!, V_cols]))
    spread_mv = (maximum(V_mat, dims=2) .- minimum(V_mat, dims=2)) .* 1000.0
    
    t_days_continuous = (df.Time_s .- df.Time_s[1]) ./ (24 * 3600.0)
    max_day = ceil(Int, t_days_continuous[end])
    
    daily_drift = zeros(max_day)
    day_indices = collect(1:max_day)
    
    for day in 1:max_day
        # Filter for indices within the current 24-hour window
        idx = findall(x -> (x >= day - 1) && (x < day), t_days_continuous)
        if !isempty(idx)
            daily_drift[day] = maximum(spread_mv[idx])
        else
            daily_drift[day] = day > 1 ? daily_drift[day-1] : 0.0
        end
    end
    
    return day_indices, daily_drift
end

"""
    plot_v2g_dashboard(cfg)

Render lifecycle figures and print terminal summary.
"""
function plot_v2g_dashboard(cfg)
    dir_fcr = joinpath(pwd(), "results", cfg.folder_fcr)
    dir_rest = joinpath(pwd(), "results", cfg.folder_rest)
    
    if !isdir(dir_fcr) || !isdir(dir_rest)
        error("[!] Cannot find one or both folders:\n$dir_fcr\n$dir_rest\nPlease check your folder names in CONFIG.")
    end

    dense_fcr_path = joinpath(dir_fcr, "dense_v2g_24h_FCR.csv")
    dense_rest_path = joinpath(dir_rest, "dense_v2g_24h_REST.csv")
    
    if !isfile(dense_fcr_path) Base.invokelatest(generate_dense_24h_csv, cfg, dense_fcr_path, true) end
    if !isfile(dense_rest_path) Base.invokelatest(generate_dense_24h_csv, cfg, dense_rest_path, false) end
    
    # 5-year trunk from the REST folder, V2G phases from their respective folders
    path_5yr = joinpath(dir_rest, "aging_5_years_master.csv")
    path_1yr_rest = joinpath(dir_rest, "v2g_master_data.csv")
    path_1yr_fcr = joinpath(dir_fcr, "v2g_master_data.csv")

    println("\n=======================================================")
    println(">>> 6-YEAR LIFECYCLE V2G SUMMARY")
    println("=======================================================")
    
    header = @sprintf("%-15s | %-15s | %-12s | %-15s", "Scenario", "End SoH (%)", "1-Yr EFC", "Peak drift (mV)")
    println(header)
    println(repeat("-", length(header)))
    
    if !isfile(path_5yr)
        println("[!] Could not find 5-year aging CSV at: $path_5yr")
        return
    end
    
    df_5yr = CSV.read(path_5yr, DataFrame)
    soh_cols = [Symbol("SoH_Cell_$i") for i in 1:cfg.rows*cfg.cols]
    avg_soh_5yr = mean(Float64.(Matrix(df_5yr[!, soh_cols])), dims=2)[:, 1]
    time_5yr_yrs = df_5yr.Time_s ./ (365 * 24 * 3600.0)
    
    # Base Aging Trunk (Black)
    p_soh = plot(time_5yr_yrs, avg_soh_5yr, label="Base aging (0-5 yrs)", color=:black, lw=2, xlabel="Time [years]", ylabel="Pack average SoH [%]", title="6-year lifecycle degradation", legend=:bottomleft)
    
    p_drift = plot(title="Daily maximum voltage drift (FCR phase)", xlabel="Time [days]", ylabel="Max cell spread [mV]", legend=:topleft)

    # Branch 1: REST (Blue)
    if isfile(path_1yr_rest)
        df_rest = CSV.read(path_1yr_rest, DataFrame)
        final_soh, efc, drift = Base.invokelatest(extract_v2g_metrics, df_rest, cfg.rows*cfg.cols, cfg.cols*cfg.cell_cap_ah)
        println(@sprintf("%-15s | %14.2f%% | %12.1f | %14.2f", "Rest", final_soh, efc, drift))
        
        avg_soh_rest = mean(Float64.(Matrix(df_rest[!, soh_cols])), dims=2)[:, 1]
        
        # Removed the t_offset addition. The checkpointed time natively continues from 5 years.
        plot!(p_soh, df_rest.Time_s ./ (365 * 24 * 3600.0), avg_soh_rest, label="Rest", color=MATLAB_COLORS[1], lw=3, ls=:dash)
        
        # Extract and plot daily max drift
        t_days_rest, drift_rest = extract_daily_voltage_drift(df_rest, cfg.rows*cfg.cols)
        plot!(p_drift, t_days_rest, drift_rest, label="Rest", color=MATLAB_COLORS[1], lw=2)
    end
    
    # Branch 2: FCR (Orange)
    if isfile(path_1yr_fcr)
        df_fcr = CSV.read(path_1yr_fcr, DataFrame)
        final_soh, efc, drift = Base.invokelatest(extract_v2g_metrics, df_fcr, cfg.rows*cfg.cols, cfg.cols*cfg.cell_cap_ah)
        println(@sprintf("%-15s | %14.2f%% | %12.1f | %14.2f", "FCR", final_soh, efc, drift))
        
        avg_soh_fcr = mean(Float64.(Matrix(df_fcr[!, soh_cols])), dims=2)[:, 1]
        
        # Removed the t_offset addition.
        plot!(p_soh, df_fcr.Time_s ./ (365 * 24 * 3600.0), avg_soh_fcr, label="FCR", color=MATLAB_COLORS[2], lw=2)
        
        # Extract and plot daily max drift
        t_days_fcr, drift_fcr = extract_daily_voltage_drift(df_fcr, cfg.rows*cfg.cols)
        plot!(p_drift, t_days_fcr, drift_fcr, label="FCR", color=MATLAB_COLORS[2], lw=2)
        
        savefig(p_drift, "v2g_voltage_drift.svg")
        println("\n>>> Exported Daily Voltage Drift plot to 'v2g_voltage_drift.svg'")
    end
    
    savefig(p_soh, "v2g_6year_lifecycle.svg")
    println(">>> Exported 6-Year Lifecycle plot to 'v2g_6year_lifecycle.svg'")
    
    # ----------------------------------------------------------------------
    # Render Dense Dashboard (FCR vs REST Traces overlay)
    # ----------------------------------------------------------------------
    if isfile(dense_fcr_path) && isfile(dense_rest_path)
        df_dense_fcr = CSV.read(dense_fcr_path, DataFrame)
        df_dense_rest = CSV.read(dense_rest_path, DataFrame)
        
        t_hrs_fcr = df_dense_fcr.Time_s ./ 3600.0
        t_hrs_rest = df_dense_rest.Time_s ./ 3600.0
        
        # 1. Voltage Plot (FCR solid orange first, then REST dashed blue on top)
        p_v = plot(title="Pack voltage", ylabel="Voltage [V]", xlabel="", legend=:bottomright)
        plot!(p_v, t_hrs_fcr, df_dense_fcr.Pack_Voltage_V, label="FCR", color=MATLAB_COLORS[2], lw=2)
        plot!(p_v, t_hrs_rest, df_dense_rest.Pack_Voltage_V, label="Rest", color=MATLAB_COLORS[1], lw=2.5, ls=:dash)
        
        # 2. Current Plot
        p_i = plot(title="Pack current", ylabel="Current [A]", xlabel="", legend=:bottomright)
        plot!(p_i, t_hrs_fcr, df_dense_fcr.Pack_Current_A, label="FCR", color=MATLAB_COLORS[2], lw=2)
        plot!(p_i, t_hrs_rest, df_dense_rest.Pack_Current_A, label="Rest", color=MATLAB_COLORS[1], lw=2.5, ls=:dash)
        
        # 3. Temperature Plot
        p_t = plot(title="Maximum cell temperature", ylabel="Temperature [°C]", xlabel="Time [hours]", legend=:topright)
        plot!(p_t, t_hrs_fcr, df_dense_fcr.Max_Temp_C, label="FCR", color=MATLAB_COLORS[2], lw=2)
        plot!(p_t, t_hrs_rest, df_dense_rest.Max_Temp_C, label="Rest", color=MATLAB_COLORS[1], lw=2.5, ls=:dash)
        
        dash_dense = plot(p_v, p_i, p_t, layout=(3,1), size=(800, 900), margin=6Plots.mm)
        out_name = "v2g_dense_comparison_24h.svg"
        savefig(dash_dense, out_name)
        println(">>> Exported 24h dense V2G comparison traces to '$out_name'")
    end
end

Base.invokelatest(plot_v2g_dashboard, CONFIG)