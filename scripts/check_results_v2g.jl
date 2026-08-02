# ==============================================================================
# check_results_v2g.jl
# Dedicated dashboard for analyzing 1-Year V2G degradation and thermal performance.
# Automatically isolates the V2G phase from multi-year CSV files.
# Contains:
# 1. smart_decimate: Drop redundant steady state points to compress array size
# 2. analyze_v2g_scenario: Isolate V2G phase and calculate degradation metrics
# ==============================================================================

using CSV
using DataFrames
using Plots
using Statistics
using Printf

# Configure dashboard setup variables
CONFIG = (
    # Assign path to master results dataset
    target_path = "results/lifecycle_run_1782613422/v2g_1_year_master.csv", 
    
    # Assign thermal management thresholds for plot visualisation
    tms_t_high = 35.0,
    tms_t_low = 32.0,
    
    export_svg = true
)

"""
    smart_decimate(df::DataFrame, dt_thresh, dV_thresh, dI_thresh, dT_thresh)

Drop redundant steady state points to compress array size.

# Arguments
- `df::DataFrame`: Base dataset containing full time sequence
- `dt_thresh`: Threshold bound for elapsed time tracking
- `dV_thresh`: Threshold bound for tracked voltage dynamics
- `dI_thresh`: Threshold bound for tracked current dynamics
- `dT_thresh`: Threshold bound for tracked temperature dynamics

# Returns
- Decimated output data frame
"""
function smart_decimate(df::DataFrame, dt_thresh, dV_thresh, dI_thresh, dT_thresh)
    if nrow(df) <= 2 return df end

    keep_idx = Int[1]
    sizehint!(keep_idx, div(nrow(df), 10))

    # Extract data frame columns to native vectors to ensure type stability and accelerate execution
    t_vec = df.Time_s
    v_vec = df.Pack_Voltage_V
    i_vec = df.Pack_Current_A
    temp_vec = df.Max_Temp_C

    last_t = t_vec[1]
    last_v = v_vec[1]
    last_i = i_vec[1]
    last_temp = temp_vec[1]

    # Evaluate absolute thresholds against sequential vector data arrays
    for i in 2:(nrow(df)-1)
        t = t_vec[i]
        v = v_vec[i]
        cur = i_vec[i]
        temp = temp_vec[i]

        if (t - last_t) >= dt_thresh || abs(v - last_v) >= dV_thresh || abs(cur - last_i) >= dI_thresh || abs(temp - last_temp) >= dT_thresh
            push!(keep_idx, i)
            last_t = t; last_v = v; last_i = cur; last_temp = temp
        end
    end
    
    push!(keep_idx, nrow(df)) 
    return df[keep_idx, :]
end

"""
    analyze_v2g_scenario(cfg)

Isolate V2G phase and calculate degradation metrics.

# Arguments
- `cfg`: Mapped tuple containing target paths and bounds

# Returns
- Nothing
"""
function analyze_v2g_scenario(cfg)
    # Announce initialisation block and pull external CSV payload
    println("\n=======================================================")
    println(">>> INITIATING V2G DATA ANALYSIS")
    println("=======================================================")
    
    println("[1/6] Loading CSV Data from disk... (This may take a moment)")
    t0 = time()
    df_raw = CSV.read(cfg.target_path, DataFrame)
    
    if nrow(df_raw) == 0
        error("The loaded CSV file is entirely empty.")
    end
    @printf("      -> Loaded %d rows in %.2f seconds\n", nrow(df_raw), time() - t0)

    # Isolate macroscopic data phase
    println("[2/6] Isolating V2G Phase and filtering 24-Hour window...")
    t1 = time()
    max_time = maximum(df_raw.Time_s)
    v2g_start_time = max(0.0, max_time - (365 * 24 * 3600.0))
    
    df_v2g = filter(row -> row.Time_s >= v2g_start_time, df_raw)
    
    # Normalise time array to begin at zero
    df_v2g.Time_s .-= df_v2g.Time_s[1] 
    
    # Isolate microscopic data phase
    v2g_max_time = maximum(df_v2g.Time_s)
    cutoff_24h = v2g_max_time - (24 * 3600.0)
    df_24h = filter(row -> row.Time_s >= cutoff_24h, df_v2g)
    @printf("      -> Processed in %.2f seconds\n", time() - t1)
    
    # Generate decimated sets for macro and micro inspection tracking
    println("[3/6] Smart Decimating 1-Year Macroscopic Data...")
    t2 = time()
    df_v2g_dec = smart_decimate(df_v2g, 3600.0, 0.5, 2.0, 0.5)
    @printf("      -> Compressed %d rows to %d rows in %.2f seconds\n", nrow(df_v2g), nrow(df_v2g_dec), time() - t2)

    println("[4/6] Smart Decimating 24-Hour Microscopic Data...")
    t3 = time()
    df_24h_dec = smart_decimate(df_24h, 60.0, 0.01, 0.1, 0.05)
    @printf("      -> Compressed %d rows to %d rows in %.2f seconds\n", nrow(df_24h), nrow(df_24h_dec), time() - t3)

    # Calculate macroscopic metrics over full duration
    println("[5/6] Calculating Analytical Metrics and Spreads...")
    t4 = time()
    soh_cols = filter(name -> startswith(name, "SoH_Cell_"), names(df_v2g_dec))
    soh_matrix = Matrix(df_v2g_dec[:, soh_cols])
    
    avg_soh = mean(soh_matrix, dims=2)[:]
    spread_soh = (maximum(soh_matrix, dims=2) .- minimum(soh_matrix, dims=2))[:]
    
    p1 = plot(df_v2g_dec.Time_s ./ 86400.0, avg_soh, 
        title="Average Pack SoH (V2G Year)", ylabel="SoH [%]", xlabel="Time [days]", 
        color=:purple, lw=2, legend=false)
        
    p2 = plot(df_v2g_dec.Time_s ./ 86400.0, spread_soh, 
        title="Maximum SoH Spread (V2G Year)", ylabel="Δ SoH [%]", xlabel="Time [days]", 
        color=:orange, lw=2, legend=false)
        
    # Calculate daily average temperature
    df_v2g_dec.Day = floor.(Int, df_v2g_dec.Time_s ./ 86400.0)
    daily_temp_df = combine(groupby(df_v2g_dec, :Day), :Max_Temp_C => mean => :Avg_Temp_C)
    
    p3 = plot(daily_temp_df.Day, daily_temp_df.Avg_Temp_C, 
        title="Daily Average Temperature", ylabel="Temp [°C]", xlabel="Time [days]", 
        color=:firebrick, lw=2, linetype=:steppost, legend=false)

    # Calculate microscopic metrics for final period
    p4 = plot(df_24h_dec.Time_s ./ 3600.0, df_24h_dec.Pack_Voltage_V, 
        label="Voltage", color=:blue, lw=2,
        ylabel="Voltage [V]", title="Electrical Traces (Last 24h)", xlabel="Time [h]", legend=:topleft)
    
    p4_twin = twinx()
    plot!(p4_twin, df_24h_dec.Time_s ./ 3600.0, df_24h_dec.Pack_Current_A, 
        label="Current", color=:red, lw=1.5, linestyle=:solid, 
        ylabel="Current [A]", legend=:bottomright)
        
    p5 = plot(df_24h_dec.Time_s ./ 3600.0, df_24h_dec.Max_Temp_C, 
        label="Max Core Temp", color=:firebrick, lw=2,
        ylabel="Temperature [°C]", title="Thermal Traces (Last 24h)", xlabel="Time [h]")
        
    hline!(p5, [cfg.tms_t_high], color=:black, linestyle=:dash, label="TMS High ($(cfg.tms_t_high)°C)")
    hline!(p5, [cfg.tms_t_low], color=:gray, linestyle=:dash, label="TMS Low ($(cfg.tms_t_low)°C)")
    
    # Calculate cell voltage spread as state of charge proxy
    v_cols = filter(name -> startswith(name, "V_Cell_"), names(df_24h_dec))
    v_matrix = Matrix(df_24h_dec[:, v_cols])
    v_spread_mv = (maximum(v_matrix, dims=2) .- minimum(v_matrix, dims=2)) .* 1000.0
    
    p6 = plot(df_24h_dec.Time_s ./ 3600.0, v_spread_mv, 
        title="Cell Voltage Spread (Last 24h)", ylabel="Spread [mV]", xlabel="Time [h]", 
        color=:teal, lw=2, legend=false)
    @printf("      -> Calculated in %.2f seconds\n", time() - t4)

    # Render dashboard and print terminal metrics
    println("[6/6] Rendering Dashboard and Exporting...")
    t5 = time()
    dashboard = plot(p1, p2, p3, p4, p5, p6, layout=(2, 3), size=(1800, 800), margin=5Plots.mm)
    display(dashboard)
    
    if cfg.export_svg
        out_name = "v2g_analysis_dashboard.svg"
        savefig(dashboard, out_name)
        @printf("      -> Exported to %s\n", out_name)
    end
    @printf("      -> Rendering complete in %.2f seconds\n", time() - t5)
    
    # Print final simulation metrics
    final_pack_soc = round(df_24h.Pack_SoC[end], digits=2)
    final_v_spread = round((maximum(Array(df_24h[end, v_cols])) - minimum(Array(df_24h[end, v_cols]))) * 1000.0, digits=2)
    
    println("\n=======================================================")
    println(">>> FINAL V2G METRICS (DAY 365)")
    println("=======================================================")
    println(" Final Pack SoC        : $final_pack_soc %")
    println(" Final Voltage Spread  : $final_v_spread mV")
    println(" Final Average SoH     : $(round(avg_soh[end], digits=2)) %")
    println(" Final SoH Spread      : $(round(spread_soh[end], digits=3)) %")
    println("=======================================================\n")
end

Base.invokelatest(analyze_v2g_scenario, CONFIG)