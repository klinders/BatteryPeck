# ==============================================================================
# check_results.jl
# Plotting script for ablation experiments from run_experiments.jl
# Contains:
# 1. smart_decimate: Drop redundant steady state points to collapse SVG file sizes
# 2. inspect_scenario: Generate short cycle sanity check dashboard for ablation scenarios
# ==============================================================================

using CSV
using DataFrames
using Plots

# Define folder to inspect and plot window controls
TARGET_FOLDER = "results_run_a_s6" 
MAX_DAYS_TO_PLOT = 90.0   
ZOOM_HOURS = 12.0         
EXPORT_SVG = true         

"""
    smart_decimate(df::DataFrame, dt_thresh, dV_thresh, dI_thresh, dT_thresh)

Drop redundant steady state points to collapse SVG file sizes.

# Arguments
- `df::DataFrame`: Input data frame
- `dt_thresh`: Time threshold for retention
- `dV_thresh`: Voltage change threshold
- `dI_thresh`: Current change threshold
- `dT_thresh`: Temperature change threshold

# Returns
- Decimated data frame
"""
function smart_decimate(df::DataFrame, dt_thresh, dV_thresh, dI_thresh, dT_thresh)
    # Preallocate index array to store retained rows
    keep_idx = Int[1]
    sizehint!(keep_idx, div(nrow(df), 10))

    last_t = df.Time_s[1]
    last_v = df.Pack_Voltage_V[1]
    last_i = df.Pack_Current_A[1]
    last_temp = df.Max_Temp_C[1]

    # Iterate through data frame and evaluate physical transient thresholds
    for i in 2:(nrow(df)-1)
        t = df.Time_s[i]
        v = df.Pack_Voltage_V[i]
        cur = df.Pack_Current_A[i]
        temp = df.Max_Temp_C[i]

        if (t - last_t) >= dt_thresh || 
           abs(v - last_v) >= dV_thresh || 
           abs(cur - last_i) >= dI_thresh || 
           abs(temp - last_temp) >= dT_thresh
           
            push!(keep_idx, i)
            last_t = t
            last_v = v
            last_i = cur
            last_temp = temp
        end
    end
    
    # Explicitly keep final point and return filtered data frame
    push!(keep_idx, nrow(df)) 
    return df[keep_idx, :]
end

"""
    inspect_scenario(folder_name::String)

Generate short cycle sanity check dashboard for ablation scenarios.

# Arguments
- `folder_name::String`: Target folder containing master results CSV

# Returns
- Nothing
"""
function inspect_scenario(folder_name::String)
    # Load raw data from target directory
    folder_path = joinpath(pwd(), folder_name)
    master_file = joinpath(folder_path, "master_results.csv")
    
    if !isfile(master_file)
        error("Could not find master_results.csv in $folder_path")
    end
    
    println("\n[!] Loading data from: $folder_name...")
    df_raw = CSV.read(master_file, DataFrame)
    
    # Truncate dataset to maximum allowed plot duration
    max_allowed_time = MAX_DAYS_TO_PLOT * 24 * 3600.0
    df_raw = filter(row -> row.Time_s <= max_allowed_time, df_raw)
    
    # Slice target window from end of dataset to generate zoomed view
    max_time = maximum(df_raw.Time_s)
    cutoff_time = max_time - (ZOOM_HOURS * 3600.0)
    df_zoomed_raw = filter(row -> row.Time_s >= cutoff_time, df_raw)

    # Apply dual track decimation for physical parameters and state of health to optimise SVG export
    println("[!] Decimating data points to optimize SVG export...")
    
    df_vit_full = smart_decimate(df_raw, 1800.0, 0.05, 0.5, 0.1)
    df_vit_zoomed = smart_decimate(df_zoomed_raw, 60.0, 0.01, 0.1, 0.02)
    
    df_soh_full = smart_decimate(df_raw, 12.0 * 3600.0, Inf, Inf, Inf) 
    df_soh_zoomed = smart_decimate(df_zoomed_raw, 300.0, Inf, Inf, Inf) 
    
    # Generate electrical plot with twin axes for voltage and current
    p1 = plot(df_vit_zoomed.Time_s ./ 3600.0, df_vit_zoomed.Pack_Voltage_V, 
        label="Pack voltage", color=:blue, lw=2,
        ylabel="Voltage [V]", title="Electrical traces (last $ZOOM_HOURS hours)", xlabel="Time [h]")
    
    p1_twin = twinx()
    plot!(p1_twin, df_vit_zoomed.Time_s ./ 3600.0, df_vit_zoomed.Pack_Current_A, 
        label="Pack current", color=:red, lw=1.5, linestyle=:dash, 
        ylabel="Current [A]", legend=:bottomright)
    
    # Generate thermal plot displaying maximum core temperature and coolant velocity
    p2 = plot(df_vit_full.Time_s ./ (24*3600), df_vit_full.Max_Temp_C, 
        label="Max core temp", color=:firebrick, lw=2,
        ylabel="Temperature [°C]", title="Thermal history ($MAX_DAYS_TO_PLOT days)", xlabel="Time [days]",
        legend=:topleft)
        
    hline!(p2, [45.0], color=:black, linestyle=:dot, label="High limit")
    hline!(p2, [42.0], color=:gray, linestyle=:dot, label="Low limit")
    
    p2_twin = twinx()
    plot!(p2_twin, df_vit_full.Time_s ./ (24*3600), df_vit_full.Velocity_ms, 
        label="Coolant velocity", color=:cyan, lw=1.5, linetype=:steppost,
        ylabel="Velocity [m/s]", legend=:topright, ylims=(0.0, 0.4))
    
    # Generate degradation plot showing state of health history for all cells
    p3 = plot(title="SoH history ($MAX_DAYS_TO_PLOT days)", ylabel="SoH [%]", xlabel="Time [days]", legend=false)
    
    soh_cols = filter(name -> startswith(name, "SoH_Cell_"), names(df_raw))
    colors = [:purple, :orange, :green, :magenta, :teal, :navy, :crimson]
    
    for (i, col) in enumerate(soh_cols)
        c = colors[mod1(i, length(colors))]
        plot!(p3, df_soh_full.Time_s ./ (24*3600), df_soh_full[!, col], lw=2, color=c)
    end
    
    # Generate zoomed degradation plot with dynamic padding
    min_soh = minimum([minimum(df_soh_zoomed[!, col]) for col in soh_cols])
    max_soh = maximum([maximum(df_soh_zoomed[!, col]) for col in soh_cols])
    pad = (max_soh - min_soh) == 0 ? 0.005 : (max_soh - min_soh) * 0.1
    
    p4 = plot(title="Zoomed SoH (last $ZOOM_HOURS hours)", ylabel="SoH [%]", xlabel="Time [h]", legend=false, ylims=(min_soh - pad, max_soh + pad))
    for (i, col) in enumerate(soh_cols)
        c = colors[mod1(i, length(colors))]
        plot!(p4, df_soh_zoomed.Time_s ./ 3600.0, df_soh_zoomed[!, col], lw=2, color=c)
    end
    
    # Combine plots into layout and export to SVG if configured
    dashboard = plot(p2, p1, p3, p4, layout=(2, 2), size=(1400, 800), left_margin=15Plots.mm, right_margin=15Plots.mm, bottom_margin=10Plots.mm)
    display(dashboard)
    
    if EXPORT_SVG
        out_name = folder_name * ".svg"
        savefig(dashboard, out_name)
        println("[!] Dashboard exported to: $out_name")
    else
        println("[!] Dashboard generated successfully.")
    end
end

inspect_scenario(TARGET_FOLDER)