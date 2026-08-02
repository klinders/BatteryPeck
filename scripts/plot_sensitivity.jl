# ==============================================================================
# plot_sensitivity.jl
# Render Tornado Dashboard for the Sensitivity Analysis
# ==============================================================================

using CSV
using DataFrames
using Plots
using Printf

CONFIG = (
    csv_path = "results/sensitivity_run_1783434663/sensitivity_results.csv",
    x_limits = (-10.0, 10.0) # Locks all 3 plots to the same scale for visual comparison
)

# Standard MATLAB Colour Palette
MATLAB_COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]

function plot_sensitivity_dashboard(cfg)
    println("\n=======================================================")
    println(">>> SENSITIVITY ANALYSIS (TORNADO DASHBOARD)")
    println("=======================================================")
    
    if !isfile(cfg.csv_path)
        error("[!] Could not find CSV at: $(cfg.csv_path)")
    end
    
    df = CSV.read(cfg.csv_path, DataFrame)
    param_names = df.Parameter
    y_pos = 1:length(param_names)
    y_ticks = (y_pos, param_names)
    
    # Print formatted terminal table
    header = @sprintf("%-20s | %11s | %11s | %11s | %11s | %11s | %11s", "Parameter", "Peak T Low", "Peak T High", "Grad Low", "Grad High", "Dose Low", "Dose High")
    println(header)
    println(repeat("-", length(header)))
    
    for i in 1:nrow(df)
        line = @sprintf("%-20s | %10.2f%% | %10.2f%% | %10.2f%% | %10.2f%% | %10.2f%% | %10.2f%%", 
            df.Parameter[i], 
            df.Peak_T_Low_Pct[i], df.Peak_T_High_Pct[i],
            df.Max_Grad_Low_Pct[i], df.Max_Grad_High_Pct[i],
            df.Dose_Low_Pct[i], df.Dose_High_Pct[i]
        )
        println(line)
    end
    
    # Tightly crop the y-axis to eliminate the blank space at the bottom
    y_lim_tight = (0.5, length(param_names) + 0.5)
    
    # Render synced Tornado Dashboard
    p1 = bar(y_pos, df.Peak_T_Low_Pct, orientation=:h, yticks=y_ticks, label="-10% variance", color=MATLAB_COLORS[1], alpha=0.85, title="Peak temperature", xlabel="% change", legend=:bottomright, xlims=cfg.x_limits, ylims=y_lim_tight)
    bar!(p1, y_pos, df.Peak_T_High_Pct, orientation=:h, label="+10% variance", color=MATLAB_COLORS[2], alpha=0.85)
    
    p2 = bar(y_pos, df.Max_Grad_Low_Pct, orientation=:h, yticks=y_ticks, label="", color=MATLAB_COLORS[1], alpha=0.85, title="Max temperature gradient", xlabel="% change", xlims=cfg.x_limits, ylims=y_lim_tight)
    bar!(p2, y_pos, df.Max_Grad_High_Pct, orientation=:h, label="", color=MATLAB_COLORS[2], alpha=0.85)
    
    p3 = bar(y_pos, df.Dose_Low_Pct, orientation=:h, yticks=y_ticks, label="", color=MATLAB_COLORS[1], alpha=0.85, title="Thermal dose", xlabel="% change", xlims=cfg.x_limits, ylims=y_lim_tight)
    bar!(p3, y_pos, df.Dose_High_Pct, orientation=:h, label="", color=MATLAB_COLORS[2], alpha=0.85)
    
    dashboard = plot(p1, p2, p3, layout=(3,1), size=(1000, 1200), margin=8Plots.mm)
    display(dashboard)
    savefig(dashboard, "sensitivity_tornado_dashboard.svg")
    println("\n>>> Dashboard exported to 'sensitivity_tornado_dashboard.svg'")
end

plot_sensitivity_dashboard(CONFIG)