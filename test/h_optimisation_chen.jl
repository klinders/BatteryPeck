# ==============================================================================
# h_optimisation_chen.jl
# Performs 1D grid search to find optimal convective heat transfer coefficient (h)
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using CSV
using DataFrames
using Plots
using Plots.Measures
using DataInterpolations
using Statistics

using Revise
using BatteryPeck

Revise.revise()

function optimise_h_conv()
    # Load experimental data
    data_file = joinpath(@__DIR__, "..", "data", "Chen2020", "LGM50_cell03.csv")
    df = CSV.read(data_file, DataFrame, skipto=15, header=14)

    t_exp = df[!, "Test Time [s]"]
    T_surf_exp = df[!, "Temperature Cell [degC]"]
    
    T_chamber_exp = df[!, "Temperature Chamber [degC]"]
    T_amb_start_K = T_chamber_exp[1] + 273.15

    exp_duration = t_exp[end]

    # Setup parameters and initial state
    p = Chen2020()
    soc_init = 0.063
    p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
    p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
    p.Vmin = 2.4
    p.Vmax = 4.3

    # Define single experiment
    exp = Experiment([CurrentDriveStep(data_file, exp_duration)])

    # Sweep h from 10 to 30 W/m²K in steps of 0.5
    h_range = 10.0:0.5:30.0
    rmse_results = Float64[]

    println("Starting 1D grid search for optimal h_conv...")

    for h_test in h_range
        println("Evaluating h_conv = $h_test W/m²K ...")
        
        # Build system with current h_test
        @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=h_test, T_ambient=T_amb_start_K)
        
        # Run simulation
        sol = simulate(sys, exp, saveat=10.0)
        
        # Extract results
        t_sim = sol.t
        T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15
        
        # Interpolate simulated temperature curve
        sim_temp_interp = LinearInterpolation(T_shell_sim, t_sim)
        
        # Map simulation to exact experimental timestamps
        T_sim_mapped = sim_temp_interp.(clamp.(t_exp, t_sim[1], t_sim[end]))
        
        # Calculate Root Mean Squared Error (RMSE) against experimental data
        rmse = sqrt(mean((T_sim_mapped .- T_surf_exp).^2))
        
        println("  -> RMSE = $(round(rmse, digits=3)) °C")
        push!(rmse_results, rmse)
    end
    
    # Return raw data to save in global workspace
    return collect(h_range), rmse_results
end

# Run computation once to obtain data
h_range, rmse_results = optimise_h_conv()

# Plot results
best_idx = argmin(rmse_results)
best_h = h_range[best_idx]
best_rmse = rmse_results[best_idx]

println("\n Grid search complete.")
println("Optimal h_conv: $best_h W/m²K")
println("Lowest RMSE:    $(round(best_rmse, digits=3)) °C")

# Plot loss landscape
p_loss = plot(h_range, rmse_results, marker=:circle, lw=2, color="#D95319",
              xlabel="Convective heat transfer coefficient, h (W/m²K)", 
              ylabel="Temperature RMSE (°C)",
              title="1D grid search: thermal validation error",
              legend=false, margin=5mm)
              
# Add point to highlight minimum
scatter!(p_loss, [best_h], [best_rmse], color=:red, markersize=6, 
         annotations=(best_h, best_rmse + 5, text("Optimal: $best_h", :center, :bottom, 10)))
         
display(p_loss)