# ==============================================================================
# h_optimisation_chen.jl
# Perform grid search to find optimal convective heat transfer coefficient
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using CSV
using DataFrames
using Plots
using Plots.Measures
using DataInterpolations
using Statistics
using Printf

using Revise
using BatteryToolkit

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

    # Set up parameters and initial state
    p = Chen2020()
    soc_init = 0.063
    p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
    p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
    p.Vmin = 2.4
    p.Vmax = 4.3

    println("Building and compiling mathematical system...")
    # Compile mathematical system outside loop
    @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=15.0, T_ambient=T_amb_start_K)
    
    # Define single experiment
    exp = Experiment([CurrentDriveStep(data_file, exp_duration)])
    
    println("Compiling base ordinary differential equation problem...")
    # Compile base ordinary differential equation problem outside loop
    prob_base = ODEProblem(sys, [sys.Pin => 0.0, sys.Iin => 0.0], (0.0, exp.tend))

    # Define sweep range
    h_range = 17.0:0.5:21.0
    rmse_results = Float64[]
    A_cell = 0.0053
    total_evals = length(h_range)

    println("\nStarting grid search for optimal convective heat transfer coefficient...")
    t_start_total = time()

    for (idx, h_test) in enumerate(h_range)
        # Calculate new thermal resistance
        R_new = 1.0 / (h_test * A_cell)
        
        # Remake problem with new parameter instead of recompiling
        prob_test = remake(prob_base, p=[sys.R_conv.R => R_new])
        
        t_start_solve = time()
        
        # Bypass solvers to allow using remade problem directly
        integrator = init(prob_test, QNDF(); tstops=exp.tstops, saveat=10.0, save_everystep=false, reltol=1e-4, abstol=1e-7)
        for step in exp.steps
            # Force Julia to use custom package function instead of standard SciML step function
            BatteryToolkit.step!(integrator, sys, step)
            if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated
                break
            end
        end
        SciMLBase.savevalues!(integrator, true)
        sol = integrator.sol
        
        solve_time = time() - t_start_solve
        
        # Extract results
        t_sim = sol.t
        T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15
        
        # Interpolate simulated temperature curve
        sim_temp_interp = LinearInterpolation(T_shell_sim, t_sim)
        
        # Map simulation to exact experimental timestamps
        T_sim_mapped = sim_temp_interp.(clamp.(t_exp, t_sim[1], t_sim[end]))
        
        # Calculate root mean square error against experimental data
        rmse = sqrt(mean((T_sim_mapped .- T_surf_exp).^2))
        
        # Print formatted output with counter and solve time
        @printf("  [%d/%d] h_conv = %.1f W/m²K | RMSE = %.3f °C | Solve time: %.3fs\n", 
                idx, total_evals, h_test, rmse, solve_time)
        
        push!(rmse_results, rmse)
    end
    
    total_time = time() - t_start_total
    println("\nSweep completed in $(round(total_time, digits=2)) seconds.")
    
    # Plot results internally to prevent global scope binding errors
    best_idx = argmin(rmse_results)
    best_h = h_range[best_idx]
    best_rmse = rmse_results[best_idx]

    println("\nGrid search complete.")
    println("Optimal convective heat transfer coefficient: $best_h W/m²K")
    println("Lowest RMSE:    $(round(best_rmse, digits=3)) °C")

    # Plot loss landscape
    p_loss = plot(h_range, rmse_results, marker=:circle, lw=2, color="#D95319",
                  xlabel="Convective heat transfer coefficient, h (W/m²K)", 
                  ylabel="Temperature RMSE (°C)",
                  title="Grid search: thermal validation error",
                  legend=false, margin=5mm)
                  
    # Add point to highlight minimum
    scatter!(p_loss, [best_h], [best_rmse], color=:red, markersize=6, 
             annotations=(best_h, best_rmse + 0.5, text("Opt: $best_h", :center, :bottom, 10)))
             
    display(p_loss)
end

# Run computation safely bypassing world age warnings
Base.invokelatest(optimise_h_conv)