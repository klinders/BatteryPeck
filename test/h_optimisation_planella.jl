# ==============================================================================
# h_optimisation_planella.jl
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

# Generate data for all temperatures
function optimise_h_conv_planella()
    temps = [0, 10, 25]
    files = ["Cell786_0p5C_0degC.csv", "Cell786_0p5C_10degC.csv", "Cell786_0p5C_25degC.csv"]
    
    # Custom search ranges for each temperature
    h_ranges = [
        12:0.5:14,
        14:0.5:17,
        27:0.5:31
    ]

    # Matrix for Planella pre-fitted parameters
    # Columns array temperature setting, negative electrode diffusion coefficient, ambient temperature
    fitted_params = [
         0.0  0.22e-14   0.02;
        10.0  0.40e-14   9.80;
        25.0  0.90e-14  24.45
    ]
    
    # Read voltage to state of charge look up table for initialisation
    lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
    lut = CSV.read(lut_file, DataFrame)
    v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)

    # Store results for each temperature as tuple containing range and results
    all_results = Dict{Int, Tuple{Any, Vector{Float64}}}()

    println("Starting grid search for Planella datasets...")
    t_start_total = time()

    for (T_celsius, file_name, h_range) in zip(temps, files, h_ranges)
        println("\nProcessing $(T_celsius)°C dataset ($file_name)...")

        # Extract row from matrix matching current temperature
        row_idx = findfirst(x -> x == T_celsius, fitted_params[:, 1])
        D_n_opt = fitted_params[row_idx, 2]
        T_amb_opt_C = fitted_params[row_idx, 3]
        
        # Load dataset
        data_file = joinpath(@__DIR__, "..", "data", "BrosaPlanella2021", file_name)
        df = CSV.read(data_file, DataFrame, skipto=18, header=16)

        t_exp = df[!, "Prog Time"]
        v_exp = df[!, "Voltage"]
        T_surf_exp = df[!, "LogTemp001"]
        p_exp = df[!, "Watt"]
        
        # Extract initial ambient temperature from matrix
        T_amb_K = T_amb_opt_C + 273.15

        # Calculate time steps and clean empty rows
        raw_dt = diff(t_exp)
        valid_indices = findall(x -> x > 0.0, raw_dt)
        
        clean_dt = raw_dt[valid_indices]
        clean_p = p_exp[valid_indices]
        clean_T_amb = fill(T_amb_K, length(clean_dt))
        
        exp_duration = t_exp[end]
        soc_init = v_to_soc_interp(v_exp[1])
        
        # Set up parameters and override with pre-fitted values
        p = Chen2020()
        p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
        p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
        # Wrap diffusion parameter in function of concentration to match model architecture
        p.n.Dₖ = c -> D_n_opt
        p.Vmin = 2.0
        p.Vmax = 4.3

        println("Building and compiling mathematical system...")
        # Compile mathematical system per temperature
        @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=15.0, T_ambient=T_amb_K)
        
        # Define experiment using true power array
        step = BatteryToolkit.DriveStep(Any[clean_dt, clean_p, clean_T_amb], exp_duration)
        exp = Experiment([step])

        println("Compiling base ordinary differential equation problem...")
        # Compile base ordinary differential equation problem per temperature
        prob_base = ODEProblem(sys, [sys.Pin => 0.0, sys.Iin => 0.0], (0.0, exp.tend), sparse=true)

        rmse_results = Float64[]
        A_cell = 0.0053
        total_evals = length(h_range)

        for (idx, h_test) in enumerate(h_range)
            # Calculate new thermal resistance
            R_new = 1.0 / (h_test * A_cell)
            
            # Remake problem with new parameter
            prob_test = remake(prob_base, p=[sys.R_conv.R => R_new])
            
            t_start_solve = time()
            
            # Bypass solvers to allow using remade problem directly
            integrator = init(prob_test, QNDF(); tstops=exp.tstops, saveat=10.0, save_everystep=false, reltol=1e-4, abstol=1e-7)
            for step in exp.steps
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
            
            # Filter experimental data to match simulation timeframe
            valid_t_idx = findall(t -> t_sim[1] <= t <= t_sim[end], t_exp)
            t_exp_valid = t_exp[valid_t_idx]
            T_surf_exp_valid = T_surf_exp[valid_t_idx]
            
            # Map simulation to valid experimental timestamps
            T_sim_mapped = sim_temp_interp.(t_exp_valid)
            
            # Calculate root mean square error against valid experimental data
            rmse = sqrt(mean((T_sim_mapped .- T_surf_exp_valid).^2))
            
            # Print formatted output with counter and solve time
            @printf("  [%d/%d] h_conv = %.1f W/m²K | RMSE = %.3f °C | Solve time: %.3fs\n", 
                    idx, total_evals, h_test, rmse, solve_time)
            
            push!(rmse_results, rmse)
        end
        
        all_results[T_celsius] = (h_range, rmse_results)
    end
    
    total_time = time() - t_start_total
    println("\nAll sweeps completed in $(round(total_time, digits=2)) seconds.")
    
    # Plot results
    colors = ["#0072BD", "#D95319", "#EDB120"]
    
    # Initialise plot
    p_loss = plot(xlabel="Convective heat transfer coefficient, h (W/m²K)", 
                  ylabel="Temperature RMSE (°C)",
                  title="Grid search: thermal validation error",
                  legend=:topright, margin=5mm)

    println("\nGrid search complete.")

    # Add curves and optimum points for each temperature
    for (idx, T_celsius) in enumerate(temps)
        local h_range, rmse_data = all_results[T_celsius]
        
        local best_idx = argmin(rmse_data)
        local best_h = h_range[best_idx]
        local best_rmse = rmse_data[best_idx]
        
        println("Optimal convective heat transfer coefficient at $(T_celsius)°C: $best_h W/m²K (RMSE: $(round(best_rmse, digits=3)) °C)")
        
        plot!(p_loss, h_range, rmse_data, marker=:circle, lw=2, color=colors[idx], label="$(T_celsius)°C dataset")
        
        scatter!(p_loss, [best_h], [best_rmse], color=:red, markersize=6, label=false,
                 annotations=(best_h, best_rmse + 0.3, text("Opt: $best_h", :center, :bottom, 10)))
    end

    display(p_loss)
end

# Run computation safely bypassing world age warnings
Base.invokelatest(optimise_h_conv_planella)