# ==============================================================================
# h_optimisation_planella.jl
# Performs 1D grid search to find optimal convective heat transfer coefficient 
# for Planella 2021 datasets at 0, 10, and 25°C
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

# Generates data for all temperatures
function optimise_h_conv_planella()
    temps = [25]
    files = ["Cell789_1C_25degC.csv"]
    
    # Read V-SoC LUT for initialisation
    lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
    lut = CSV.read(lut_file, DataFrame)
    v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)

    h_range = 30.0:1:40.0
    all_rmse_results = Dict{Int, Vector{Float64}}()

    println("Starting grid search for Planella datasets...")

    for (T_celsius, file_name) in zip(temps, files)
        println("\nProcessing $(T_celsius)°C dataset...")
        
        # Load dataset
        data_file = joinpath(@__DIR__, "..", "data", "BrosaPlanella2021", file_name)
        df = CSV.read(data_file, DataFrame, skipto=18, header=16)

        t_exp = df[!, "Prog Time"]
        v_exp = df[!, "Voltage"]
        T_surf_exp = df[!, "LogTempMid"]
        raw_i_exp = df[!, "Current"]
        status = df[!, "Status"]
        
        T_amb_K = T_celsius + 273.15

        # Process current directions
        i_exp = zeros(Float64, length(raw_i_exp))
        for j in eachindex(raw_i_exp)
            mode_str = strip(String(status[j]))
            if mode_str == "CHA"
                i_exp[j] = -abs(raw_i_exp[j])
            elseif mode_str == "DCH"
                i_exp[j] = abs(raw_i_exp[j])
            else
                i_exp[j] = 0.0
            end
        end

        # Calculate time steps and clean zero dt rows
        raw_dt = diff(t_exp)
        valid_indices = findall(x -> x > 0.0, raw_dt)
        
        clean_dt = raw_dt[valid_indices]
        clean_i = i_exp[valid_indices]
        clean_T_amb = fill(T_amb_K, length(clean_dt))
        
        exp_duration = t_exp[end]
        soc_init = v_to_soc_interp(v_exp[1])
        
        p = Chen2020()
        p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
        p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
        p.Vmin = 2.0
        p.Vmax = 4.3

        step = BatteryPeck.CurrentDriveStep(Any[clean_dt, clean_i, clean_T_amb], exp_duration)
        exp = Experiment([step])

        rmse_results = Float64[]

        for h_test in h_range
            println("Evaluating h_conv = $h_test W/m²K ...")
            
            @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=h_test, T_ambient=T_amb_K)
            sol = simulate(sys, exp, saveat=10.0)
            
            t_sim = sol.t
            T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15
            
            sim_temp_interp = LinearInterpolation(T_shell_sim, t_sim)
            T_sim_mapped = sim_temp_interp.(clamp.(t_exp, t_sim[1], t_sim[end]))
            
            rmse = sqrt(mean((T_sim_mapped .- T_surf_exp).^2))
            push!(rmse_results, rmse)
        end
        
        all_rmse_results[T_celsius] = rmse_results
    end
    
    return collect(h_range), all_rmse_results
end

# Run computation once to obtain data
h_range, all_rmse_results = optimise_h_conv_planella()

# Plot results
colors = ["#0072BD", "#D95319", "#EDB120"]
temps = [25]

# Initialise plot
p_loss = plot(xlabel="Convective heat transfer coefficient, h (W/m²K)", 
              ylabel="Temperature RMSE (°C)",
              title="1D grid search: thermal validation error (Planella)",
              legend=:topright, margin=5mm)

println("\nGrid search complete.")

# Add curves and optimum points for each temperature
for (idx, T_celsius) in enumerate(temps)
    rmse_data = all_rmse_results[T_celsius]
    
    best_idx = argmin(rmse_data)
    best_h = h_range[best_idx]
    best_rmse = rmse_data[best_idx]
    
    println("Optimal h_conv at $(T_celsius)°C: $best_h W/m²K (RMSE: $(round(best_rmse, digits=3)) °C)")
    
    plot!(p_loss, h_range, rmse_data, marker=:circle, lw=2, color=colors[idx], label="$(T_celsius)°C dataset")
    
    scatter!(p_loss, [best_h], [best_rmse], color=:red, markersize=6, label=false,
             annotations=(best_h, best_rmse + 0.5, text("Opt: $best_h", :center, :bottom, 10)))
end

display(p_loss)