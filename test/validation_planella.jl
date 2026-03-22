# =====================================================================================================================
# validation_planella.jl
# Compare SingleCellCoreShellPack simulation to Planella dataset
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using CSV
using DataFrames
using Plots
using Plots.Measures
using DataInterpolations
using Statistics
using Printf
using Logging

using Revise
using BatteryToolkit

Revise.revise()

function validate_planella(target_crate="All")
    # Read voltage to state of charge look up table for initialisation
    lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
    lut = CSV.read(lut_file, DataFrame)
    v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)

    # Define master groups by temperature and convective heat transfer coefficient
    master_temp_groups = [
        (0.0, 13.0, [
            ("Cell786_0p5C_0degC.csv", "0.5C"), 
            ("Cell792_1C_0degC.csv", "1C"),
            ("Cell795_2C_0degC.csv", "2C")
        ]),
        (10.0, 16.0, [
            ("Cell786_0p5C_10degC.csv", "0.5C"), 
            ("Cell792_1C_10degC.csv", "1C"),
            ("Cell795_2C_10degC.csv", "2C")
        ]),
        (25.0, 29.5, [
            ("Cell786_0p5C_25degC.csv", "0.5C"), 
            ("Cell792_1C_25degC.csv", "1C"),
            ("Cell795_2C_25degC.csv", "2C")
        ])
    ]

    # Filter files based on target selection
    temp_groups = []
    for (T_celsius, h_opt, files) in master_temp_groups
        if target_crate == "All"
            filtered_files = files
        else
            filtered_files = filter(f -> f[2] == target_crate, files)
        end
        push!(temp_groups, (T_celsius, h_opt, filtered_files))
    end

    # Dictionary for Planella pre fitted parameters
    fitted_params = Dict(
        0.0 => Dict(
            "0.5C" => (0.22e-14, 0.02),
            "1C"   => (0.55e-14, 0.35),
            "2C"   => (1.50e-14, -0.30)
        ),
        10.0 => Dict(
            "0.5C" => (0.40e-14, 9.80),
            "1C"   => (1.00e-14, 10.10),
            "2C"   => (3.00e-14, 9.60)
        ),
        25.0 => Dict(
            "0.5C" => (0.90e-14, 24.45),
            "1C"   => (2.00e-14, 24.68),
            "2C"   => (6.00e-14, 24.30)
        )
    )

    test_labels = String[]
    
    # Calculate total required rows for data matrix
    total_simulations = sum(length(files) for (_, _, files) in temp_groups)
    data_rows = Matrix{Any}(undef, total_simulations, 4)
    solve_times = Float64[]
    row_idx = 1
    
    # Array to hold subplots
    master_plots = []
    
    # Define unique colours for temperature groups
    temp_colours = ["#0072BD", "#D95319", "#EDB120"]

    for (temp_idx, (T_celsius, h_opt, files)) in enumerate(temp_groups)
        println("\nProcessing $(T_celsius)°C group...")

        # Initialise plot subplots per temperature forcing left margin to prevent clipping
        m = 7mm
        p_volt = plot(ylabel="Voltage (V)", legend=:bottomleft, margin=m, left_margin=15mm)
        p_temp = plot(xlabel="Time (s)", ylabel="Temperature (°C)", legend=:topleft, margin=m, left_margin=15mm)
        
        # Select colour for current temperature group
        current_colour = temp_colours[temp_idx]
        
        # Arrays to store group metrics for subplot title
        group_rmse_v = Float64[]
        group_rmse_T = Float64[]

        for (file_name, Crate) in files
            label_base = "$(T_celsius)°C, $(Crate)"
            println("  Processing $(label_base)...")
            push!(test_labels, label_base)

            # Load dataset setting header and skipping units
            data_file = joinpath(@__DIR__, "..", "data", "BrosaPlanella2021", file_name)
            df = CSV.read(data_file, DataFrame, skipto=18, header=16)

            # Extract raw columns
            t_exp = df[!, "Prog Time"]
            v_exp = df[!, "Voltage"]
            T_surf_exp = df[!, "LogTemp001"]
            p_exp = df[!, "Watt"]
            status_col = hasproperty(df, :status) ? df[!, "status"] : df[!, "Status"]

            # Identify first segment starting exactly at discharge
            start_idx = 1
            for i in 2:nrow(df)
                prev_stat = strip(String(status_col[i-1]))
                curr_stat = strip(String(status_col[i]))
                if prev_stat != "DCH" && curr_stat == "DCH"
                    start_idx = i - 1
                    break
                end
            end

            end_idx = nrow(df)
            for i in (start_idx + 1):nrow(df)
                curr_stat = strip(String(status_col[i]))
                if curr_stat == "CHA"
                    end_idx = i - 1
                    break
                end
            end

            # Extract segment data
            t_seg = t_exp[start_idx:end_idx]
            v_seg = v_exp[start_idx:end_idx]
            T_surf_seg = T_surf_exp[start_idx:end_idx]
            p_seg = p_exp[start_idx:end_idx]

            # Local time vector starting from zero for simulator and plotting
            t_sim_input = t_seg .- t_seg[1]

            # Extract parameters from dictionary matching current temperature and rate
            D_n_opt, T_amb_opt = fitted_params[T_celsius][Crate]

            # Extract initial ambient temperature from dictionary or CSV
            if isnothing(T_amb_opt)
                T_amb_K = T_surf_seg[1] + 273.15
            else
                T_amb_K = T_amb_opt + 273.15
            end

            # Read initial voltage from segment to calculate state of charge
            v_init = v_seg[1]
            soc_init = v_to_soc_interp(v_init)

            # Set up parameters
            p = Chen2020()
            p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
            p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
            
            # Wrap diffusion parameter in function of concentration to match model architecture
            p.n.Dₖ = c -> D_n_opt
            
            # Loosen voltage constraints to prevent solver crashing at high rates
            p.Vmin = 1.0
            p.Vmax = 5.5

            println("    Building and compiling MTK system...")
            
            # Build MTK system using optimised heat transfer coefficient and specific parameters
            @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=h_opt, T_ambient=T_amb_K)

            # Calculate time steps and clean empty time rows
            raw_dt = diff(t_sim_input)
            valid_indices = findall(x -> x > 0.0, raw_dt)
            
            clean_dt = raw_dt[valid_indices]
            clean_p = p_seg[valid_indices]
            clean_T_amb = fill(T_amb_K, length(clean_dt))
            
            exp_duration = t_sim_input[end]

            # Define experiment using cleaned data array
            step = BatteryToolkit.DriveStep(Any[clean_dt, clean_p, clean_T_amb], exp_duration)
            exp = Experiment([step])

            # Run simulation and capture solve time while silencing internal package warnings
            t_start = time()
            local sol = with_logger(NullLogger()) do
                simulate(sys, exp; saveat=10.0, verbose=false)
            end
            push!(solve_times, time() - t_start)

            # Extract results directly from zero based solution
            t_sim = sol.t
            v_sim = sol[sys.cell.v]
            T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15

            # Interpolate simulated curves to match experimental timestamps
            v_sim_interp = LinearInterpolation(v_sim, t_sim)
            T_sim_interp = LinearInterpolation(T_shell_sim, t_sim)
            
            # Match lengths in case simulation terminates early
            valid_t_idx = findall(t -> t_sim[1] <= t <= t_sim[end], t_sim_input)
            t_eval_valid = t_sim_input[valid_t_idx]
            v_seg_valid = v_seg[valid_t_idx]
            T_surf_seg_valid = T_surf_seg[valid_t_idx]

            # Map simulation to zero based segment experimental timestamps
            v_sim_mapped = v_sim_interp.(t_eval_valid)
            T_sim_mapped = T_sim_interp.(t_eval_valid)

            # Calculate coefficient of determination helper function
            function calculate_r2(exp_data, sim_mapped)
                ss_res = sum((exp_data .- sim_mapped).^2)
                ss_tot = sum((exp_data .- mean(exp_data)).^2)
                return ss_tot == 0 ? 0.0 : 1.0 - (ss_res / ss_tot)
            end

            # Calculate voltage metrics
            rmse_v_mv = sqrt(mean((v_sim_mapped .- v_seg_valid).^2)) * 1000
            r2_v = calculate_r2(v_seg_valid, v_sim_mapped)
            
            # Calculate temperature metrics
            rmse_T = sqrt(mean((T_sim_mapped .- T_surf_seg_valid).^2))
            r2_T = calculate_r2(T_surf_seg_valid, T_sim_mapped)

            # Store in matrix
            data_rows[row_idx, 1] = round(rmse_v_mv, digits=2)
            data_rows[row_idx, 2] = round(r2_v, digits=4)
            data_rows[row_idx, 3] = round(rmse_T, digits=2)
            data_rows[row_idx, 4] = round(r2_T, digits=4)
            row_idx += 1
            
            # Store in group arrays for subplot title
            push!(group_rmse_v, rmse_v_mv)
            push!(group_rmse_T, rmse_T)

            # Generate labels
            label_exp = "Exp. " * Crate
            label_sim = "Sim. " * Crate
            
            # Add to plots using unified colour for temperature group plotting from zero
            plot!(p_volt, t_sim_input, v_seg, label=label_exp, lw=2, color=:black, linestyle=:dash)
            plot!(p_volt, t_sim, v_sim, label=label_sim, lw=2, color=current_colour, linestyle=:solid)

            plot!(p_temp, t_sim_input, T_surf_seg, label=label_exp, lw=2, color=:black, linestyle=:dash)
            plot!(p_temp, t_sim, T_shell_sim, label=label_sim, lw=2, color=current_colour, linestyle=:solid)
        end
        
        # Calculate average root mean square error for temperature group
        avg_rmse_v = mean(group_rmse_v)
        avg_rmse_T = mean(group_rmse_T)
        
        # Add informative title to top subplot
        title!(p_volt, "$(T_celsius)°C\nAvg RMSE: $(round(avg_rmse_v, digits=1)) mV | $(round(avg_rmse_T, digits=2)) °C")

        # Combine plot layouts for temperature group and add to array
        l = @layout [a; b]
        p_combined = plot(p_volt, p_temp, layout=l)
        push!(master_plots, p_combined)
    end

    # Combine all temperature groups into master figure with unified title
    p_master = plot(master_plots..., layout=(1, length(master_plots)), size=(600 * length(master_plots), 750), 
                    margin=5mm, plot_title="Voltage and temperature validation against Planella2021", plot_titlefontsize=20)
    display(p_master)

    # Output final metrics table
    println("\nFinal Model Performance Metrics")
    
    # Print formatted headers
    @printf("%-18s | %-18s | %-10s | %-14s | %-10s\n", 
            "Test", "Voltage RMSE (mV)", "Voltage R²", "Temp RMSE (°C)", "Temp R²")

    # Print data rows
    for i in 1:size(data_rows, 1)
        @printf("%-18s | %-18.2f | %-10.4f | %-14.2f | %-10.4f\n", 
                test_labels[i], data_rows[i, 1], data_rows[i, 2], data_rows[i, 3], data_rows[i, 4])
    end

    println("Average solve time: $(round(mean(solve_times), digits=2)) seconds\n")
end

# Run validation sequence safely
# Pass "0.5C", "1C", or "2C" to run a specific rate
Base.invokelatest(validate_planella, "0.5C")