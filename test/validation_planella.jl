# =====================================================================================================================
# validation_planella.jl
# Compare SingleCellCoreShellPack simulation to Planella 2021 dataset at various temperatures and C-rates
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

using Revise
using BatteryPeck

Revise.revise()

function validate_planella()
    # Initialise plot subplots
    m = 7mm
    p_volt = plot(ylabel="Voltage (V)", legend=:bottomleft, margin=m)
    p_temp = plot(ylabel="Temperature (°C)", legend=:topleft, margin=m)
    p_curr = plot(xlabel="Time (s)", ylabel="Current (A)", legend=:bottomleft, margin=m)

    # Define full test matrix with temperature-specific h_conv values
    files = [
        ("Cell781_0p1C_0degC.csv",   0.0, "0.1C", 20.5),
        ("Cell781_0p1C_10degC.csv", 10.0, "0.1C", 27.0),
        ("Cell781_0p1C_25degC.csv", 25.0, "0.1C", 30.0),
        ("Cell789_1C_0degC.csv",     0.0, "1C",   20.5),
        ("Cell789_1C_10degC.csv",   10.0, "1C",   27.0),
        ("Cell789_1C_25degC.csv",   25.0, "1C",   30.0),
        ("Cell793_2C_0degC.csv",     0.0, "2C",   20.5),
        ("Cell793_2C_10degC.csv",   10.0, "2C",   27.0),
        ("Cell793_2C_25degC.csv",   25.0, "2C",   30.0)
    ]

    # Initialise metrics collection
    test_labels = String[]
    data_rows = Matrix{Any}(undef, length(files), 4) 
    solve_times = Float64[]

    # Read V-SoC LUT for initialisation
    lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
    lut = CSV.read(lut_file, DataFrame)
    v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)

    for (idx, (file_name, T_celsius, Crate, h_opt)) in enumerate(files)
        label_base = "$(T_celsius)°C, $(Crate)"
        println("\nProcessing $(label_base)...")
        push!(test_labels, label_base)
        
        # Load dataset setting header to row 16 and skipping units on row 17
        data_file = joinpath(@__DIR__, "..", "data", "BrosaPlanella2021", file_name)
        df = CSV.read(data_file, DataFrame, skipto=18, header=16)

        # Extract raw columns
        t_exp = df[!, "Prog Time"]
        v_exp = df[!, "Voltage"]
        T_surf_exp = df[!, "LogTempMid"]
        raw_i_exp = df[!, "Current"]
        status = df[!, "Status"]
        
        T_amb_K = T_celsius + 273.15

        # Process current directions based on status string
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

        # Dynamically calculate soc_init
        soc_init = v_to_soc_interp(v_exp[1])
        
        # Setup parameters
        p = Chen2020()
        p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
        p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
        
        # Loosen voltage constraints to prevent solver crashing at 2C
        # High resistance at 0°C and 10°C causes huge voltage spikes during charging phase
        p.Vmin = 1.0
        p.Vmax = 5.5

        # Build system using optimised h coefficient
        @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=h_opt, T_ambient=T_amb_K)
        
        # Define experiment using cleaned data array
        step = BatteryPeck.CurrentDriveStep(Any[clean_dt, clean_i, clean_T_amb], exp_duration)
        exp = Experiment([step])

        # Run simulation and capture solve time
        println("  Running simulation...")
        t_start = time()
        sol = simulate(sys, exp, saveat=10.0)
        push!(solve_times, time() - t_start)

        # Extract results
        t_sim = sol.t
        v_sim = sol[sys.cell.v]
        i_sim = sol[sys.I]
        T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15

        # Interpolate simulated curves to match experimental timestamps
        v_sim_interp = LinearInterpolation(v_sim, t_sim)
        T_sim_interp = LinearInterpolation(T_shell_sim, t_sim)
        
        # Map simulation to exact experimental timestamps using clamped times
        t_clamped = clamp.(t_exp, t_sim[1], t_sim[end])
        v_sim_mapped = v_sim_interp.(t_clamped)
        T_sim_mapped = T_sim_interp.(t_clamped)

        # Calculate R2 helper function
        function calculate_r2(exp_data, sim_mapped)
            ss_res = sum((exp_data .- sim_mapped).^2)
            ss_tot = sum((exp_data .- mean(exp_data)).^2)
            return 1.0 - (ss_res / ss_tot)
        end

        # Calculate voltage metrics
        rmse_v_mv = sqrt(mean((v_sim_mapped .- v_exp).^2)) * 1000
        r2_v = calculate_r2(v_exp, v_sim_mapped)
        
        # Calculate temperature metrics
        rmse_T = sqrt(mean((T_sim_mapped .- T_surf_exp).^2))
        r2_T = calculate_r2(T_surf_exp, T_sim_mapped)

        # Store in matrix
        data_rows[idx, 1] = round(rmse_v_mv, digits=2)
        data_rows[idx, 2] = round(r2_v, digits=4)
        data_rows[idx, 3] = round(rmse_T, digits=2)
        data_rows[idx, 4] = round(r2_T, digits=4)

        # Generate labels
        label_exp = "Exp. " * label_base
        label_sim = "Sim. " * label_base
        
        # Assign color cycle based on idx 
        plot_c = idx 

        # Add to plots
        plot!(p_volt, t_exp, v_exp, label=label_exp, lw=2, color=plot_c, linestyle=:dash)
        plot!(p_volt, t_sim, v_sim, label=label_sim, lw=2, color=plot_c, linestyle=:solid)

        plot!(p_temp, t_exp, T_surf_exp, label=label_exp, lw=2, color=plot_c, linestyle=:dash)
        plot!(p_temp, t_sim, T_shell_sim, label=label_sim, lw=2, color=plot_c, linestyle=:solid)

        plot!(p_curr, t_exp, i_exp, label=label_exp, lw=2, color=plot_c, linestyle=:dash)
        plot!(p_curr, t_sim, i_sim, label=label_sim, lw=2, color=plot_c, linestyle=:solid)
    end

    # Output final metrics table
    println("\n" * "="^78)
    println("Final Model Performance Metrics")
    println("-"^78)
    
    # Print formatted headers
    @printf("%-18s | %-18s | %-10s | %-14s | %-10s\n", 
            "Test", "Voltage RMSE (mV)", "Voltage R²", "Temp RMSE (°C)", "Temp R²")
    println("-"^78)

    # Print data rows
    for i in 1:size(data_rows, 1)
        @printf("%-18s | %-18.2f | %-10.4f | %-14.2f | %-10.4f\n", 
                test_labels[i], data_rows[i, 1], data_rows[i, 2], data_rows[i, 3], data_rows[i, 4])
    end

    println("-"^78)
    println("Average Solve Time: $(round(mean(solve_times), digits=2)) seconds")
    println("="^78 * "\n")

    # Combine plot layouts
    l = @layout [a; b; c]
    p_combined = plot(p_volt, p_temp, p_curr, layout=l, size=(1000, 1000), plot_title="Planella 2021 multi-rate thermal validation")
    display(p_combined)
end

# Run validation sequence
validate_planella()