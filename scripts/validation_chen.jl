# =====================================================================================================================
# validation_chen.jl
# Compare SingleCellCoreShellPack simulation to experimental Chen dataset across 4 segments
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using CSV
using DataFrames
using Plots
using Plots.Measures
using DataInterpolations
using Logging
using Statistics
using Printf

using Revise
using BatteryToolkit

Revise.revise()

# Set toggle to isolate and simulate R to D to R segments sequentially
simulate_segments_only = true

# Load experimental data
data_file = joinpath(@__DIR__, "..", "data", "Chen2020", "LGM50_cell03.csv")

# Skip metadata and load dataframe to extract plotting arrays
df = CSV.read(data_file, DataFrame, skipto=15, header=14)

t_exp = df[!, "Test Time [s]"]
v_exp = df[!, "Voltage [V]"]
T_surf_exp = df[!, "Temperature Cell [degC]"]
T_chamber_exp = df[!, "Temperature Chamber [degC]"]
T_amb_K = T_chamber_exp[1] + 273.15

# Read voltage SOC lookup table
lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
lut = CSV.read(lut_file, DataFrame)

# Create interpolator mapping voltage to SOC
v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)

# Process experimental current for plotting
raw_i_exp = df[!, "Current [A]"]
md_exp = df[!, "Md"]
i_exp = zeros(Float64, length(raw_i_exp))

for j in eachindex(raw_i_exp)
    mode_str = strip(String(md_exp[j]))
    if mode_str == "C"
        i_exp[j] = -raw_i_exp[j]
    elseif mode_str == "D"
        i_exp[j] = raw_i_exp[j]
    else
        i_exp[j] = 0.0
    end
end

# Extract isolated R to D to R segments from dataframe
function extract_discharge_segments(df)
    segments = DataFrame[]
    n = nrow(df)
    
    i = 1
    while i <= n
        if strip(df.Md[i]) == "D"
            # Find start of preceding R block
            start_idx = i
            while start_idx > 1 && strip(df.Md[start_idx - 1]) == "R"
                start_idx -= 1
            end
            
            # Find end of D block
            end_idx = i
            while end_idx < n && strip(df.Md[end_idx + 1]) == "D"
                end_idx += 1
            end
            
            # Find end of succeeding R block
            while end_idx < n && strip(df.Md[end_idx + 1]) == "R"
                end_idx += 1
            end
            
            push!(segments, df[start_idx:end_idx, :])
            i = end_idx
        end
        i += 1
    end
    
    # Ignore segments that do not start at full charge
    filter!(seg -> seg[1, "Voltage [V]"] >= 4.1, segments)
    
    return segments
end

# Load parameter set
p = Chen2020()

p.Vmin = 1.0
p.Vmax = 5.0

m = 7mm
xlim_range = (0, t_exp[end])
colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]

let
    if !simulate_segments_only
        # Initialise plot layouts with full CSV traces
        p_volt = plot(t_exp, v_exp, label="Exp. voltage", ylabel="Voltage (V)", 
                    lw=2, color=:black, linestyle=:dash, legend=:bottomleft, margin=m, xlims=xlim_range)

        p_temp = plot(t_exp, T_surf_exp, label="Exp. surface temp", ylabel="Temperature (°C)", 
                    lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)
        plot!(p_temp, t_exp, T_chamber_exp, label="Exp. chamber temp", lw=2, color=:gray, linestyle=:dot)

        p_curr = plot(t_exp, i_exp, label="Exp. current", xlabel="Time (s)", ylabel="Current (A)", 
                    lw=2, color=:black, linestyle=:dash, legend=:topleft, margin=m, xlims=xlim_range)

        # Local time vector starting from zero for simulator
        t_sim_input = t_exp .- t_exp[1]

        # Initialise SOC based on very first starting voltage
        v_init = v_exp[1]
        soc_init = v_to_soc_interp(v_init) 

        # Map SOC to stoichiometries
        z_n_init = p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)
        z_p_init = p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)

        # Set initial solid concentrations
        p.n.c₀ = z_n_init * p.n.c₊  
        p.p.c₀ = z_p_init * p.p.c₊

        # Build coupled core shell system
        @mtkbuild sys = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=19.5, T_ambient=T_amb_K)

        # Format data for current drive step
        raw_dt = diff(t_sim_input)
        valid_indices = findall(x -> x > 0.0, raw_dt)
        clean_dt = raw_dt[valid_indices]
        clean_i = i_exp[valid_indices]
        clean_T_amb = T_chamber_exp[valid_indices .+ 1] .+ 273.15

        # Define experiment for entire dataset
        exp_full = Experiment([
            CurrentDriveStep(Any[clean_dt, clean_i, clean_T_amb], t_sim_input[end])
        ])

        println("Running continuous validation simulation...\n")

        # Run simulation
        sol = with_logger(NullLogger()) do
            simulate(sys, exp_full; saveat=10.0, verbose=false)
        end

        # Extract simulation results and shift time back to global timeframe
        t_sim = sol.t .+ t_exp[1]
        v_sim = sol[sys.cell.v]
        i_sim = sol[sys.I]
        T_shell_sim = sol[sys.thermal.shell_cap.T] .- 273.15

        # Interpolate for metrics calculation against experimental data points
        v_sim_interp = LinearInterpolation(v_sim, t_sim)
        T_shell_sim_interp = LinearInterpolation(T_shell_sim, t_sim)

        t_clamped = clamp.(t_exp, t_sim[1], t_sim[end])
        v_sim_mapped = v_sim_interp.(t_clamped)
        T_sim_mapped = T_shell_sim_interp.(t_clamped)

        # Calculate RMSE
        rmse_v = sqrt(mean((v_sim_mapped .- v_exp).^2))
        rmse_T = sqrt(mean((T_sim_mapped .- T_surf_exp).^2))

        # Calculate R squared
        ss_res_v = sum((v_exp .- v_sim_mapped).^2)
        ss_tot_v = sum((v_exp .- mean(v_exp)).^2)
        r2_v = ss_tot_v == 0 ? 0.0 : 1.0 - (ss_res_v / ss_tot_v)

        ss_res_T = sum((T_surf_exp .- T_sim_mapped).^2)
        ss_tot_T = sum((T_surf_exp .- mean(T_surf_exp)).^2)
        r2_T = ss_tot_T == 0 ? 0.0 : 1.0 - (ss_res_T / ss_tot_T)

        # Print overall metrics
        println("Continuous Simulation Results:")
        @printf("  Starting voltage | %.4f V\n", v_init)
        @printf("  Voltage          | RMSE: %7.2f mV | R²: %.4f\n", rmse_v * 1000, r2_v)
        @printf("  Temperature      | RMSE: %7.2f °C | R²: %.4f\n\n", rmse_T, r2_T)

        # Plot continuous traces
        plot!(p_volt, t_sim, v_sim, label="Sim. voltage", lw=2, color=colours[1])
        plot!(p_temp, t_sim, T_shell_sim, label="Sim. surface temp", lw=2, color=colours[2])
        plot!(p_curr, t_sim, i_sim, label="Sim. current", lw=2, color=colours[3])

    else
        # Initialise empty plots for segmented traces
        p_volt = plot(ylabel="Voltage (V)", legend=:bottomleft, margin=m, xlims=xlim_range)
        p_temp = plot(ylabel="Temperature (°C)", legend=:topleft, margin=m, xlims=xlim_range)
        p_curr = plot(xlabel="Time (s)", ylabel="Current (A)", legend=:topleft, margin=m, xlims=xlim_range)

        # Extract discharge segments from dataframe
        segments = extract_discharge_segments(df)
        println("Found $(length(segments)) discharge segments to simulate.\n")

        for (idx, seg) in enumerate(segments)
            println("Simulating segment $idx...")

            # Extract segment time and convert to local time starting from zero
            t_seg_exp = seg[!, "Test Time [s]"]
            t_sim_input = t_seg_exp .- t_seg_exp[1]

            # Extract experimental variables for segment
            v_seg_exp = seg[!, "Voltage [V]"]
            T_chamber_seg = seg[!, "Temperature Chamber [degC]"]
            T_surf_seg = seg[!, "Temperature Cell [degC]"]
            raw_i_seg = seg[!, "Current [A]"]
            md_seg = seg[!, "Md"]

            # Process segment current
            i_seg = zeros(Float64, length(raw_i_seg))
            for j in eachindex(raw_i_seg)
                mode_str = strip(String(md_seg[j]))
                if mode_str == "C"
                    i_seg[j] = -raw_i_seg[j]
                elseif mode_str == "D"
                    i_seg[j] = raw_i_seg[j]
                else
                    i_seg[j] = 0.0
                end
            end

            # Read starting voltage to set accurate initial SOC
            v_init = v_seg_exp[1]
            soc_init = v_to_soc_interp(v_init)
            
            # Map SOC to stoichiometries
            z_n_init = p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)
            z_p_init = p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)

            # Set initial solid concentrations
            p.n.c₀ = z_n_init * p.n.c₊  
            p.p.c₀ = z_p_init * p.p.c₊

            # Build coupled core shell system with segment parameters
            T_amb_seg_K = T_chamber_seg[1] + 273.15
            @mtkbuild sys_seg = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=19.5, T_ambient=T_amb_seg_K)

            # Format data for current drive step
            raw_dt = diff(t_sim_input)
            valid_indices = findall(x -> x > 0.0, raw_dt)
            clean_dt = raw_dt[valid_indices]
            clean_i = i_seg[valid_indices]
            clean_T_amb = T_chamber_seg[valid_indices .+ 1] .+ 273.15 

            # Define experiment for segment
            exp_seg = Experiment([
                CurrentDriveStep(Any[clean_dt, clean_i, clean_T_amb], t_sim_input[end])
            ])

            # Run simulation
            sol = with_logger(NullLogger()) do
                simulate(sys_seg, exp_seg; saveat=10.0, verbose=false)
            end

            # Verify solver stability and calculate metrics
            if sol.t[end] > 10.0
                # Shift time back to global timeframe
                t_shift = t_seg_exp[1]
                t_sim = sol.t .+ t_shift
                v_sim = sol[sys_seg.cell.v]
                i_sim = sol[sys_seg.I]
                T_shell_sim = sol[sys_seg.thermal.shell_cap.T] .- 273.15

                # Calculate metrics for segment
                v_sim_interp = LinearInterpolation(v_sim, t_sim)
                T_shell_sim_interp = LinearInterpolation(T_shell_sim, t_sim)

                t_clamped = clamp.(t_seg_exp, t_sim[1], t_sim[end])
                v_sim_mapped = v_sim_interp.(t_clamped)
                T_sim_mapped = T_shell_sim_interp.(t_clamped)

                rmse_v = sqrt(mean((v_sim_mapped .- v_seg_exp).^2))
                rmse_T = sqrt(mean((T_sim_mapped .- T_surf_seg).^2))

                ss_res_v = sum((v_seg_exp .- v_sim_mapped).^2)
                ss_tot_v = sum((v_seg_exp .- mean(v_seg_exp)).^2)
                r2_v = ss_tot_v == 0 ? 0.0 : 1.0 - (ss_res_v / ss_tot_v)

                ss_res_T = sum((T_surf_seg .- T_sim_mapped).^2)
                ss_tot_T = sum((T_surf_seg .- mean(T_surf_seg)).^2)
                r2_T = ss_tot_T == 0 ? 0.0 : 1.0 - (ss_res_T / ss_tot_T)

                # Print segment metrics
                println("Segment $idx Results:")
                @printf("  Starting voltage | %.4f V\n", v_init)
                @printf("  Voltage          | RMSE: %7.2f mV | R²: %.4f\n", rmse_v * 1000, r2_v)
                @printf("  Temperature      | RMSE: %7.2f °C | R²: %.4f\n\n", rmse_T, r2_T)

                # Plot experimental segment traces
                plot!(p_volt, t_seg_exp, v_seg_exp, label=idx==1 ? "Exp. voltage" : "", lw=2, color=:black, linestyle=:dash)
                plot!(p_temp, t_seg_exp, T_surf_seg, label=idx==1 ? "Exp. surface temp" : "", lw=2, color=:black, linestyle=:dash)
                plot!(p_temp, t_seg_exp, T_chamber_seg, label=idx==1 ? "Exp. chamber temp" : "", lw=2, color=:gray, linestyle=:dot)
                plot!(p_curr, t_seg_exp, i_seg, label=idx==1 ? "Exp. current" : "", lw=2, color=:black, linestyle=:dash)

                # Plot simulation segment traces using a single fixed colour
                plot!(p_volt, t_sim, v_sim, label=idx==1 ? "Sim. voltage" : "", lw=2, color=colours[1])
                scatter!(p_volt, [t_shift], [v_init], shape=:cross, markersize=6, color=:black, label="")

                plot!(p_temp, t_sim, T_shell_sim, label=idx==1 ? "Sim. surface temp" : "", lw=2, color=colours[2])
                plot!(p_curr, t_sim, i_sim, label=idx==1 ? "Sim. current" : "", lw=2, color=colours[3])
            else
                println("Warning: Simulation $idx aborted early.\n")
            end
        end
    end

    # Combine layouts
    l = @layout [a; b; c]
    p_combined = plot(p_volt, p_temp, p_curr, layout=l, size=(1000, 1000), plot_title="Coupled core-shell thermal SPMe validation")
    
    # Show the plot in the VS Code Plot Pane (defaults to PNG for speed)
    display(p_combined)

    # Set toggle to automatically export a vector graphics file of the final plot
    export_vector_image = true

    # Automatically save a vector version if the toggle is true
    if export_vector_image
        output_filename = "validation_chen_plot.svg" 
        savefig(p_combined, output_filename)
        println("Vector image successfully exported to: $output_filename")
    end
end