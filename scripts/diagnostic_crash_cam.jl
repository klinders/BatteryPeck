# ==============================================================================
# diagnostic_crash_cam.jl
# Full-Pack V2G Autopsy: 90-Day Adaptive BMS Learning Loop.
# ==============================================================================

using Plots
using Revise
using BatteryToolkit
using CSV, DataFrames
using Serialization
using SciMLBase
using Printf

Revise.revise()

const CHECKPOINT_5_YEAR = joinpath(pwd(), "aging_5_years_scn_5.jls")
const CHECKPOINT_V2G_MORNING = joinpath(pwd(), "v2g_morning_checkpoint.jls")

function save_checkpoint(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)
    data = Dict(
        "therm_u" => cosim.therm_integrator.u,
        "therm_p" => cosim.therm_integrator.p,
        "therm_t" => cosim.therm_integrator.t,
        "cells_u" => [int.u for int in cosim.cell_integrators],
        "cells_p" => [int.p for int in cosim.cell_integrators],
        "cells_t" => [int.t for int in cosim.cell_integrators]
    )
    serialize(filepath, data)
    if !silent println("[!] Checkpoint successfully saved to: $filepath") end
end

function load_checkpoint!(cosim::ExplicitPackSimulator, filepath::String; silent::Bool=false)
    data = deserialize(filepath)
    SciMLBase.reinit!(cosim.therm_integrator, data["therm_u"]; t0=data["therm_t"], reset_dt=true)
    cosim.therm_integrator.p = data["therm_p"]

    for i in 1:cosim.num_cells
        SciMLBase.reinit!(cosim.cell_integrators[i], data["cells_u"][i]; t0=data["cells_t"][i], reset_dt=true)
        cosim.cell_integrators[i].p = data["cells_p"][i]
    end
    if !silent println("[!] Checkpoint successfully loaded into CoSim architecture.") end
end

function balance_pack!(cosim::ExplicitPackSimulator, target_soc::Float64, max_duration_hrs::Float64)
    println("\n>>> [PHASE 1.5] INITIATING WORKSHOP ACTIVE BALANCING...")
    
    num_cells = cosim.num_cells
    dt = 60.0 
    duration = max_duration_hrs * 3600.0
    t_bal = 0.0
    start_t_abs = cosim.therm_integrator.t
    
    K_prop = 2.0           
    max_bal_current = 0.2  
    
    sys = cosim.sys_cell_template
    socs = zeros(num_cells)
    
    I_prev = zeros(num_cells)
    for i in 1:num_cells
        I_prev[i] = cosim.cell_integrators[i].ps[cosim.I_sym]
    end
    
    any_crashed = false

    while t_bal < duration - 1e-6
        for i in 1:num_cells
            socs[i] = cosim.cell_integrators[i][sys.cell.soc]
        end
        
        if mod(t_bal, 3600.0) == 0.0 || t_bal == 0.0
            min_soc = minimum(socs) * 100.0
            max_soc = maximum(socs) * 100.0
            spread = max_soc - min_soc
            @printf("\r    -> Time: %5.1f h / %5.1f h | Target: %.1f%% | Spread: %5.3f%% | Min: %5.2f%% | Max: %5.2f%%        ", 
                    t_bal/3600.0, max_duration_hrs, target_soc*100.0, spread, min_soc, max_soc)
            flush(stdout)
            
            if spread < 0.05 && abs(max_soc - target_soc*100.0) < 0.05
                print("\n    -> Pack fully synchronized early!\n")
                break
            end
        end
        
        if t_bal + dt > duration
            dt = duration - t_bal
        end
        
        any_crashed = false
        
        for i in 1:num_cells
            int_c = cosim.cell_integrators[i]
            if int_c.sol.retcode != SciMLBase.ReturnCode.Success && int_c.sol.retcode != SciMLBase.ReturnCode.Default
                continue
            end
            
            err = target_soc - socs[i]
            I_target = clamp(K_prop * err, -max_bal_current, max_bal_current)
            I_curr = I_prev[i]
            dI_dt = (I_target - I_curr) / dt
            
            int_c.ps[cosim.I_sym] = I_curr
            int_c.ps[cosim.dI_dt_sym] = dI_dt
            int_c.ps[cosim.t_start_sym] = start_t_abs + t_bal
            
            t_target = start_t_abs + t_bal + dt
            SciMLBase.add_tstop!(int_c, t_target)
            while int_c.t < t_target
                SciMLBase.step!(int_c)
                if int_c.sol.retcode != SciMLBase.ReturnCode.Success && int_c.sol.retcode != SciMLBase.ReturnCode.Default
                    any_crashed = true
                    break
                end
            end
            I_prev[i] = I_target
        end
        
        if any_crashed break end
        
        for i in 1:num_cells
            cosim.therm_integrator.ps[cosim.Q_syms[i]] = 0.0
            cosim.therm_integrator.ps[cosim.dQ_dt_syms[i]] = 0.0
        end
        cosim.therm_integrator.ps[cosim.t_start_therm_sym] = start_t_abs + t_bal
        
        t_target_therm = start_t_abs + t_bal + dt
        SciMLBase.add_tstop!(cosim.therm_integrator, t_target_therm)
        while cosim.therm_integrator.t < t_target_therm
            SciMLBase.step!(cosim.therm_integrator)
        end
        
        t_bal += dt
    end
    print("\n>>> WORKSHOP BALANCING COMPLETE.\n")
end

function load_fcr_data(duration)
    t_data = Float64[]
    f_data = Float64[]
    try
        df = CSV.read(joinpath("data", "V2G", "RTE_Frequence_2024", "RTE_Frequence_2024_02.txt"), DataFrame; delim=';')
        freq_col = names(df)[occursin.(r"freq"i, names(df))][1]
        raw_f = df[!, freq_col]
        f_data = map(raw_f) do val
            if ismissing(val) return 50.0 
            elseif val isa Number return Float64(val)
            elseif val isa AbstractString return tryparse(Float64, replace(strip(val), "," => ".")) === nothing ? 50.0 : tryparse(Float64, replace(strip(val), "," => "."))
            else return 50.0 end
        end
        f_data = convert(Vector{Float64}, f_data)
        t_data = collect(0.0 : 10.0 : (length(f_data)-1)*10.0)
    catch
        t_data = collect(0.0 : 10.0 : duration)
        f_data = 50.0 .+ 0.15 .* sin.(t_data ./ 300.0) .+ 0.05 .* randn(length(t_data))
    end
    return t_data, f_data
end

function run_diagnostic_v2g()
    TARGET_DAYS = 90
    println(">>> INITIATING $TARGET_DAYS-DAY ADAPTIVE V2G DIAGNOSTIC...")
    
    # 1. Architecture Setup
    elec_params = Chen2020()
    elec_params.Vmin = 2.5; elec_params.Vmax = 4.2
    
    tms_geom = TMSGeometry(channel_width=0.002, channel_height=0.050, number_of_channels=7, wall_thickness=0.001)

    nema_params = PackParameters(
        fluid = get_coolant_properties(:water_nema2026), pipe_wall = get_solid_properties(:aluminium_nema2026),
        potting_material = get_solid_properties(:bergquist_tgf_1500), casing_material = get_solid_properties(:aluminium_nema2026),
        tms_geometry = tms_geom, cell_gap_thickness = 0.002, axial_potting_thickness = 0.005, casing_thickness = 0.01,
        ambient_temperature = 298.15, inlet_temperature = 298.15, mass_flow_rate = 0.01, ambient_convection_coefficient = 5.0 
    )
    tms = ReactiveTMS(T_high = 40.0, T_low = 37.0)
    
    geom = build_pack_geometry(rows=4, cols=7, cell_pitch=0.025)
    cosim = build_pack_simulator(geom, elec_params, nema_params, 4, 7; verbose=false, geom_sigma=0.005)
    
    println("\n>>> [PHASE 1] LOADING 5-YEAR CHECKPOINT AND BALANCING...")
    load_checkpoint!(cosim, CHECKPOINT_5_YEAR; silent=true)
    balance_pack!(cosim, 0.65, 48.0)
    
    # Save the starting point of Phase 2
    save_checkpoint(cosim, CHECKPOINT_V2G_MORNING; silent=true)
    
    # 2. V2G Setup (1-Day Schedule)
    pack_1C_amps = 5.0 * 7
    pack_v_max = 4.2 * 4
    t_fcr, f_fcr = load_fcr_data(13.5 * 3600.0)
    max_v2g_power_w = 2.0 * pack_1C_amps * (3.7 * 4)
    wltp_path = joinpath("data", "V2G", "driving_power_wltp.csv")
    if isfile(wltp_path)
        f_wltp = CSV.File(wltp_path) |> Tables.matrix
        target_max_power = 2.25 * pack_1C_amps * (3.7 * 4)
        wltp_scaler = target_max_power / maximum(abs.(f_wltp[:,2]))
        t_wltp = f_wltp[:,1]; p_wltp = -f_wltp[:,2] .* wltp_scaler; dt_wltp = diff(t_wltp)
        tend = findfirst(t_wltp .>= 1800.0)
        if isnothing(tend) tend = length(dt_wltp) end
        wltp_step = DriveStep(Any[dt_wltp[1:tend], p_wltp[1:tend]], 1800.0)
    else
        wltp_step = CurrentStep(-0.5 * pack_1C_amps, 1800.0)
    end

    v2g_daily_schedule = Step[
        wltp_step,                                                              
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 8.5 * 3600.0), 
        wltp_step,                                                              
        FCRStep(t_fcr, f_fcr, max_v2g_power_w, 13.5 * 3600.0),                  
        TargetCurrentStep(0.25 * pack_1C_amps, 0.61, pack_v_max, 1.0 * 3600.0)  
    ]
    
    exp_v2g_day = Experiment(v2g_daily_schedule)
    total_t_v2g_day = sum([s.period for s in v2g_daily_schedule])

    println("\n>>> [PHASE 2] EXECUTING $TARGET_DAYS-DAY ADAPTIVE LEARNING LOOP...")
    
    master_df = DataFrame()
    current_derate_exp = 1.0
    day = 1
    
    while day <= TARGET_DAYS
        retries = 0
        day_success = false
        
        while !day_success
            prefix = "[Day $(lpad(day, 2))/$TARGET_DAYS | Derate ^$(round(current_derate_exp, digits=1))] "
            
            try
                df_day = simulate_pack!(cosim, exp_v2g_day, nema_params, tms; 
                    total_time=total_t_v2g_day, save_csv=false, 
                    derate_exponent=current_derate_exp, print_prefix=prefix,
                    force_even_current=false, is_isothermal=false, geom_sigma=0.005,
                    m_active=0.01, m_passive=0.01, heat_multiplier=15.0,
                    opt_w_ratio=Dict(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
                    opt_alpha_elec=Dict(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
                    opt_dt_max_elec=Dict(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
                    opt_alpha_therm=Dict(:smooth => 0.1, :drive => 2.59, :grid => 5.9, :aging => 0.5, :rest => 0.1),
                    opt_dt_max_therm=Dict(:smooth => 158.33, :drive => 6.88, :grid => 57.92, :aging => 500.0, :rest => 3600.0),
                    dense_logging=false, verbose=true
                )
                
                # Check if any cell threw a fatal solver error (Unstable/Terminated)
                crashed_cell = 0
                for i in 1:cosim.num_cells
                    rc = cosim.cell_integrators[i].sol.retcode
                    if rc != SciMLBase.ReturnCode.Success && rc != SciMLBase.ReturnCode.Default
                        crashed_cell = i
                        break
                    end
                end
                
                if crashed_cell > 0
                    error("Cell $crashed_cell crashed.")
                end
                
                # Success!
                append!(master_df, df_day)
                save_checkpoint(cosim, CHECKPOINT_V2G_MORNING; silent=true)
                day_success = true
                day += 1
                
            catch e
                retries += 1
                print("\n[!] Crash detected on Day $day. Reloading morning checkpoint...")
                load_checkpoint!(cosim, CHECKPOINT_V2G_MORNING; silent=true)
                
                # If we crashed twice on the same day with this exponent, increment it!
                if retries >= 2
                    current_derate_exp += 0.1
                    print(" Failed twice. Increasing Derate Exponent to $(round(current_derate_exp, digits=1)).\n")
                    retries = 0 
                else
                    print(" Retrying...\n")
                end
            end
        end
    end
    
    println("\n>>> DIAGNOSTIC COMPLETE. SAVING MASTER CSV...")
    CSV.write("diagnostic_90_days_master.csv", master_df)
    println(">>> Saved to: diagnostic_90_days_master.csv")
end

Base.invokelatest(run_diagnostic_v2g)