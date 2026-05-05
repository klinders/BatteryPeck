# ==============================================================================
# anticipative_tms.jl
# Compares Baseline Reactive TMS vs. Anticipative Predictive TMS
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq
using Plots
using Revise
using Logging
using BatteryToolkit
using CSV
using DataFrames
using DataInterpolations
using DiffEqCallbacks

Revise.revise()

"""
Safely assign guesses to abstract arrays or scalar variables.
"""
function safe_guess!(dict, var, val)
    try
        arr = collect(var)
        if arr isa AbstractArray
            for el in arr dict[el] = val end
        else
            dict[var] = val
        end
    catch
        dict[var] = val
    end
end

"""
Execute full multi-cell diagnostic comparing reactive and anticipative controllers.
"""
function run_comparison_test()
    println("Initialising V2X coupled diagnostic...")

    elec_params = Chen2020()
    elec_params.Vmin = 1.0 
    elec_params.Vmax = 5.5
    
    lut_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
    lut = CSV.read(lut_file, DataFrame)
    v_to_soc_interp = LinearInterpolation(lut.SoC, lut.Voltage)
    soc_init = v_to_soc_interp(4.0) 
    
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊  
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊

    therm_params = build_pack_parameters()
    a_c = therm_params.tms_geometry.channel_width * therm_params.tms_geometry.channel_height * therm_params.tms_geometry.number_of_channels
    
    m_passive = velocity_to_mass_flow(0.1, therm_params.fluid.density, a_c)
    m_active = velocity_to_mass_flow(0.3, therm_params.fluid.density, a_c) 
    
    therm_params = PackParameters(
        fluid = therm_params.fluid, pipe_wall = therm_params.pipe_wall,
        potting_material = therm_params.potting_material, casing_material = therm_params.casing_material,
        tms_geometry = therm_params.tms_geometry, cell_gap_thickness = therm_params.cell_gap_thickness,
        axial_potting_thickness = therm_params.axial_potting_thickness, casing_thickness = therm_params.casing_thickness,
        ambient_temperature = 298.15, 
        inlet_temperature = 298.15,
        mass_flow_rate = m_passive, 
        ambient_convection_coefficient = 5.0 
    )

    geom = build_pack_geometry(rows=1, cols=4, pattern=:square, tms_routing=:single_row, tms_encasement=:start_bottom, cell_pitch=0.025)

    println("Building symbolic system...")
    sys = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        structural_simplify(CoupledMultiCellPack(name=:sys, geom=geom, therm_params=therm_params, elec_params=elec_params))
    end

    my_guesses = Dict()
    for i in 1:4
        c = getproperty(sys, i==1 ? :cell : Symbol("spme_$i"))
        if hasproperty(c, :sei)
            safe_guess!(my_guesses, c.sei.J.u, 1e-6)
            safe_guess!(my_guesses, c.sei.Δϕₛ.u, 0.1)
        end
    end

    exp = Experiment([
        RestStep(600.0),
        CurrentStep(2.5, 4000.0),
        RestStep(10.0),
        CurrentStep(-2.5, 4000.0),
        RestStep(10.0),
        CurrentStep(5.0, 2000.0),
        RestStep(10.0),
        CurrentStep(-5.0, 2000.0),
        RestStep(10.0),
        CurrentStep(10.0, 500.0),
        RestStep(10.0),
        CurrentStep(-10.0, 500.0),
        RestStep(10.0),
        CurrentStep(10.0, 500.0),
        RestStep(10.0),
        CurrentStep(-10.0, 500.0),
        RestStep(10.0),
        CurrentStep(2.5, 4000.0),
        RestStep(10.0),
        CurrentStep(-2.5, 4000.0),
        RestStep(10.0),
        CurrentStep(10.0, 500.0),
        RestStep(10.0),
        CurrentStep(-10.0, 500.0),
        RestStep(10.0),
        CurrentStep(10.0, 500.0),
        RestStep(600.0)
    ])

    total_time = exp.tend

    println("Compiling ODE problem...")
    prob = ODEProblem(sys, [sys.Pin => 0.0, sys.Iin => 0.0], (0.0, total_time); guesses=my_guesses, sparse=true, jac=false)

    println("Executing Run 1: Reactive Baseline...")
    cb_baseline = build_baseline_callback(sys.thermal_pack, 4, m_active, m_passive, 288.15, 298.15)
    
    saved_values_base = SavedValues(Float64, Tuple{Float64, Float64})
    cb_save_base = SavingCallback(
        (u, t, integrator) -> (
            integrator.ps[sys.thermal_pack.fluid_inlet.m_flow_in], 
            integrator.ps[sys.Iin]
        ), 
        saved_values_base;
        saveat = 1.0:1.0:total_time, 
        save_everystep = false       
    )
    cb_all_base = CallbackSet(cb_baseline, cb_save_base)

    sol_base = simulate(sys, prob, exp, QNDF(); saveat=1.0, dtmax=5.0, callback=cb_all_base)
    
    t_base = sol_base.t
    T_max_base = get_max_cell_temp(sol_base, sys, 4)
    V_pack_base = sol_base[sys.V]
    
    flow_base = [val[1] for val in saved_values_base.saveval]
    I_applied = [-val[2] for val in saved_values_base.saveval]
    t_scraped_base = saved_values_base.t

    println("Executing Run 2: Anticipative Controller...")
    cb_anticipative = build_anticipative_callback(sys.thermal_pack, 4, exp.steps, m_active, m_passive, 5.0, 300.0, 288.15, 298.15)
    
    saved_values_anti = SavedValues(Float64, Tuple{Float64, Float64})
    cb_save_anti = SavingCallback(
        (u, t, integrator) -> (
            integrator.ps[sys.thermal_pack.fluid_inlet.m_flow_in], 
            integrator.ps[sys.Iin]
        ), 
        saved_values_anti;
        saveat = 1.0:1.0:total_time, 
        save_everystep = false       
    )
    cb_all_anti = CallbackSet(cb_anticipative, cb_save_anti)

    sol_anti = simulate(sys, prob, exp, QNDF(); saveat=1.0, dtmax=5.0, callback=cb_all_anti)

    t_anti = sol_anti.t
    T_max_anti = get_max_cell_temp(sol_anti, sys, 4)
    V_pack_anti = sol_anti[sys.V]
    
    flow_anti = [val[1] for val in saved_values_anti.saveval]
    t_scraped_anti = saved_values_anti.t

    println("Rendering results...")
    
    # Create the 4 distinct plot panels stacked vertically
    p1 = plot(t_scraped_base, I_applied, label="Applied Current", ylabel="[A]", title="Dynamic Load Profile", lw=2, color=:black, fill=(0, 0.2, :black), linetype=:steppost)
    
    p2 = plot(t_base, V_pack_base, label="Reactive Voltage", ylabel="[V]", title="Pack Voltage", lw=2, color=:red, linestyle=:dash)
    plot!(p2, t_anti, V_pack_anti, label="Anticipative Voltage", lw=2, color=:blue, alpha=0.7)
    
    p3 = plot(t_scraped_base, flow_base, label="Reactive Pump", ylabel="[kg/s]", title="Coolant Mass Flow", lw=2, color=:red, linestyle=:dash, linetype=:steppost)
    plot!(p3, t_scraped_anti, flow_anti, label="Anticipative Pump", lw=2, color=:blue, alpha=0.7, linetype=:steppost)
    
    # Added xlabel here for the bottom-most plot
    p4 = plot(t_base, T_max_base, label="Reactive Peak Temp", ylabel="[C]", xlabel="Time [s]", title="Maximum Cell Core Temperature", lw=2, color=:red, linestyle=:dash)
    plot!(p4, t_anti, T_max_anti, label="Anticipative Peak Temp", lw=2, color=:blue)
    hline!(p4, [35.0], label="35C Safety Limit", color=:gray, linestyle=:dot, lw=2)

    # Display using a 4x1 vertical stack layout. Height increased to 1200 for breathing room.
    display(plot(p1, p2, p3, p4, layout=(4,1), size=(900,1200)))
    println("Dynamic Comparison Complete!")
end

"""
Extract maximum cell core temperature across all pack elements.
"""
function get_max_cell_temp(sol, sys, num_cells)
    temps = hcat([sol[getproperty(sys.thermal_pack, Symbol("cell_$i")).core_cap.T] .- 273.15 for i in 1:num_cells]...)
    return maximum(temps, dims=2)[:, 1]
end

Base.invokelatest(run_comparison_test)