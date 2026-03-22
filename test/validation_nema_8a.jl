# ==============================================================================
# validation_nema_8a.jl     (passive cooling)
# Validates passive LPTN thermal network against Nema et al. (2026) Fig 8a
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq 
using DataInterpolations
using DelimitedFiles
using Plots
using Logging

using ModelingToolkitStandardLibrary.Thermal
using Revise
using BatteryToolkit

Revise.revise()

"""
    load_and_preprocess(filepath)

Loads CSV data and filters out backward-pass digitisation artifacts.
Ensures strictly increasing time array for interpolation stability.
"""
function load_and_preprocess(filepath)
    raw_data = readdlm(filepath, ',')
    if typeof(raw_data[1,1]) <: AbstractString
        raw_data = raw_data[2:end, :]
    end
    raw_data = Float64.(raw_data)
    
    t_raw = raw_data[:, 1]
    y_raw = raw_data[:, 2]
    
    max_idx = argmax(t_raw)
    t_clean = t_raw[1:max_idx]
    y_clean = y_raw[1:max_idx]
    
    valid_idx = [1]
    for i in 2:length(t_clean)
        if t_clean[i] > t_clean[valid_idx[end]]
            push!(valid_idx, i)
        end
    end
    
    t_final = t_clean[valid_idx]
    y_final = y_clean[valid_idx]
    
    t_shifted = t_final .- minimum(t_final)
    return t_shifted, y_final
end

"""
    compute_metrics(t_sim, y_sim, t_val, y_val)

Calculates RMSE (Root Mean Square Error) and R² (Coefficient of Determination).
"""
function compute_metrics(t_sim, y_sim, t_val, y_val)
    sim_interp = LinearInterpolation(y_sim, t_sim)
    valid_idx = findall(t -> t <= t_sim[end], t_val)
    t_val_trimmed = t_val[valid_idx]
    y_val_trimmed = y_val[valid_idx]
    
    y_sim_at_val = [sim_interp(t) for t in t_val_trimmed]
    rmse = sqrt(sum((y_sim_at_val .- y_val_trimmed).^2) / length(y_val_trimmed))
    y_mean = sum(y_val_trimmed) / length(y_val_trimmed)
    ss_tot = sum((y_val_trimmed .- y_mean).^2)
    ss_res = sum((y_val_trimmed .- y_sim_at_val).^2)
    r2 = 1.0 - (ss_res / ss_tot)
    
    return rmse, r2
end

data_dir = joinpath(@__DIR__, "..", "data", "Nema2026")

t_hg_1c, y_hg_1c = load_and_preprocess(joinpath(data_dir, "7a_298K_HG.csv"))
const hg_interp_1c = LinearInterpolation(y_hg_1c, t_hg_1c)

function q_gen_1c(t)
    if t < t_hg_1c[1] return y_hg_1c[1]
    elseif t > t_hg_1c[end] return y_hg_1c[end]
    else return hg_interp_1c(t) end
end
@register_symbolic q_gen_1c(t)

t_hg_5c, y_hg_5c = load_and_preprocess(joinpath(data_dir, "7b_298K_HG.csv"))
const hg_interp_5c = LinearInterpolation(y_hg_5c, t_hg_5c)

function q_gen_5c(t)
    if t < t_hg_5c[1] return y_hg_5c[1]
    elseif t > t_hg_5c[end] return y_hg_5c[end]
    else return hg_interp_5c(t) end
end
@register_symbolic q_gen_5c(t)

val_params = build_pack_parameters()
geom = build_pack_geometry(
    rows = 4, cols = 7, pattern = :hexagonal, 
    tms_routing = :single_row, tms_encasement = :full, cell_pitch = 0.025
)
num_cells = length(geom.cell_coords)

@parameters t

"""
    build_passive_system(name::Symbol, geom, params::PackParameters)

Constructs passive natural convection system bypassing all liquid cooling nodes.
Prevents zero-flow singularities by isolating active fluid physics.
"""
function build_passive_system(name::Symbol, geom, params::PackParameters)
    cells = [CoreShellCell(name=Symbol("cell_$i"), T_start=params.ambient_temperature) 
             for i in 1:num_cells]
             
    @named ambient_temp = FixedTemperature(T = params.ambient_temperature)
    
    casing_area = num_cells * (0.021 * 0.021) 
    casing_vol = casing_area * params.casing_thickness
    C_casing = casing_vol * params.casing_material.density * params.casing_material.specific_heat
    @named casing_mass = HeatCapacitor(C = C_casing, T = params.ambient_temperature)
    
    R_conv_val = 1.0 / (params.ambient_convection_coefficient * casing_area)
    @named casing_convection = ThermalResistor(R = R_conv_val)
    
    eqs = Equation[]
    
    push!(eqs, connect(casing_mass.port, casing_convection.port_a))
    push!(eqs, connect(casing_convection.port_b, ambient_temp.port))
    
    cell_face_area = pi * (params.tms_geometry.channel_width / 2.0)^2
    R_axial_val = params.axial_potting_thickness / (params.potting_material.thermal_conductivity * (2 * cell_face_area))
    axial_resistors = [ThermalResistor(name=Symbol("R_ax_$i"), R=R_axial_val) for i in 1:num_cells]
    
    for i in 1:num_cells
        push!(eqs, connect(cells[i].port_shell, axial_resistors[i].port_a))
        push!(eqs, connect(axial_resistors[i].port_b, casing_mass.port))
    end
    
    R_radial_val = params.cell_gap_thickness / (params.potting_material.thermal_conductivity * (0.065 * 0.021)) 
    gap_resistors = [ThermalResistor(name=Symbol("R_gap_$idx"), R=R_radial_val) for idx in 1:length(geom.cell_edges)]
    
    for (idx, edge) in enumerate(geom.cell_edges)
        c1, c2 = edge
        push!(eqs, connect(cells[c1].port_shell, gap_resistors[idx].port_a))
        push!(eqs, connect(gap_resistors[idx].port_b, cells[c2].port_shell))
    end
    
    systems = Any[ambient_temp, casing_mass, casing_convection]
    append!(systems, cells)
    append!(systems, axial_resistors)
    append!(systems, gap_resistors)

    return ODESystem(eqs, t, [], []; systems=systems, name=name)
end

"""
    run_passive_validation(rate_name, q_func, t_hg_end, val_filename)

Builds passive natural convection system, runs solver, calculates metrics, and plots results.
"""
function run_passive_validation(rate_name, q_func, t_hg_end, val_filename)
    println("\nStarting validation for: $rate_name discharge (Figure 8a)")
    
    sys_simplified, cells = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        sys_base = build_passive_system(Symbol("sys_$rate_name"), geom, val_params)
        cls = [getproperty(sys_base, Symbol("cell_$i")) for i in 1:num_cells]
        forcing_eqs = [cls[i].Q_volumetric_in.u ~ q_func(t) for i in 1:num_cells]
        full_sys = ODESystem(forcing_eqs, t; systems=[sys_base], name=Symbol("full_$rate_name"))
        (structural_simplify(full_sys), cls) 
    end

    prob = ODEProblem(sys_simplified, [], (0.0, t_hg_end), jac=true, sparse=true)
    @time sol = solve(prob, QNDF(), saveat=1.0, dtmax=1.0) 

    t_sim = sol.t
    
    T_cells_K = [sol[cells[i].core_cap.T] for i in 1:num_cells]
    T_max_sim_K = maximum(hcat(T_cells_K...), dims=2)[:, 1]
    T_max_sim_C = T_max_sim_K .- 273.15

    t_val, y_val_K = load_and_preprocess(joinpath(data_dir, val_filename))
    y_val_C = y_val_K .- 273.15
    
    rmse, r2 = compute_metrics(t_sim, T_max_sim_C, t_val, y_val_C)
    
    println("  -> RMSE: $(round(rmse, digits=3)) °C")
    println("  -> R²:   $(round(r2, digits=4))")

    plt = plot(
        t_sim, T_max_sim_C, 
        label="LPTN Passive", 
        linewidth=2, color=:green,
        xlabel="Time (s)", ylabel="Temperature (°C)",
        title="Nema et al. (2026) - Fig 8a ($rate_name No Cooling)",
        legend=:bottomright, grid=true
    )
    scatter!(
        plt, t_val, y_val_C, 
        label="Nema 3D CFD Data", 
        markershape=:circle, color=:red, markersize=4, alpha=0.7
    )
    annotate!(plt, [(t_hg_end*0.1, maximum(y_val_C)*0.95, text("RMSE: $(round(rmse, digits=2)) °C\nR²: $(round(r2, digits=3))", 10, :left))])
    
    return plt
end

p_1c = run_passive_validation("1C", q_gen_1c, t_hg_1c[end], "8a_1C.csv")
p_5c = run_passive_validation("5C", q_gen_5c, t_hg_5c[end], "8a_5C.csv")

display(plot(p_1c, p_5c, layout=(1, 2), size=(1000, 400)))