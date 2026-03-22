# ==============================================================================
# validation_nema.jl
# Validates 1D/2D LPTN thermal network against Nema et al. (2026) 3D CFD data
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq 
using DataInterpolations
using DelimitedFiles
using Plots
using Logging

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

Equations:
RMSE = sqrt( sum((y_sim - y_val)^2) / N )
R² = 1 - ( sum((y_val - y_sim)^2) / sum((y_val - y_mean)^2) )
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

"""
    q_gen_1c(t)

Returns volumetric heat generation for 1C discharge rate.
Interpolates Nema et al. (2026) Figure 7a data.
"""
function q_gen_1c(t)
    if t < t_hg_1c[1] return y_hg_1c[1]
    elseif t > t_hg_1c[end] return y_hg_1c[end]
    else return hg_interp_1c(t) end
end
@register_symbolic q_gen_1c(t)

t_hg_5c, y_hg_5c = load_and_preprocess(joinpath(data_dir, "7b_298K_HG.csv"))
const hg_interp_5c = LinearInterpolation(y_hg_5c, t_hg_5c)

"""
    q_gen_5c(t)

Returns volumetric heat generation for 5C discharge rate.
Interpolates Nema et al. (2026) Figure 7b data.
"""
function q_gen_5c(t)
    if t < t_hg_5c[1] return y_hg_5c[1]
    elseif t > t_hg_5c[end] return y_hg_5c[end]
    else return hg_interp_5c(t) end
end
@register_symbolic q_gen_5c(t)

# Sets physical environment and baseline boundaries
# Editable values:
# `v_coolant`: Alters fluid velocity, changing convective heat transfer and pressure drop
# `W_channel`: Alters width of cooling channel, affecting total mass flow and cooling area
# `H_channel`: Alters height of cooling channel
v_coolant = 0.3 
rho_water = 998.0
W_channel = 0.050 
H_channel = 0.004 
A_cross = W_channel * H_channel
m_dot_baseline = rho_water * A_cross * v_coolant 

tms_geom_nema = TMSGeometry(
    channel_width = W_channel, 
    channel_height = H_channel, 
    number_of_channels = 1, 
    wall_thickness = 0.001  
)

val_params = PackParameters(
    fluid = get_coolant_properties(:water_nema2026),
    pipe_wall = get_solid_properties(:aluminium_nema2026),
    potting_material = get_solid_properties(:bergquist_tgf_1500),
    casing_material = get_solid_properties(:aluminium_nema2026),
    tms_geometry = tms_geom_nema,
    
    cell_gap_thickness = (0.025 - 0.021) / 2.0, 
    axial_potting_thickness = 0.005,
    casing_thickness = 0.003,
    
    ambient_temperature = 298.15,
    inlet_temperature = 298.15,
    mass_flow_rate = m_dot_baseline, 
    ambient_convection_coefficient = 5.0
)

geom = build_pack_geometry(
    rows = 4, cols = 7, pattern = :hexagonal, 
    tms_routing = :single_row, tms_encasement = :full, cell_pitch = 0.025
)
num_cells = length(geom.cell_coords)

@parameters t

"""
    run_and_validate(rate_name, q_func, t_hg_end, val_filename)

Builds LPTN system, runs solver, calculates metrics, and plots validation results.
Extracts absolute peak centre temperature (T_max) from volume-averaged core node (T_avg).

Equations:
T_max = T_avg + (Q_watts * R_offset)
where R_offset (0.535 K/W) is difference between centre-to-surface and average-to-surface radial resistance.
"""
function run_and_validate(rate_name, q_func, t_hg_end, val_filename)
    println("\nStarting validation for: $rate_name discharge")
    
    sys_simplified, cells = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        sys_base = build_pack_system(Symbol("sys_$rate_name"), geom, val_params)
        cls = [getproperty(sys_base, Symbol("cell_$i")) for i in 1:num_cells]
        forcing_eqs = [cls[i].Q_volumetric_in.u ~ q_func(t) for i in 1:num_cells]
        full_sys = ODESystem(forcing_eqs, t; systems=[sys_base], name=Symbol("full_$rate_name"))
        (structural_simplify(full_sys), cls) 
    end

    prob = ODEProblem(sys_simplified, [], (0.0, t_hg_end), jac=true, sparse=true)
    @time sol = solve(prob, QNDF(), saveat=1.0, dtmax=1.0) 

    t_sim = sol.t
    
    T_avg_cells_K = [sol[cells[i].core_cap.T] for i in 1:num_cells]
    T_avg_max_K = maximum(hcat(T_avg_cells_K...), dims=2)[:, 1]
    T_avg_max_C = T_avg_max_K .- 273.15

    base_interp = rate_name == "1C" ? hg_interp_1c : hg_interp_5c
    V_jellyroll = 2.13e-5
    
    T_max_sim_C = Float64[]
    for i in 1:length(t_sim)
        current_t = t_sim[i]
        
        q_volumetric = current_t <= base_interp.t[end] ? base_interp(current_t) : base_interp.u[end]
        q_watts = q_volumetric * V_jellyroll 
        
        t_max_instant = T_avg_max_C[i] + (q_watts * 0.535)
        push!(T_max_sim_C, t_max_instant)
    end

    t_val, y_val_K = load_and_preprocess(joinpath(data_dir, val_filename))
    y_val_C = y_val_K .- 273.15
    
    rmse, r2 = compute_metrics(t_sim, T_max_sim_C, t_val, y_val_C)
    
    println("  -> RMSE: $(round(rmse, digits=3)) °C")
    println("  -> R²:   $(round(r2, digits=4))")

    plt = plot(
        t_sim, T_max_sim_C, 
        label="LPTN T_max (Centre)", 
        linewidth=2, color=:blue,
        xlabel="Time (s)", ylabel="Temperature (°C)",
        title="Nema et al. (2026) - $rate_name validation",
        legend=:bottomright, grid=true
    )
    plot!(
        plt, t_sim, T_avg_max_C,
        label="LPTN T_avg (Volume)",
        linewidth=2, color=:cyan, linestyle=:dash
    )
    scatter!(
        plt, t_val, y_val_C, 
        label="Nema 3D CFD Data", 
        markershape=:circle, color=:red, markersize=4, alpha=0.7
    )
    annotate!(plt, [(t_hg_end*0.1, maximum(y_val_C)*0.95, text("RMSE: $(round(rmse, digits=2)) °C\nR²: $(round(r2, digits=3))", 10, :left))])
    
    return plt
end

p_1c = run_and_validate("1C", q_gen_1c, t_hg_1c[end], "8b_1C.csv")
p_5c = run_and_validate("5C", q_gen_5c, t_hg_5c[end], "8b_5C.csv")

display(plot(p_1c, p_5c, layout=(1, 2), size=(1000, 400)))