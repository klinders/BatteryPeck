# ==============================================================================
# validation_nema_11.jl     (hydraulics)
# Validates the LPTN pressure drop and pump power against Nema 2026 (Fig 11a/11b)
#
# FORENSIC DISCOVERIES & CORRECTIONS VS. NEMA 2026:
#
# 1. Pump Power Typo: Nema reports pump power in mW (e.g., 19.6 mW at 0.5 m/s), 
#    but raw hydraulic physics (W = Q * dP) based on their own reported pressure 
#    drops yield a power exactly 10x higher (~196 mW). The paper contains a 
#    post-processing decimal/unit conversion error.
#
# 2. Viscosity Mismatch: Nema Table 1 claims a dynamic viscosity of 0.002 Pa.s. 
#    However, at low velocities (0.1 m/s), matching their 0.19 kPa pressure drop 
#    mathematically requires using standard water properties (~0.00089 Pa.s).
#
# 3. 1D Minor Losses: To match the 3D CFD's quadratic (v^2) pressure drop curve, 
#    a minor loss coefficient (K ≈ 35.0 per meter) must be added to the 1D LPTN 
#    to account for form drag caused by the 180-degree serpentine bends.
#
# 4. Fluid Gap Geometry: To prevent linear viscous friction from severely 
#    overestimating pressure drop, the hydraulic diameter must be calculated 
#    using the actual cross-section of the fluid (not the bounding box).
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

Loads CSV data and extracts time and value columns.
"""
function load_and_preprocess(filepath)
    raw_data = readdlm(filepath, ',')
    if typeof(raw_data[1,1]) <: AbstractString
        raw_data = raw_data[2:end, :]
    end
    raw_data = Float64.(raw_data)
    
    return raw_data[:, 1], raw_data[:, 2]
end

data_dir = joinpath(@__DIR__, "..", "data", "Nema2026")

v_val_dp_40, dp_val_40 = load_and_preprocess(joinpath(data_dir, "11a_PressureDrop.csv"))
v_val_W_40, W_val_40   = load_and_preprocess(joinpath(data_dir, "11a_PumpPower.csv"))

v_val_dp_50, dp_val_50 = load_and_preprocess(joinpath(data_dir, "11b_PressureDrop.csv"))
v_val_W_50, W_val_50   = load_and_preprocess(joinpath(data_dir, "11b_PumpPower.csv"))

@parameters t

"""
    build_hydraulic_system(name::Symbol, geom, params::PackParameters)

Constructs isothermal 1D fluid network to evaluate pressure drop across serpentine channels.
Includes minor loss components to account for form drag in bends.

Editable values:
`K_per_metre`: Adjusts form drag multiplier for serpentine bends. Directly affects quadratic pressure drop curve.
"""
function build_hydraulic_system(name::Symbol, geom, params::PackParameters)
    num_tms_nodes = length(geom.tms_coords)
    
    total_length = 0.0
    for edge in geom.flow_edges
        t1, t2 = edge
        c1 = geom.tms_coords[t1]
        c2 = geom.tms_coords[t2]
        total_length += sqrt((c1[1] - c2[1])^2 + (c1[2] - c2[2])^2)
    end
    node_length = total_length / num_tms_nodes
    
    tms_nodes = [TMSNode(name=Symbol("tms_$i"), params=params, length=node_length) for i in 1:num_tms_nodes]
    
    K_per_metre = 33.0
    minor_losses = []
    
    for (idx, edge) in enumerate(geom.flow_edges)
        t1, t2 = edge
        c1 = geom.tms_coords[t1]
        c2 = geom.tms_coords[t2]
        dist = sqrt((c1[1] - c2[1])^2 + (c1[2] - c2[2])^2)
        
        node_K = K_per_metre * dist
        push!(minor_losses, MinorLoss(name=Symbol("bend_$idx"), params=params, K=node_K))
    end
    
    eqs = Equation[]
    
    for (idx, edge) in enumerate(geom.flow_edges)
        t1, t2 = edge
        push!(eqs, connect(tms_nodes[t1].port_b, minor_losses[idx].port_a))
        push!(eqs, connect(minor_losses[idx].port_b, tms_nodes[t2].port_a))
    end
    
    first_tms = tms_nodes[geom.flow_edges[1][1]]
    last_tms  = tms_nodes[geom.flow_edges[end][2]]
    
    @named fluid_inlet = FluidSource(m_flow_val = params.mass_flow_rate, T_val = params.inlet_temperature)
    @named fluid_outlet = FluidSink(p_val = 101325.0)
    
    push!(eqs, connect(fluid_inlet.port, first_tms.port_a))
    push!(eqs, connect(last_tms.port_b, fluid_outlet.port))
    
    all_systems = vcat(tms_nodes, minor_losses, [fluid_inlet, fluid_outlet])
    
    return ODESystem(eqs, t, [], []; systems=all_systems, name=name)
end

"""
    run_hydraulics(channel_width, velocities)

Simulates pressure drop and calculates pump power for various inlet velocities.

Equations:
Mass flow rate: m_dot = ρ * A_cross * v
Pump power: W = Q * ΔP = (m_dot / ρ) * ΔP

Editable values:
`velocities`: Array defining fluid simulation speeds.
"""
function run_hydraulics(channel_width, velocities)
    
    channel_height = 0.004 
    A_cross = channel_width * channel_height
    rho_water = 998.2 
    
    tms_geom = TMSGeometry(
        channel_width = channel_width, 
        channel_height = channel_height, 
        number_of_channels = 1, 
        wall_thickness = 0.001
    )
    
    geom = build_pack_geometry(
        rows = 4, cols = 7, pattern = :hexagonal, 
        tms_routing = :single_row, tms_encasement = :full, cell_pitch = 0.025
    )
    
    dp_results = Float64[]
    power_results = Float64[]
    power_results_corrected = Float64[]
    
    for v in velocities
        
        m_dot = rho_water * A_cross * v
        
        params = PackParameters(
            fluid = get_coolant_properties(:water_nema2026),
            pipe_wall = get_solid_properties(:aluminium_nema2026),
            potting_material = get_solid_properties(:bergquist_tgf_1500),
            casing_material = get_solid_properties(:aluminium_nema2026),
            tms_geometry = tms_geom,
            
            cell_gap_thickness = (0.025 - 0.021) / 2.0,
            axial_potting_thickness = 0.005,
            casing_thickness = 0.003,
            
            ambient_temperature = 298.15,
            inlet_temperature = 298.15,
            mass_flow_rate = m_dot,
            ambient_convection_coefficient = 5.0
        )
        
        sys_simplified = with_logger(ConsoleLogger(stderr, Logging.Error)) do
            raw_sys = build_hydraulic_system(:hyd_sys, geom, params)
            structural_simplify(raw_sys)
        end
        
        prob = ODEProblem(sys_simplified, [], (0.0, 1.0), jac=true, sparse=true)
        sol = solve(prob, QNDF(), saveat=1.0)
        
        p_in  = sol[sys_simplified.fluid_inlet.port.p][end]
        p_out = sol[sys_simplified.fluid_outlet.port.p][end]
        dp_Pa = p_in - p_out
        
        dp_kPa = dp_Pa / 1000.0
        push!(dp_results, dp_kPa)
        
        Q_m3_s = m_dot / rho_water
        W_watts = Q_m3_s * dp_Pa
        W_mW = W_watts * 1000.0
        
        push!(power_results_corrected, W_mW)
        
        push!(power_results, W_mW / 10.0) 
    end
    
    return dp_results, power_results, power_results_corrected
end

velocities = 0.1:0.05:0.5

results_dp = Dict()
results_W = Dict()
results_W_corrected = Dict()

for w in [0.040, 0.050]
    println("Running hydraulics for width = $(w*1000) mm...")
    dp, w_p, w_c = run_hydraulics(w, velocities)
    results_dp[w] = dp
    results_W[w] = w_p
    results_W_corrected[w] = w_c
end

p1 = plot(
    velocities, results_dp[0.040], label="LPTN (40mm)", 
    color=:blue, linewidth=2, marker=:circle,
    xlabel="Velocity (m/s)", ylabel="Pressure Drop (kPa)",
    title="Fig 11a: DP (Corrected CFD)", legend=:topleft, grid=true
)
scatter!(p1, v_val_dp_40, dp_val_40, label="Nema CFD", color=:red, markersize=5)

p2 = plot(
    velocities, results_W[0.040], label="LPTN (40mm)", 
    color=:green, linewidth=2, marker=:circle,
    xlabel="Velocity (m/s)", ylabel="Power (mW)",
    title="Fig 11a: Pump Power", legend=:topleft, grid=true
)
scatter!(p2, v_val_W_40, W_val_40, label="Nema CSV (Typo)", color=:orange, markersize=4)
scatter!(p2, velocities, results_W_corrected[0.040], label="CSV Corrected", color=:purple, markershape=:star5, markersize=6)

p3 = plot(
    velocities, results_dp[0.050], label="LPTN (50mm)", 
    color=:blue, linewidth=2, marker=:circle,
    xlabel="Velocity (m/s)", ylabel="Pressure Drop (kPa)",
    title="Fig 11b: DP (Flawed CFD)", legend=:topleft, grid=true
)
scatter!(p3, v_val_dp_50, dp_val_50, label="Nema CFD", color=:red, markersize=5)

p4 = plot(
    velocities, results_W[0.050], label="LPTN (50mm)", 
    color=:green, linewidth=2, marker=:circle,
    xlabel="Velocity (m/s)", ylabel="Power (mW)",
    title="Fig 11b: Pump Power", legend=:topleft, grid=true
)
scatter!(p4, v_val_W_50, W_val_50, label="Nema CSV (Typo)", color=:orange, markersize=4)
scatter!(p4, velocities, results_W_corrected[0.050], label="CSV Corrected", color=:purple, markershape=:star5, markersize=6)

display(plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800)))