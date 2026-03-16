# ==============================================================================
# test_geometry.jl
# Render generated pack geometry and verify multi physics graph connections
# ==============================================================================

using Plots
using Plots.Measures

using Revise
using BatteryToolkit

Revise.revise()

# Draw circles in physical data coordinates
function circle_shape(x, y, r)
    θ = LinRange(0, 2π, 72)
    return x .+ r .* cos.(θ), y .+ r .* sin.(θ)
end

let
    n_rows = 5
    n_cols = 6

    println("Generating geometry map...")
    
    geom = Base.invokelatest(() -> build_pack_geometry(
        rows=n_rows, cols=n_cols, 
        pattern=:square, 
        tms_routing=:double_row, 
        tms_encasement=:start_bottom,
        cell_pitch=0.023
    ))

    println("Calculating unified plot boundaries...")
    
    # Extract coordinates for bounding box
    all_x = [c[1] for c in geom.cell_coords] ∪ [t[1] for t in geom.tms_coords]
    all_y = [c[2] for c in geom.cell_coords] ∪ [t[2] for t in geom.tms_coords]
    
    min_x, max_x = minimum(all_x), maximum(all_x)
    min_y, max_y = minimum(all_y), maximum(all_y)
    
    cell_radius = 0.021 / 2.0 
    
    # Calculate padding to prevent clipping
    pad_x = (max_x - min_x) * 0.05 + cell_radius * 1.5
    pad_y = (max_y - min_y) * 0.05 + cell_radius * 1.5
    
    shared_xlims = (min_x - pad_x, max_x + pad_x)
    shared_ylims = (min_y - pad_y, max_y + pad_y)

    println("Rendering geometry subplots...")

    # Initialise thermal subplot
    p_therm = plot(title="Thermal graph (cell conduction, convection, and fluid flow)", 
                   aspect_ratio=:equal, legend=:outerright, margin=0mm, 
                   framestyle=:none, xlims=shared_xlims, ylims=shared_ylims)

    # Draw battery cells
    for (idx, c) in enumerate(geom.cell_coords)
        lbl = idx == 1 ? "Battery cell" : false
        plot!(p_therm, circle_shape(c[1], c[2], cell_radius), seriestype=:shape, 
              color=:lightgray, linecolor=:gray, label=lbl)
    end

    # Draw cell conduction edges
    for (idx, (c_in, c_out, dist)) in enumerate(geom.cell_edges)
        x_coords = [geom.cell_coords[c_in][1], geom.cell_coords[c_out][1]]
        y_coords = [geom.cell_coords[c_in][2], geom.cell_coords[c_out][2]]
        lbl = idx == 1 ? "Cell conduction" : false
        plot!(p_therm, x_coords, y_coords, color=:gray, alpha=0.5, linewidth=2, label=lbl)
    end

    # Draw TMS convection edges
    for (idx, (c_idx, t_idx, dist)) in enumerate(geom.convection_edges)
        x_coords = [geom.cell_coords[c_idx][1], geom.tms_coords[t_idx][1]]
        y_coords = [geom.cell_coords[c_idx][2], geom.tms_coords[t_idx][2]]
        lbl = idx == 1 ? "TMS convection" : false
        plot!(p_therm, x_coords, y_coords, color=:magenta, alpha=0.6, linewidth=3, label=lbl)
    end

    # Draw fluid flow path
    for (idx, (t_in, t_out, dist)) in enumerate(geom.flow_edges)
        x_coords = [geom.tms_coords[t_in][1], geom.tms_coords[t_out][1]]
        y_coords = [geom.tms_coords[t_in][2], geom.tms_coords[t_out][2]]
        lbl = idx == 1 ? "Fluid flow path" : false
        plot!(p_therm, x_coords, y_coords, color=:blue, linewidth=3, label=lbl)
    end

    # Draw TMS nodes
    tms_x = [t[1] for t in geom.tms_coords]
    tms_y = [t[2] for t in geom.tms_coords]
    scatter!(p_therm, tms_x, tms_y, color=:white, markersize=4, markerstrokecolor=:blue, markerstrokewidth=2, label="TMS node")

    # Initialise electrical subplot
    p_elec = plot(title="Electrical graph (series and parallel busbars)", 
                  aspect_ratio=:equal, legend=:outerright, margin=0mm,
                  framestyle=:none, xlims=shared_xlims, ylims=shared_ylims)

    # Draw battery cells
    for (idx, c) in enumerate(geom.cell_coords)
        lbl = idx == 1 ? "Battery cell" : false
        plot!(p_elec, circle_shape(c[1], c[2], cell_radius), seriestype=:shape, 
              color=:lightgray, linecolor=:gray, label=lbl)
    end

    # Draw busbar connections
    labeled_parallel, labeled_series = false, false
    for (c_in, c_out, conn_type) in geom.electrical_edges
        x_coords = [geom.cell_coords[c_in][1], geom.cell_coords[c_out][1]]
        y_coords = [geom.cell_coords[c_in][2], geom.cell_coords[c_out][2]]
        if conn_type == :parallel
            lbl = !labeled_parallel ? "Parallel busbar" : false
            plot!(p_elec, x_coords, y_coords, color=:orange, alpha=0.8, linewidth=4, label=lbl)
            labeled_parallel = true
        elseif conn_type == :series
            lbl = !labeled_series ? "Series busbar" : false
            plot!(p_elec, x_coords, y_coords, color=:purple, alpha=0.6, linewidth=3, label=lbl)
            labeled_series = true
        end
    end

    # Initialise ambient boundary subplot
    p_ambient = plot(title="Ambient boundary (lateral exposed perimeter fraction)", 
                     aspect_ratio=:equal, margin=0mm, legend=false,
                     framestyle=:none, xlims=shared_xlims, ylims=shared_ylims)

    # Add invisible scatter to force colour bar rendering
    scatter!(p_ambient, [min_x], [min_y], zcolor=[0], clims=(0, 1), color=:inferno, 
             label="", markersize=0, markerstrokewidth=0, colorbar=true)

    # Draw heat map cells
    for (idx, (c_idx, exp_frac)) in enumerate(geom.casing_edges)
        c = geom.cell_coords[c_idx]
        col = cgrad(:inferno)[exp_frac] 
        plot!(p_ambient, circle_shape(c[1], c[2], cell_radius), seriestype=:shape, 
              color=col, linecolor=:black, linewidth=1, label="")
    end

    # Combine subplots into dashboard
    p_dashboard = plot(p_therm, p_elec, p_ambient, layout=(3, 1), size=(900, 1000))
    display(p_dashboard)
    
    # Find bounds dynamically
    idx_external = 1 
    idx_internal = clamp(length(geom.cell_coords) ÷ 2 + 1, 1, length(geom.cell_coords))
    
    # Translate linear index to row/col coordinates
    r_ext = (idx_external - 1) ÷ n_cols + 1
    c_ext = (idx_external - 1) % n_cols + 1
    
    r_int = (idx_internal - 1) ÷ n_cols + 1
    c_int = (idx_internal - 1) % n_cols + 1
    
    val_ext = geom.casing_edges[idx_external][2]
    val_int = geom.casing_edges[idx_internal][2]

    println("\n[Verification Report]")
    println("- Lateral exposure (external corner cell at row $r_ext, col $c_ext): $(round(val_ext * 100, digits=2))%")
    println("- Lateral exposure (internal core cell at row $r_int, col $c_int): $(round(val_int * 100, digits=2))%")
end