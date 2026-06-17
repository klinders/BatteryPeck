# ==============================================================================
# PackGeometry.jl
# Define battery pack spatial layout and auto generate multi-physics graph edges
# Contains:
# 1. PackGeometry: Store spatial coordinates and connecting edges for thermal, fluid, and electrical graphs
# 2. build_pack_geometry: Generate spatial layout and multi-physics graph edges for battery pack
# ==============================================================================

export PackGeometry, build_pack_geometry

"""
    PackGeometry

Store spatial coordinates and connecting edges for thermal, fluid, and electrical graphs.
"""
struct PackGeometry
    # Define physical node locations
    cell_coords::Vector{Tuple{Float64, Float64}}
    tms_coords::Vector{Tuple{Float64, Float64}}
    
    # Define thermal network routing
    cell_edges::Vector{Tuple{Int, Int, Float64}}
    convection_edges::Vector{Tuple{Int, Int, Float64}}
    casing_edges::Vector{Tuple{Int, Float64}} 
    
    # Define fluid and electrical network routing
    flow_edges::Vector{Tuple{Int, Int, Float64}}
    electrical_edges::Vector{Tuple{Int, Int, Symbol}}
end

"""
    build_pack_geometry(; rows, cols, pattern, tms_routing, tms_encasement, cell_pitch)

Generate spatial layout and multi-physics graph edges for battery pack.

# Arguments
- `rows::Int`: Quantity of cell rows
- `cols::Int`: Quantity of cell columns
- `pattern::Symbol`: Spatial layout pattern
- `tms_routing::Symbol`: Cooling pipe routing style
- `tms_encasement::Symbol`: Boundary condition for cooling pipe
- `cell_pitch::Float64`: Centre-to-centre distance between adjacent cells

# Returns
- Instantiated PackGeometry struct containing multi-physics graph definitions
"""
function build_pack_geometry(; 
        rows::Int, 
        cols::Int, 
        pattern::Symbol = :hexagonal, 
        tms_routing::Symbol = :single_row,
        tms_encasement::Symbol = :full,
        cell_pitch::Float64 = 0.023
    )
    
    # Initialise coordinate array and calculate vertical spacing based on grid pattern
    cell_coords = Tuple{Float64, Float64}[]
    D = cell_pitch
    y_step = pattern == :hexagonal ? D * (sqrt(3)/2) : D

    # Iterate through grid limits to plot staggered or flush cell positions
    for r in 1:rows
        for c in 1:cols
            x = (c - 1) * D
            y = (r - 1) * y_step
            
            if pattern == :hexagonal && iseven(r)
                x += D / 2.0
            end
            push!(cell_coords, (x, y))
        end
    end

    # Apply routing and encasement rule engine to determine active vertical fluid levels
    tms_passes = Int[]
    
    if tms_routing == :single_row
        if tms_encasement == :full
            tms_passes = collect(0:rows)
        elseif tms_encasement == :start_bottom
            tms_passes = collect(0:rows-1)
        elseif tms_encasement == :start_top
            tms_passes = collect(1:rows)
        elseif tms_encasement == :internal
            tms_passes = collect(1:rows-1)
        else
            error("Invalid TMS encasement for :single_row.")
        end
    elseif tms_routing == :double_row
        if iseven(rows)
            if tms_encasement == :full
                tms_passes = collect(0:2:rows)
            elseif tms_encasement == :internal
                tms_passes = collect(1:2:(rows-1))
            else
                error("For even number of rows with :double_row routing, only :full and :internal encasements are allowed.")
            end
        else
            if tms_encasement == :start_bottom
                tms_passes = push!(collect(1:2:rows), 0)
                sort!(tms_passes)
            elseif tms_encasement == :start_top
                tms_passes = push!(collect(0:2:(rows-1)), rows)
                sort!(tms_passes)
            else
                error("For uneven number of rows with :double_row routing, only :start_bottom and :start_top encasements are allowed.")
            end
        end
    else
        error("Invalid TMS routing. Choose :single_row or :double_row.")
    end

    # Precalculate spatial coordinates for channel bends between active passes
    u_turns = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    direction = 1
    
    for i in firstindex(tms_passes):(lastindex(tms_passes)-1)
        pass_level = tms_passes[i]
        next_pass = tms_passes[i+1]
        turn_pts = Tuple{Float64, Float64}[]
        
        row_start_block = pass_level + 1
        row_end_block = next_pass
        
        for r in row_start_block : row_end_block
            if direction == 1 
                cx, cy = cell_coords[(r - 1) * cols + cols]
                
                if pattern == :hexagonal
                    if iseven(r) 
                        push!(turn_pts, (cx + (D/2)*cos(-pi/3), cy + (D/2)*sin(-pi/3)))
                        push!(turn_pts, (cx + D/2, cy))
                        push!(turn_pts, (cx + (D/2)*cos(pi/3), cy + (D/2)*sin(pi/3)))
                    else 
                        push!(turn_pts, (cx + D/2, cy))
                    end
                else 
                    if r > row_start_block
                        push!(turn_pts, (cx + D/2, cy - D/2)) 
                    end
                    push!(turn_pts, (cx + D/2, cy))
                    if r < row_end_block
                        push!(turn_pts, (cx + D/2, cy + D/2)) 
                    end
                end
            else 
                cx, cy = cell_coords[(r - 1) * cols + 1]
                
                if pattern == :hexagonal
                    if !iseven(r) 
                        push!(turn_pts, (cx + (D/2)*cos(4pi/3), cy + (D/2)*sin(4pi/3)))
                        push!(turn_pts, (cx - D/2, cy))
                        push!(turn_pts, (cx + (D/2)*cos(2pi/3), cy + (D/2)*sin(2pi/3)))
                    else 
                        push!(turn_pts, (cx - D/2, cy))
                    end
                else 
                    if r > row_start_block
                        push!(turn_pts, (cx - D/2, cy - D/2))
                    end
                    push!(turn_pts, (cx - D/2, cy))
                    if r < row_end_block
                        push!(turn_pts, (cx - D/2, cy + D/2))
                    end
                end
            end
        end
        
        sort!(turn_pts, by = p -> p[2])
        u_turns[i] = turn_pts
        direction *= -1
    end
    
    # Calculate continuous straight channel segments tracing cell boundaries
    all_pts = Tuple{Float64, Float64}[]
    direction = 1
    
    for i in eachindex(tms_passes)
        pass_level = tms_passes[i]
        y_tms = pass_level * y_step - (y_step / 2.0)
        
        pass_min_x = -Inf
        pass_max_x = Inf
        
        if direction == 1 
            if i > 1
                pass_min_x = u_turns[i-1][end][1] 
            end
            if i < lastindex(tms_passes)
                pass_max_x = u_turns[i][1][1]    
            end
        else 
            if i > 1
                pass_max_x = u_turns[i-1][end][1] 
            end
            if i < lastindex(tms_passes)
                pass_min_x = u_turns[i][1][1]    
            end
        end
        
        straight_pts = Tuple{Float64, Float64}[]
        
        for c in 1:cols
            if pass_level > 0
                r = pass_level
                cx, cy = cell_coords[(r - 1) * cols + c]
                if pattern == :hexagonal
                    push!(straight_pts, (cx + (D/2)*cos(pi/3), cy + (D/2)*sin(pi/3)))
                    push!(straight_pts, (cx + (D/2)*cos(2pi/3), cy + (D/2)*sin(2pi/3)))
                else
                    push!(straight_pts, (cx, cy + D/2))
                end
            end
            
            if pass_level < rows
                r = pass_level + 1
                cx, cy = cell_coords[(r - 1) * cols + c]
                if pattern == :hexagonal
                    push!(straight_pts, (cx + (D/2)*cos(-pi/3), cy + (D/2)*sin(-pi/3)))
                    push!(straight_pts, (cx + (D/2)*cos(-2pi/3), cy + (D/2)*sin(-2pi/3)))
                else
                    push!(straight_pts, (cx, cy - D/2))
                end
            end
        end
        
        # Clip straight segments precisely to bend boundaries and order array sequentially
        filter!(p -> (p[1] >= pass_min_x - 1e-4) && (p[1] <= pass_max_x + 1e-4), straight_pts)
        straight_pts = unique(map(p -> (round(p[1], digits=5), round(p[2], digits=5)), straight_pts))
        
        if direction == 1
            sort!(straight_pts, by = p -> p[1])
        else
            sort!(straight_pts, by = p -> p[1], rev=true)
        end
        
        # Extend absolute inlet and outlet nodes mathematically beyond pack boundary
        if i == 1
            if direction == 1
                pushfirst!(straight_pts, (straight_pts[1][1] - D, y_tms))
            else
                pushfirst!(straight_pts, (straight_pts[1][1] + D, y_tms))
            end
        end
        if i == lastindex(tms_passes)
            if direction == 1
                push!(straight_pts, (straight_pts[end][1] + D, y_tms))
            else
                push!(straight_pts, (straight_pts[end][1] - D, y_tms))
            end
        end
        
        append!(all_pts, straight_pts)
        
        if i < lastindex(tms_passes)
            append!(all_pts, u_turns[i])
        end
        
        direction *= -1
    end
    
    # Flatten coordinate array and filter overlapping points to finalise cooling node positions
    tms_coords = Tuple{Float64, Float64}[]
    for p in all_pts
        if isempty(tms_coords) || !isapprox(tms_coords[end][1], p[1], atol=1e-4) || !isapprox(tms_coords[end][2], p[2], atol=1e-4)
            push!(tms_coords, p)
        end
    end
    
    # Link sequential fluid nodes to formulate continuous physical pipe routing
    flow_edges = Tuple{Int, Int, Float64}[]
    for i in firstindex(tms_coords):(lastindex(tms_coords)-1)
        dist = sqrt((tms_coords[i][1] - tms_coords[i+1][1])^2 + (tms_coords[i][2] - tms_coords[i+1][2])^2)
        push!(flow_edges, (i, i+1, dist))
    end

    # Identify adjacent cells within expected pitch distance to build transverse thermal conduction network
    cell_edges = Tuple{Int, Int, Float64}[]
    tolerance = 1e-4
    for i in eachindex(cell_coords)
        for j in (i+1):lastindex(cell_coords)
            dist = sqrt((cell_coords[i][1] - cell_coords[j][1])^2 + (cell_coords[i][2] - cell_coords[j][2])^2)
            
            if abs(dist - D) < tolerance
                row_i = round(Int, cell_coords[i][2] / y_step) + 1
                row_j = round(Int, cell_coords[j][2] / y_step) + 1
                blocked_by_tms = false
                
                # Verify thermal path avoids physical obstruction by cooling ribbon
                if row_i != row_j
                    gap_level = min(row_i, row_j)
                    if gap_level in tms_passes
                        blocked_by_tms = true
                    end
                end
                
                if !blocked_by_tms
                    push!(cell_edges, (i, j, dist))
                end
            end
        end
    end

    # Identify cells within physical contact range of cooling channel to register convective bridges
    convection_edges = Tuple{Int, Int, Float64}[]
    convection_threshold = 0.55 * D 
    for c_idx in eachindex(cell_coords)
        for t_idx in eachindex(tms_coords)
            dist = sqrt((cell_coords[c_idx][1] - tms_coords[t_idx][1])^2 + (cell_coords[c_idx][2] - tms_coords[t_idx][2])^2)
            
            if dist <= convection_threshold
                push!(convection_edges, (c_idx, t_idx, dist))
            end
        end
    end

    # Calculate dynamic angular shadows cast by adjacent geometry to determine ambient exposure fraction
    casing_edges = Tuple{Int, Float64}[]
    shadow_width = pattern == :hexagonal ? pi / 3.0 : pi / 2.0 
    
    for c_idx in eachindex(cell_coords)
        cx, cy = cell_coords[c_idx]
        arcs = Tuple{Float64, Float64}[]
        
        for (i, j, dist) in cell_edges
            if i == c_idx || j == c_idx
                other = i == c_idx ? j : i
                ox, oy = cell_coords[other]
                ang = atan(oy - cy, ox - cx)
                push!(arcs, (ang - shadow_width/2, ang + shadow_width/2))
            end
        end
        
        for (c, t, dist) in convection_edges
            if c == c_idx
                tx, ty = tms_coords[t]
                ang = atan(ty - cy, tx - cx)
                push!(arcs, (ang - shadow_width/2, ang + shadow_width/2))
            end
        end
        
        split_arcs = Tuple{Float64, Float64}[]
        for (s, e) in arcs
            s_mod, e_mod = mod(s, 2pi), mod(e, 2pi)
            if s_mod > e_mod
                push!(split_arcs, (s_mod, 2pi), (0.0, e_mod))
            else
                push!(split_arcs, (s_mod, e_mod))
            end
        end
        sort!(split_arcs, by = x -> x[1])
        
        # Merge overlapping shadow arcs into continuous blocked regions
        covered_length = 0.0
        if !isempty(split_arcs)
            merged = [split_arcs[1]]
            for k in 2:lastindex(split_arcs)
                curr_s, curr_e = split_arcs[k]
                prev_s, prev_e = merged[end]
                if curr_s <= prev_e + 1e-5
                    merged[end] = (prev_s, max(prev_e, curr_e))
                else
                    push!(merged, (curr_s, curr_e))
                end
            end
            
            covered_length = sum(e - s for (s, e) in merged)
        end
        
        exposed_fraction = max(0.0, 1.0 - (covered_length / (2pi)))
        push!(casing_edges, (c_idx, exposed_fraction))
    end

    # Wire parallel and series connections sequentially to generate electrical topology graph
    electrical_edges = Tuple{Int, Int, Symbol}[]
    get_idx(r, c) = (r - 1) * cols + c
    
    for r in 1:rows
        for c in 1:cols
            idx = get_idx(r, c)
            if c < cols
                push!(electrical_edges, (idx, get_idx(r, c+1), :parallel))
            end
            if r < rows
                push!(electrical_edges, (idx, get_idx(r+1, c), :series))
            end
        end
    end

    return PackGeometry(cell_coords, tms_coords, cell_edges, convection_edges, casing_edges, flow_edges, electrical_edges)
end