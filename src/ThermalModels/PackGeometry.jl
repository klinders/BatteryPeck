# ==============================================================================
# PackGeometry.jl
# Define battery pack spatial layout and auto generate multi-physics graph edges
# ==============================================================================

export PackGeometry, build_pack_geometry

struct PackGeometry
    # Stores spatial coordinates of cell centres
    cell_coords::Vector{Tuple{Float64, Float64}}
    # Stores spatial coordinates of cooling channel nodes
    tms_coords::Vector{Tuple{Float64, Float64}}
    
    # Defines thermal conduction paths between adjacent cells
    cell_edges::Vector{Tuple{Int, Int, Float64}}
    # Defines thermal bridges between cells and cooling ribbon
    convection_edges::Vector{Tuple{Int, Int, Float64}}
    # Defines exposed surface area fraction for ambient cooling
    casing_edges::Vector{Tuple{Int, Float64}} 
    
    # Defines sequential fluid routing between cooling nodes
    flow_edges::Vector{Tuple{Int, Int, Float64}}
    
    # Defines series and parallel electrical connections
    electrical_edges::Vector{Tuple{Int, Int, Symbol}}
end

"""
    build_pack_geometry(; rows, cols, pattern, tms_routing, tms_encasement, cell_pitch)

Generates spatial layout and multi-physics graph edges for battery pack.

# Arguments
- `rows::Int`: Quantity of cell rows
- `cols::Int`: Quantity of cell columns
- `pattern::Symbol`: Spatial layout pattern (`:hexagonal` or `:square`)
- `tms_routing::Symbol`: Cooling pipe routing style (`:single_row` or `:double_row`)
- `tms_encasement::Symbol`: Boundary condition for cooling pipe (`:full`, `:internal`, `:start_bottom`, or `:start_top`)
- `cell_pitch::Float64`: Centre-to-centre distance between adjacent cells [m]

# Returns
- `PackGeometry`: Struct containing node coordinates and connecting edges for thermal, fluid, and electrical graphs
"""
function build_pack_geometry(; 
        rows::Int, 
        cols::Int, 
        pattern::Symbol = :hexagonal, 
        tms_routing::Symbol = :single_row,
        tms_encasement::Symbol = :full,
        cell_pitch::Float64 = 0.023
    )
    
    # Initialise cell coordinate array
    cell_coords = Tuple{Float64, Float64}[]
    D = cell_pitch
    
    # Calculate vertical row spacing based on layout pattern
    y_step = pattern == :hexagonal ? D * (sqrt(3)/2) : D

    # Iterate through rows and columns to plot cell positions
    for r in 1:rows
        for c in 1:cols
            x = (c - 1) * D
            y = (r - 1) * y_step
            
            # Offset even rows for hexagonal packing to nest cells
            if pattern == :hexagonal && iseven(r)
                x += D / 2.0
            end
            push!(cell_coords, (x, y))
        end
    end

    # Initialise array to store vertical levels of cooling ribbon
    tms_passes = Int[]
    
    # Apply routing and encasement rule engine to determine active fluid levels
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

    # Precalculate spatial coordinates for channel bends
    u_turns = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    
    # Left to right flow is represented by positive direction
    direction = 1
    
    # Iterate through passes to generate turn geometry between levels
    for i in firstindex(tms_passes):(lastindex(tms_passes)-1)
        pass_level = tms_passes[i]
        next_pass = tms_passes[i+1]
        turn_pts = Tuple{Float64, Float64}[]
        
        row_start_block = pass_level + 1
        row_end_block = next_pass
        
        for r in row_start_block : row_end_block
            if direction == 1 
                cx, cy = cell_coords[(r - 1) * cols + cols]
                
                # Map specific angular points to wrap ribbon around right side cells
                if pattern == :hexagonal
                    if iseven(r) 
                        # Plot three points at -60, 0, and +60 degrees to trace curved right edge of staggered cell
                        push!(turn_pts, (cx + (D/2)*cos(-pi/3), cy + (D/2)*sin(-pi/3)))
                        push!(turn_pts, (cx + D/2, cy))
                        push!(turn_pts, (cx + (D/2)*cos(pi/3), cy + (D/2)*sin(pi/3)))
                    else 
                        # Plot single apex point at 0 degrees for flush cell
                        push!(turn_pts, (cx + D/2, cy))
                    end
                else 
                    # Trace right-angled corners for square grid layout
                    if r > row_start_block
                        # Plot bottom right corner of square bounding box
                        push!(turn_pts, (cx + D/2, cy - D/2)) 
                    end
                    # Plot centre right edge
                    push!(turn_pts, (cx + D/2, cy))
                    if r < row_end_block
                        # Plot top right corner of square bounding box
                        push!(turn_pts, (cx + D/2, cy + D/2)) 
                    end
                end
            else 
                cx, cy = cell_coords[(r - 1) * cols + 1]
                
                # Map specific angular points to wrap ribbon around left side cells
                if pattern == :hexagonal
                    if !iseven(r) 
                        # Plot three points at 240, 180, and 120 degrees to trace curved left edge of staggered cell
                        push!(turn_pts, (cx + (D/2)*cos(4pi/3), cy + (D/2)*sin(4pi/3)))
                        push!(turn_pts, (cx - D/2, cy))
                        push!(turn_pts, (cx + (D/2)*cos(2pi/3), cy + (D/2)*sin(2pi/3)))
                    else 
                        # Plot single apex point at 180 degrees for flush cell
                        push!(turn_pts, (cx - D/2, cy))
                    end
                else 
                    # Trace right-angled corners for square grid layout
                    if r > row_start_block
                        # Plot bottom left corner of square bounding box
                        push!(turn_pts, (cx - D/2, cy - D/2))
                    end
                    # Plot centre left edge
                    push!(turn_pts, (cx - D/2, cy))
                    if r < row_end_block
                        # Plot top left corner of square bounding box
                        push!(turn_pts, (cx - D/2, cy + D/2))
                    end
                end
            end
        end
        
        # Sort bend points purely vertically to ensure smooth flow path
        sort!(turn_pts, by = p -> p[2])
        u_turns[i] = turn_pts
        direction *= -1
    end
    
    # Initialise array for all continuous channel coordinates
    all_pts = Tuple{Float64, Float64}[]
    direction = 1
    
    # Iterate through passes to generate straight channel segments
    for i in eachindex(tms_passes)
        pass_level = tms_passes[i]
        y_tms = pass_level * y_step - (y_step / 2.0)
        
        pass_min_x = -Inf
        pass_max_x = Inf
        
        # Determine strict bounding box for current pass to prevent overlap
        if direction == 1 
            if i > 1
                # Anchor to top of left bend
                pass_min_x = u_turns[i-1][end][1] 
            end
            if i < lastindex(tms_passes)
                # Anchor to bottom of right bend
                pass_max_x = u_turns[i][1][1]     
            end
        else 
            if i > 1
                # Anchor to top of right bend
                pass_max_x = u_turns[i-1][end][1] 
            end
            if i < lastindex(tms_passes)
                # Anchor to bottom of left bend
                pass_min_x = u_turns[i][1][1]     
            end
        end
        
        straight_pts = Tuple{Float64, Float64}[]
        
        # Calculate straight segments tracing along top and bottom of cell centres
        for c in 1:cols
            if pass_level > 0
                r = pass_level
                cx, cy = cell_coords[(r - 1) * cols + c]
                if pattern == :hexagonal
                    # Trace upper curve of cell at 60 and 120 degrees
                    push!(straight_pts, (cx + (D/2)*cos(pi/3), cy + (D/2)*sin(pi/3)))
                    push!(straight_pts, (cx + (D/2)*cos(2pi/3), cy + (D/2)*sin(2pi/3)))
                else
                    # Trace flat top edge of square bounding box
                    push!(straight_pts, (cx, cy + D/2))
                end
            end
            
            if pass_level < rows
                r = pass_level + 1
                cx, cy = cell_coords[(r - 1) * cols + c]
                if pattern == :hexagonal
                    # Trace lower curve of cell at -60 and -120 degrees
                    push!(straight_pts, (cx + (D/2)*cos(-pi/3), cy + (D/2)*sin(-pi/3)))
                    push!(straight_pts, (cx + (D/2)*cos(-2pi/3), cy + (D/2)*sin(-2pi/3)))
                else
                    # Trace flat bottom edge of square bounding box
                    push!(straight_pts, (cx, cy - D/2))
                end
            end
        end
        
        # Clip straight segments precisely to bend boundaries
        filter!(p -> (p[1] >= pass_min_x - 1e-4) && (p[1] <= pass_max_x + 1e-4), straight_pts)
        
        # Remove rounding errors to ensure perfect overlapping points
        straight_pts = unique(map(p -> (round(p[1], digits=5), round(p[2], digits=5)), straight_pts))
        
        # Order points sequentially according to flow direction
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
        
        # Combine straight segments into main array
        append!(all_pts, straight_pts)
        
        # Append corresponding U-turn to complete pass continuity
        if i < lastindex(tms_passes)
            append!(all_pts, u_turns[i])
        end
        
        direction *= -1
    end
    
    # Flatten coordinate array and remove overlapping spatial points
    tms_coords = Tuple{Float64, Float64}[]
    for p in all_pts
        if isempty(tms_coords) || !isapprox(tms_coords[end][1], p[1], atol=1e-4) || !isapprox(tms_coords[end][2], p[2], atol=1e-4)
            push!(tms_coords, p)
        end
    end
    
    # Link sequential fluid nodes to form continuous physical pipe
    flow_edges = Tuple{Int, Int, Float64}[]
    for i in firstindex(tms_coords):(lastindex(tms_coords)-1)
        dist = sqrt((tms_coords[i][1] - tms_coords[i+1][1])^2 + (tms_coords[i][2] - tms_coords[i+1][2])^2)
        push!(flow_edges, (i, i+1, dist))
    end

    # Identify adjacent cells for transverse thermal conduction network
    cell_edges = Tuple{Int, Int, Float64}[]
    tolerance = 1e-4
    for i in eachindex(cell_coords)
        for j in (i+1):lastindex(cell_coords)
            dist = sqrt((cell_coords[i][1] - cell_coords[j][1])^2 + (cell_coords[i][2] - cell_coords[j][2])^2)
            
            # Check if physical gap between cells matches expected pitch
            if abs(dist - D) < tolerance
                row_i = round(Int, cell_coords[i][2] / y_step) + 1
                row_j = round(Int, cell_coords[j][2] / y_step) + 1
                blocked_by_tms = false
                
                # Verify thermal path is not physically blocked by cooling ribbon
                if row_i != row_j
                    gap_level = min(row_i, row_j)
                    if gap_level in tms_passes
                        blocked_by_tms = true
                    end
                end
                
                # Register valid conduction edge
                if !blocked_by_tms
                    push!(cell_edges, (i, j, dist))
                end
            end
        end
    end

    # Identify cells within physical contact range of cooling channel
    convection_edges = Tuple{Int, Int, Float64}[]
    convection_threshold = 0.55 * D 
    for c_idx in eachindex(cell_coords)
        for t_idx in eachindex(tms_coords)
            dist = sqrt((cell_coords[c_idx][1] - tms_coords[t_idx][1])^2 + (cell_coords[c_idx][2] - tms_coords[t_idx][2])^2)
            
            # Register valid convective bridge
            if dist <= convection_threshold
                push!(convection_edges, (c_idx, t_idx, dist))
            end
        end
    end

    # Assign dynamic shadow angle based on physical packing layout
    casing_edges = Tuple{Int, Float64}[]
    shadow_width = pattern == :hexagonal ? pi / 3.0 : pi / 2.0 
    
    for c_idx in eachindex(cell_coords)
        cx, cy = cell_coords[c_idx]
        arcs = Tuple{Float64, Float64}[]
        
        # Calculate angular shadow cast by adjacent cells
        for (i, j, dist) in cell_edges
            if i == c_idx || j == c_idx
                other = i == c_idx ? j : i
                ox, oy = cell_coords[other]
                ang = atan(oy - cy, ox - cx)
                push!(arcs, (ang - shadow_width/2, ang + shadow_width/2))
            end
        end
        
        # Calculate angular shadow cast by adjacent cooling nodes
        for (c, t, dist) in convection_edges
            if c == c_idx
                tx, ty = tms_coords[t]
                ang = atan(ty - cy, tx - cx)
                push!(arcs, (ang - shadow_width/2, ang + shadow_width/2))
            end
        end
        
        # Convert angles to positive values and handle circular wrap around
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
            
            # Sum total radians blocked by merged shadows
            covered_length = sum(e - s for (s, e) in merged)
        end
        
        # Calculate percentage of cell perimeter exposed to ambient air
        exposed_fraction = max(0.0, 1.0 - (covered_length / (2pi)))
        push!(casing_edges, (c_idx, exposed_fraction))
    end

    # Generate electrical topology graph for potential extension
    electrical_edges = Tuple{Int, Int, Symbol}[]
    get_idx(r, c) = (r - 1) * cols + c
    
    # Wire parallel and series connections sequentially
    for r in 1:rows
        for c in 1:cols
            idx = get_idx(r, c)
            if c < cols
                push!(electrical_edges, (idx, get_idx(r, c+1), :series))
            end
            if r < rows
                push!(electrical_edges, (idx, get_idx(r+1, c), :parallel))
            end
        end
    end

    return PackGeometry(cell_coords, tms_coords, cell_edges, convection_edges, casing_edges, flow_edges, electrical_edges)
end