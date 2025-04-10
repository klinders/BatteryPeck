module PackModel

using ..CellModel
export battery_pack_ecm!

function battery_pack_ecm!(du, u, p, t)
    num_cells = p[:num_cells]
    for i in 1:num_cells
        cell_params = Dict(:Q => p[:Q], :R0 => p[:R0], :R1 => p[:R1], :C1 => p[:C1], :I => p[:I])
        du[i:i+1] .= battery_cell_ecm!(du[i:i+1], u[i:i+1], cell_params, t)
    end
end

end  # module PackModel