module PackModel

using ..CellModel
export battery_pack_ecm!

function battery_pack_ecm!(du, u, p, t)
    for i in 1:p.num_cells       
        u_per_cell = size(u,1)/p.num_cells # equations per cell
        j = Int(u_per_cell*(i-1) + 1)
        
        u_cell = @view u[j:j+1]
        du_cell = @view du[j:j+1]

        battery_cell_ecm!(du_cell, u_cell, p.cell_params, t)
    end
end

end  # module PackModel