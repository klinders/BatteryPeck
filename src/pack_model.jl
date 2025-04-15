module PackModel

using ..CellModel
export battery_pack_ecm!

function battery_pack_ecm!(du, u, p, t)
    for i in 1:p.num_cells
        # get the equations per cell
        u_per_cell = size(u,1)/p.num_cells 
        
        # Index the state and change in state
        j = Int(u_per_cell*(i-1) + 1)
        u_cell = @view u[j:j+1]
        du_cell = @view du[j:j+1]

        p.cell_params.I = p.I

        # Get the cell model for each cell
        battery_cell_ecm!(du_cell, u_cell, p.cell_params, t)
    end
end

end  # module PackModel