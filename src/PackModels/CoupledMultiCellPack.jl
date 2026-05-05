# ==============================================================================
# CoupledMultiCellPack.jl
# Automatically wires SPMe cells and LPTN (Virtual Series)
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using BatteryToolkit
export CoupledMultiCellPack

"""
Construct fully coupled electrical and thermal multi-cell pack.
"""
function CoupledMultiCellPack(; name, geom::PackGeometry, therm_params::PackParameters, elec_params::BatteryParameters)
    
    num_cells = length(geom.cell_coords)
    
    @parameters begin
        t
        
        # Define boundary inputs as parameters
        Pin = 0.0
        Iin = 0.0
    end
    
    @variables begin
        V(t)
        I(t)
    end
    
    spme_cells = [SPMe(name=(i==1 ? :cell : Symbol("spme_$i")), params=elec_params) for i in 1:num_cells]
    thermal_pack = build_pack_system(:thermal_pack, geom, therm_params)
    therm_cells = [getproperty(thermal_pack, Symbol("cell_$i")) for i in 1:num_cells]
    
    @named ground = Ground()
    # Create an independent current source for every cell
    cell_sources = [Current(name=Symbol("source_$i")) for i in 1:num_cells]
    
    eqs = Equation[
        I ~ Iin + Pin/V
        # V_pack is strictly the mathematical sum of all cells
        V ~ sum(c.v for c in spme_cells)
    ]
    
    # Virtual Series Binding
    for i in 1:num_cells
        # Ground every cell independently to eliminate floating potentials
        push!(eqs, connect(ground.g, spme_cells[i].n))
        push!(eqs, connect(cell_sources[i].p, spme_cells[i].p))
        push!(eqs, connect(cell_sources[i].n, ground.g))
        
        # Force every independent source to draw the exact global series current
        push!(eqs, cell_sources[i].I.u ~ I)
        
        # Thermal coupling
        push!(eqs, therm_cells[i].Q_volumetric_in.u ~ spme_cells[i].Q_total)
        push!(eqs, spme_cells[i].T.u ~ therm_cells[i].core_cap.T)
    end
    
    sys = ODESystem(eqs, t, [V, I], [Pin, Iin]; systems=[spme_cells..., thermal_pack, ground, cell_sources...], name=name)
    return sys
end