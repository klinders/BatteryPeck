# =====================================================================================================================
# Electrolyte.jl
#
# Electrolyte equations for SPMe
# Source: https://doi.org/10.1016/j.apm.2022.12.009
# =====================================================================================================================

# Import package
using ModelingToolkit

function Electrolyte(; name, p::ElectrolyteParameters, g)
    # Independent variables
    @parameters begin
        t # Time
    end

    Δx = g.Δx
    Δxₗ = g.Δxₗ
    Δxᵣ = g.Δxᵣ
    x = g.x_centers
    
    # Time derivative
    Dt = Differential(t)
    
    # Constants
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant

    # Components
    @named i_app = RealInput() # Applied current density in electrolyte
    @named T = RealInput() # Ambient temperature

    # Time-dependent state variables
    @variables begin
        # Eq. 4a LHS; [Porosity] * [Electrolyte concentration]
        # Subdivide into total number of subdivisions over all regions to match FVM cells
        (ϵcₑ(t))[1:g.Nₜ] = [
            # Scalar regions (single porosity scalar value per region)
            # Initialise at c₀ (Eq. 4c)
            [p.ϵₙ*p.c₀ for _ in g.ixₙ]..., # Negative electrode
            [p.ϵₛ*p.c₀ for _ in g.ixₛ]..., # Separator
            [p.ϵₚ*p.c₀ for _ in g.ixₚ]..., # Positive electrode
        ]

        # Initialise porosity and define as symbolic state variable
        ϵₙ(t) = p.ϵₙ # Negative electrode
        ϵₛ(t) = p.ϵₛ # Separator
        ϵₚ(t) = p.ϵₚ # Positive electrode

        # Define concentration as symbolic state variable (calculated in "eqns" below)
        (cₑ(t))[1:g.Nₜ]
    end

    # Eq. 54; Spatial derivative of applied current density in electrolyte
    function iₑ(x)
        # Negative electrode
        if x <= p.Lₙ
            return i_app.u/p.Lₙ
        # Separator
        elseif x <= p.Lₙ + p.Lₛ
            return 0.0
        # Positive electrode
        else
            return -i_app.u/p.Lₚ
        end
    end

    # Concatenate porosity vectors
    ϵ = [
        [ϵₙ for _ in g.ixₙ] # Negative electrode
        [ϵₛ for _ in g.ixₛ] # Separator
        [ϵₚ for _ in g.ixₚ] # Positive electrode
    ]

    # Concatenate Bruggeman coefficients vectors
    b = [
        [p.bₙ for _ in g.ixₙ] # Negative electrode
        [p.bₛ for _ in g.ixₛ] # Separator
        [p.bₚ for _ in g.ixₚ] # Positive electrode
    ]
    
    # Retrieve FVM geometry
    Δx = g.Δx       # Cell widths
    Δxₗ = g.Δxₗ      # Left neighbour distances
    Δxᵣ = g.Δxᵣ     # Right neighbour distances
    x = g.x_centers # Cell centres
    
    # Eq. 4a RHS; Effective diffusivity (transport efficiency, a.k.a. inverse MacMullin number, calculated using Bruggeman coefficient)
    Dᵢ = [p.Dₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:g.Nₜ]
    # Face diffusivities (see "helpers.jl" for explanation and source)
    # Outer most cells have only one neighbour
    Dₗ = [nothing, [D_face(Dᵢ[i-1], Dᵢ[i], Δx[i-1], Δx[i]) for i in 2:g.Nₜ]...] # Left neighbours
    Dᵣ = [[D_face(Dᵢ[i], Dᵢ[i+1], Δx[i], Δx[i+1]) for i in 1:g.Nₜ-1]..., nothing] # Right neighbours

    # Electrolyte equations
    eqns = [
        # FVM gradient = (c[i+1]-c[i])/Δx[i]   or  (c[i]-c[i-1])/Δx[i]
        # FVM 1D divergence = (flux_right - flux_left)/width_cell

        # Porosity (assumed constant, without side reactions)
        Dt(ϵₙ) ~ 0, # Negative electrode
        Dt(ϵₛ) ~ 0, # Separator
        Dt(ϵₚ) ~ 0, # Positive electrode
        
        # Concentration in electrolyte
        [cₑ[i] ~ ϵcₑ[i]/ϵ[i] for i in 1:g.Nₜ]...,

        # Eq. 4a;
        # Left most cell (no flux from/to left neighbour)
        Dt(ϵcₑ[1]) ~ 
            (
            Dᵣ[1] * (cₑ[2]-cₑ[1]) / Δxᵣ[1] 
            + (1 - p.t₊(cₑ[1])) * iₑ(x[1]) * Δx[1] / F
            ) / Δx[1], 

        # Middle cells
        [Dt(ϵcₑ[i]) ~ 
            (
            Dᵣ[i] * (cₑ[i+1]-cₑ[i]) / Δxᵣ[i] - 
            Dₗ[i] * (cₑ[i]-cₑ[i-1]) / Δxₗ[i] + 
            (1 - p.t₊(cₑ[i])) * iₑ(x[i]) * Δx[i] / F
            ) / Δx[i]
        for i in 2:g.Nₜ-1]...,

        # Right most cell (no flux from/to right neighbour)
        Dt(ϵcₑ[end]) ~ 
            (
            -Dₗ[end] * (cₑ[end]-cₑ[end-1]) / Δxₗ[end]
            + (1 - p.t₊(cₑ[end])) * iₑ(x[end]) * Δx[end] / F
            ) / Δx[end], 
    ]

    # Construct ODESystem with equations and child components
    System(eqns, t; name=name, systems=[i_app, T])
end