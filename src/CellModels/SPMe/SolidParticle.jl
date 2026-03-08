# =====================================================================================================================
# SolidParticle.jl
#
# Particle equations for SPMe
# Source: https://doi.org/10.1016/j.apm.2022.12.009
# =====================================================================================================================

# Import package
using ModelingToolkit

function SolidParticle(; name, p::SolidParticleParameters, g)
    # Independent variables
    @parameters begin
        t # Time
    end

    # Time derivative
    Dt = Differential(t)

    # Constants
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant

    # Components
    @named J = RealInput() # Volumetric current density
    @named T = RealInput() # Ambient temperature

    # Time-dependent state variables
    @variables begin
        # Subdivide concentration into Nᵣ cells, and initialise at c₀ (Eq. 1d)
        (c(t))[1:g.Nᵣ] = repeat([p.c₀],g.Nᵣ) # Local concentration
        c_avr(t)                             # Average concentration
        c_surf(t)                            # Surface concentration
        U₀(t)                                # Open-circuit potential
        z(t)                                 # Stoichiometry
        ϕ̄ₛ(t)                                 # 
    end

    # Retrieve FVM geometry
    Δr,r,Vᵢ,Aₗ,Aᵣ = g.Δr, g.r, g.Vᵢ, g.Aₗ, g.Aᵣ

    # Face diffusivities (see "helpers.jl" for explanation and source)
    # Outer most cells have only one neighbour
    Dₗ = [nothing; [D_face(p.Dₖ(c[i-1]), p.Dₖ(c[i]), Δr, Δr) for i in 2:g.Nᵣ]] # Left neighbours
    Dᵣ = [[D_face(p.Dₖ(c[i]), p.Dₖ(c[i+1]), Δr, Δr) for i in 1:g.Nᵣ-1]; nothing] # Right neighbours

    # Particle equations
    eqns = [
        # FVM gradient = (c[i+1]-c[i])/Δr[i]   or  (c[i]-c[i-1])/Δr[i]
        # FVM 3D spherical divergence = (flux_right*area_right - flux_left*area_left)/volume_cell

        c_avr ~ sum(c)/g.Nᵣ                # Average concentration
        c_surf ~ 1.5*c[end] - 0.5*c[end-1] # Surface concentration
        z ~ c_surf/p.c₊                    # Stoichiometry
        U₀ ~ p.Uₖ(z)                        # Open-circuit potential

        # Eqs. 1a/3a;
        # Inner most cell (no flux from/to left neighbour)
        Dt(c[1]) ~ (Dᵣ[1] * (c[2]-c[1]) * Aᵣ[1] / Δr) / Vᵢ[1]

        # Middle cells
        [Dt(c[i]) ~ 
            (
            Dᵣ[i] * (c[i+1]-c[i]) * Aᵣ[i] / Δr - 
            Dₗ[i] * (c[i]-c[i-1]) * Aₗ[i] / Δr
            ) / Vᵢ[i] 
        for i in 2:g.Nᵣ-1]...

       # Eqs. 1c/3c;
        # Outer most cell
        Dt(c[end]) ~ 
            (
            -Aᵣ[end] * J.u / p.aₖ / F 
            - Dₗ[end] * (c[end]-c[end-1]) * Aₗ[end] / Δr
            ) / Vᵢ[end]
    ]

    # Construct ODESystem with equations and child components
    System(eqns, t; name=name, systems=[J, T])
end