using ModelingToolkit

"""
    SolidParticle(; name, p::SolidParticleParameters, g)

Create a ModelingToolkit system for solid-state lithium diffusion in a battery electrode.

Models radial diffusion of lithium ions within spherical electrode particles using the finite
volume method. Computes surface concentration, stoichiometry, and open-circuit potential.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SolidParticleParameters`: Electrode material parameters
- `g`: FVM geometry object with node locations and volumes

# Input Ports
- `J`: Surface current density (A/m²)
- `T`: Temperature (K)

# Output Variables
- `c`: Concentration profile across particles (mol/m³)
- `c_surf`: Surface lithium concentration (mol/m³)
- `z`: Stoichiometry = c_surf/c_max
- `U₀`: Open-circuit potential (V)

# Notes
Uses second-order accurate finite volume discretization with ghost nodes for boundary conditions.
"""
function SolidParticle(; name, p::SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant

    @named J = RealInput()
    @named T = RealInput()

    # Time derivative operator
    Dt = Differential(t)
    @show(p.c₀)

    @variables begin
        # I am adding two ghost nodes for the boundary conditions
        (c(t))[1:g.Nᵣ] = repeat([p.c₀],g.Nᵣ)
        (D(t))[1:g.Nᵣ]
        (σ(t))[1:g.Nᵣ]
        D_r(t)
        c_avr(t)
        c_r(t)
        c_surf(t)
        U₀(t)
        z(t)
        ϕ̄ₛ(t)
    end

    # Discretized equations
    Δr,r,Vᵢ,Aₗ,Aᵣ = g.Δr, g.r, g.Vᵢ, g.Aₗ, g.Aᵣ

    θ_M = p.Ω/ (R * T.u) * (2 * p.Ω * p.E) / (9 * (1 - p.ν))
    c₀_cr = 0.0

    D_f = [p.Dₖ(c[i], T.u) for i in 1:g.Nᵣ] # Diffusivities at cell centers
    Dₗ = [nothing; [D_face(D[i-1],D[i],Δr,Δr) for i in 2:g.Nᵣ]] # Left diffusivities
    Dᵣ = [[D_face(D[i], D[i+1],Δr,Δr) for i in 1:g.Nᵣ-1]; nothing] # Right diffusivities

    eqns = [
        # Diffusion with stress
        [σ[i] ~ 1 + θ_M * (c[i] - c₀_cr) for i in 1:g.Nᵣ]...
        [D[i] ~ p.Dₖ(c[i], T.u)*σ[i] for i in 1:g.Nᵣ]...
        D_r ~ sum([D[i]*Vᵢ[i] for i in 1:g.Nᵣ])/sum(Vᵢ)
        
        c_avr ~ sum(c)/g.Nᵣ
        c_r ~ sum([c[i]*Vᵢ[i] for i in 1:g.Nᵣ])/sum(Vᵢ)
        c_surf ~ 1.5*c[end] - 0.5*c[end-1] # Surface concentration
        z ~ c_surf/p.c₊ # Stoichiometry
        U₀ ~ p.Uₖ(z) # Open-circuit potential

        # Boundary condition center
        Dt(c[1]) ~ (Dᵣ[1]*Aᵣ[1]*(c[2] - c[1])/Δr)/Vᵢ[1]

        # Internal nodes
        [Dt(c[i]) ~ (Dᵣ[i]*Aᵣ[i]*(c[i+1] - c[i])/Δr - Dₗ[i]*Aₗ[i]*(c[i] - c[i-1])/Δr)/Vᵢ[i] for i in 2:g.Nᵣ-1]...

        # Boundary condition edge
        Dt(c[end]) ~ (-Aᵣ[end]*J.u/F - Dₗ[end]*Aₗ[end]*(c[end] - c[end-1])/Δr)/Vᵢ[end]

    ]

    System(eqns,t; name=name,systems=[J, T])
end