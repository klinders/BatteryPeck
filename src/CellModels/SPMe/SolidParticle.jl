# =====================================================================================================================
# SolidParticle.jl
# Particle equations for SPMe
# Source: https://doi.org/10.1016/j.apm.2022.12.009
# =====================================================================================================================

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

    Dt = Differential(t)
    R = 8.314 
    F = 96485 

    @named J = RealInput() 
    @named T = RealInput(guess=298.15) 

    @variables begin
        (c(t))[1:g.Nᵣ] = repeat([p.c₀],g.Nᵣ) 
        (D(t))[1:g.Nᵣ], [guess=ones(g.Nᵣ)*1e-14]
        (σ(t))[1:g.Nᵣ], [guess=ones(g.Nᵣ)]
        D_r(t), [guess=1e-14]
        c_avr(t), [guess=p.c₀]
        c_r(t), [guess=p.c₀]
        c_surf(t), [guess=p.c₀]
        U₀(t), [guess=p.Uₖ(p.c₀/p.c₊)]
        z(t), [guess=p.c₀/p.c₊]
        ϕ̄ₛ(t), [guess=0.0]
        ϵₛ(t) = p.ϵₛ # active material volume fraction
        aₖ(t), [guess=3*p.ϵₛ/p.Rₖ] # specific surface area
    end

    Δr,r,Vᵢ,Aₗ,Aᵣ = g.Δr, g.r, g.Vᵢ, g.Aₗ, g.Aᵣ

    θ_M = p.Ω / (R * T.u) * (2 * p.Ω * p.E) / (9 * (1 - p.ν))
    c₀_cr = 0.0

    Dₗ = [nothing; [D_face(D[i-1],D[i],Δr,Δr) for i in 2:g.Nᵣ]] 
    Dᵣ = [[D_face(D[i], D[i+1],Δr,Δr) for i in 1:g.Nᵣ-1]; nothing] 

    eqns = [
        [σ[i] ~ 1 + θ_M * (c[i] - c₀_cr) for i in 1:g.Nᵣ]...,
        # FIXED: Removed T.u from p.Dₖ call to match our 1-argument parameter structure
        [D[i] ~ p.Dₖ(c[i])*σ[i] for i in 1:g.Nᵣ]...,
        D_r ~ sum([D[i]*Vᵢ[i] for i in 1:g.Nᵣ])/sum(Vᵢ),
        
        c_avr ~ sum(c)/g.Nᵣ,
        c_r ~ sum([c[i]*Vᵢ[i] for i in 1:g.Nᵣ])/sum(Vᵢ),
        c_surf ~ 1.5*c[end] - 0.5*c[end-1], 
        z ~ c_surf/p.c₊,
        U₀ ~ p.Uₖ(z),
        aₖ ~ 3*ϵₛ/p.Rₖ,

        Dt(c[1]) ~ (Dᵣ[1]*Aᵣ[1]*(c[2] - c[1])/Δr)/Vᵢ[1],
        [Dt(c[i]) ~ (Dᵣ[i]*Aᵣ[i]*(c[i+1] - c[i])/Δr - Dₗ[i]*Aₗ[i]*(c[i] - c[i-1])/Δr)/Vᵢ[i] for i in 2:g.Nᵣ-1]...,
        Dt(c[end]) ~ (-Aᵣ[end]*J.u/F - Dₗ[end]*Aₗ[end]*(c[end] - c[end-1])/Δr)/Vᵢ[end]
    ]

    System(eqns, t; name=name, systems=[J, T])
end