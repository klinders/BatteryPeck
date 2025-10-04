using ModelingToolkit
include("../ParameterSets/Base.jl")

# Harmonic mean for face diffusivity
# From:
# http://dx.doi.org/10.1149/2.0291607jes
function D_face(Dleft, Dright, Δxleft, Δxright)
    # Harmonic mean of left and right diffusivities
    return (Δxleft + Δxright)/(Δxleft/Dleft + Δxright/Dright)
end

function SolidParticle(; name, p::SolidParticleParameters, Nᵣ=10)
    
    @parameters begin
        t
    end

    @constants begin
        R = 8.314 # Universal gas constant
        F = 96485 # Faraday's constant
        T = 298 # Temperature
    end

    @named J = RealInput()

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        # I am adding two ghost nodes for the boundary conditions
        (c(t))[1:Nᵣ] = repeat([p.c₀],Nᵣ)
        c_avr(t)
        c_surf(t)
        U₀(t)
        z(t)

    end

    # Discretized equations
    Δr = p.Rₖ/Nᵣ
    r = ([Δr*(i-0.5) for i in 1:Nᵣ]) # Radial positions
    Vᵢ = 4/3*π*[(r[i] + Δr/2)^3 - (r[i] - Δr/2)^3 for i in 1:Nᵣ]
    Aₗ = 4*π*[(r[i] - Δr/2)^2 for i in 1:Nᵣ]
    Aᵣ = 4*π*[(r[i] + Δr/2)^2 for i in 1:Nᵣ]
    Dₗ = [nothing, [D_face(p.Dₖ(c[i-1]), p.Dₖ(c[i]),Δr,Δr) for i in 2:Nᵣ]...] # Left diffusivities
    Dᵣ = [[D_face(p.Dₖ(c[i]), p.Dₖ(c[i+1]),Δr,Δr) for i in 1:Nᵣ-1]..., nothing] # Left diffusivities

    eqns = [
        c_avr ~ sum(c)/Nᵣ
        c_surf ~ c[Nᵣ] # Surface concentration
        z ~ c_surf/p.c₊ # Stoichiometry
        U₀ ~ p.Uₖ(z) # Open-circuit potential

        # Boundary condition center
        Dt(c[1]) ~ (Dᵣ[1]*Aᵣ[1]*(c[2] - c[1])/Δr)/Vᵢ[1]

        # Internal nodes
        [Dt(c[i]) ~ (Dᵣ[i]*Aᵣ[i]*(c[i+1] - c[i])/Δr - Dₗ[i]*Aₗ[i]*(c[i] - c[i-1])/Δr)/Vᵢ[i] for i in 2:Nᵣ-1]...

        # Boundary condition edge
        Dt(c[end]) ~ (-Aᵣ[end]*J.u/F/p.aₖ - Dₗ[end]*Aₗ[end]*(c[end] - c[end-1])/Δr)/Vᵢ[end]

    ]

    System(eqns,t; name=name,systems=[J])
end