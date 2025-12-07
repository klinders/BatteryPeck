using ModelingToolkit

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
    
    @variables begin
        # I am adding two ghost nodes for the boundary conditions
        (c(t))[1:g.Nᵣ] = repeat([p.c₀],g.Nᵣ)
        c_avr(t)
        c_surf(t)
        U₀(t)
        z(t)
    end

    # Discretized equations
    Δr,r,Vᵢ,Aₗ,Aᵣ = g.Δr, g.r, g.Vᵢ, g.Aₗ, g.Aᵣ
    Dₗ = [nothing, [D_face(p.Dₖ(c[i-1]), p.Dₖ(c[i]),Δr,Δr) for i in 2:g.Nᵣ]...] # Left diffusivities
    Dᵣ = [[D_face(p.Dₖ(c[i]), p.Dₖ(c[i+1]),Δr,Δr) for i in 1:g.Nᵣ-1]..., nothing] # Right diffusivities

    eqns = [
        c_avr ~ sum(c)/g.Nᵣ
        c_surf ~ 1.5*c[end] - 0.5*c[end-1] # Surface concentration
        z ~ c_surf/p.c₊ # Stoichiometry
        U₀ ~ p.Uₖ(z) # Open-circuit potential

        # Boundary condition center
        Dt(c[1]) ~ (Dᵣ[1]*Aᵣ[1]*(c[2] - c[1])/Δr)/Vᵢ[1]

        # Internal nodes
        [Dt(c[i]) ~ (Dᵣ[i]*Aᵣ[i]*(c[i+1] - c[i])/Δr - Dₗ[i]*Aₗ[i]*(c[i] - c[i-1])/Δr)/Vᵢ[i] for i in 2:g.Nᵣ-1]...

        # Boundary condition edge
        Dt(c[end]) ~ (-Aᵣ[end]*J.u/F/p.aₖ - Dₗ[end]*Aₗ[end]*(c[end] - c[end-1])/Δr)/Vᵢ[end]

    ]

    System(eqns,t; name=name,systems=[J, T])
end