using ModelingToolkit

function LithiumPlating(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named η_sei = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named cₑ = RealInputArray(nin=N)

    # Time derivative operator
    Dt = Differential(t)

    scale = 1000 # c_typical
    α_plating = 0.65 # Li plating transfer coefficient
    α_stripping = 1 - α_plating
    k_plating = 1e-11

    @variables begin
        # Plating concentration
        (c_plating(t))[1:N] = 0
        (c_dead(t))[1:N] = 0
        (j_stripping(t))[1:N]
        (ϕf(t))[1:N]

        c_plating_x(t)
        c_dead_x(t)
        j_stripping_x(t)
        ϕf_x(t)
    end

    η_stripping = [Δϕₛ.u[i] + η_sei.u[i] for i in 1:N]
    η_plating = -η_stripping

    j0_stripping = F*k_plating .*c_plating
    j0_plating = F*k_plating .*cₑ.u
    
    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_stripping[i] ~ -j0_plating[i]*exp(α_plating*F/R/T.u*η_plating[i]) for i in 1:N]...,
        
        # Irreversable
        [Dt(c_dead[i]) ~ -aₖ.u*j_stripping[i]/F for i in 1:N]...,
        [Dt(c_plating[i]) ~ 0 for i in 1:N]...,
        c_plating_x ~ sum([c_plating[i] for i in 1:N])/N,
        c_dead_x ~ sum([c_dead[i] for i in 1:N])/N,
        j_stripping_x ~ sum([j_stripping[i] for i in 1:N])/N,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ,η_sei, aₖ, cₑ])
end