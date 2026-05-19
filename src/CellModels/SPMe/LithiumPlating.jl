module LithiumPlating

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

function NoPlating(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
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

    k_plating = 1e-9
    j_strip0 = F*k_plating*1000

    @variables begin
        # Plating concentration
        (c_plating(t))[1:N] = 0
        (c_dead(t))[1:N] = 0
        (j_stripping(t))[1:N], [guess=j_strip0]
        (ϕf(t))[1:N]

        c_plating_x(t)
        c_dead_x(t)
        j_stripping_x(t)
        ϕf_x(t)
    end

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_stripping[i] ~ 0 for i in 1:N]...,
        
        # Irreversable
        [Dt(c_dead[i]) ~ 0 for i in 1:N]...,
        [Dt(c_plating[i]) ~ 0 for i in 1:N]...,
        c_plating_x ~ sum([c_plating[i] for i in 1:N])/N,
        c_dead_x ~ sum([c_dead[i] for i in 1:N])/N,
        j_stripping_x ~ sum([j_stripping[i] for i in 1:N])/N,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ,η_sei, aₖ, cₑ])
end

function IrreversiblePlating(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
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
    k_plating = 1e-9
    V̄ = 1.3e-05 # Partial molar volume of lithium [m3.mol-1]

    j_strip0 = F*k_plating*1000

    @variables begin
        # Plating concentration
        (c_plating(t))[1:N] = 0
        (c_dead(t))[1:N] = 0
        (L_plating(t))[1:N] = 0
        (L_dead(t))[1:N] = 0
        (j_stripping(t))[1:N], [guess=j_strip0]
        (ϕf(t))[1:N]
        (η_plating(t))[1:N]
        (η_stripping(t))[1:N]
        (j0_plating(t))[1:N]
        (j0_stripping(t))[1:N]

        c_plating_x(t)
        c_dead_x(t)
        L_plating_x(t)
        L_dead_x(t)
        j_stripping_x(t)
        ϕf_x(t)
        η_plating_x(t)
        η_stripping_x(t)
    end
    
    eqns = [
        [j0_stripping[i] ~ F*k_plating*c_plating[i] for i in 1:N]...,
        [j0_plating[i] ~ F*k_plating*cₑ.u[i] for i in 1:N]...,
        [η_stripping[i] ~ Δϕₛ.u[i] + η_sei.u[i] for i in 1:N]...,
        [η_plating[i] ~ -η_stripping[i] for i in 1:N]...,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_stripping[i] ~ -j0_plating[i]*exp(α_plating*F/R/T.u*0.85*η_plating[i]) for i in 1:N]...,
        
        # Irreversable
        [Dt(c_dead[i]) ~ -aₖ.u*j_stripping[i]/F for i in 1:N]...,
        [Dt(c_plating[i]) ~ 0 for i in 1:N]...,
        c_plating_x ~ sum([c_plating[i] for i in 1:N])/N,
        c_dead_x ~ sum([c_dead[i] for i in 1:N])/N,
        j_stripping_x ~ sum([j_stripping[i] for i in 1:N])/N,
        [L_plating[i] ~ c_plating[i]*V̄/aₖ.u for i in 1:N]...,
        [L_dead[i] ~ c_dead[i]*V̄/aₖ.u for i in 1:N]...,
        L_plating_x ~ sum([L_plating[i] for i in 1:N])/N,
        L_dead_x ~ sum([L_dead[i] for i in 1:N])/N,
        η_plating_x ~ sum([η_plating[i] for i in 1:N])/N,
        η_stripping_x ~ sum([η_stripping[i] for i in 1:N])/N,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ,η_sei, aₖ, cₑ])
end

end