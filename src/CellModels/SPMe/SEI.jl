module SEI

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

function NoSEI(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = 0
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N] = 0

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
    end

    η_sei = [Δϕₛ.u[i] - p.U - ϕf[i] for i in 1:N]
    
    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ 0 for i in 1:N]...,

        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        Q_sei ~ 0,

    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
end


function ReactionLimitedSEI(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = 0#scale
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_sei(t)
    end

    η_sei = [Δϕₛ.u[i] - p.U - ϕf[i] for i in 1:N]
    c_sei₀ = p.Lf₀/p.V̄*s.aₖ

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ -p.j_sei₀*exp(-p.α*F/R/T.u*η_sei[i]) for i in 1:N]...,

        [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        Q_sei ~ (c_sei_x-c_sei₀)*p.V̄*p.z*F/3600,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
end

function SolventDiffusionLimitedSEI(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)

    c_sei₀ = p.Lf₀/p.V̄*s.aₖ
    Hcc = 0.065
    Wcc = 1.58
    
    Vk = Hcc*Wcc*g.el.Ls[1] # Volume of SEI 

    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = c_sei₀
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_loss(t)
    end
    
    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ -p.D_sol*p.c_sol*F/L_sei[i] for i in 1:N]...,

        [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        Q_loss ~ (c_sei_x-c_sei₀)*Vk*p.z*F/3600,

    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
end

end