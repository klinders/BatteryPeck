using ModelingToolkit

function SEIGrowth(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = 0#p.Lf₀/p.V̄*s.aₖ
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
    end

    η_sei = [Δϕₛ.u[i] - p.U - ϕf[i] for i in 1:N]
    
    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ -p.j_sei₀*exp(-p.α*F/R/T.u*η_sei[i]) for i in 1:N]...,

        [Dt(c_sei[i]) ~ 0 for i in 1:N]..., #s.aₖ*j_sei[i]/(F/p.z)
        [L_sei[i] ~ c_sei[i]*p.V̄/s.aₖ for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]/s.aₖ/p.σ for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ])
end