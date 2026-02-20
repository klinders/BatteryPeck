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
    @named ϕₑ = RealInputArray(nin=N)
    @named ϕₛ = RealInputArray(nin=N)

    # Time derivative operator
    Dt = Differential(t)
    
    Lsei_0 = 5e-09 # Initial SEI thickness in m
    
    @variables begin
        # I am adding two ghost nodes for the boundary conditions
        (L_sei(t))[1:N] = Lsei_0
        # (η_sei(t))[1:N]
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]

        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
    end
    
    k_sei = 1e-12 # SEI growth rate constant in m/s
    Usei = 0.4 # SEI open circuit potential
    c_0 = 4541 # Initial ec concentration for SEI growth in mol*m^-3
    D_ec = 2e-19 # Diffusivity of ec in SEI in m^2*s^-1
    alpha = 0.5 # SEI transfer coefficient
    n_sei = 2 # Number of electrons transferred in SEI reaction
    ρ_sei = 1690 # Density of SEI in kg*m^-3
    M_sei = 0.162 # Molar mass of SEI in kg*mol^-1
    σ_sei = 5e-6 # SEI conductivity in S*m^-1
    Rsei = 2e5
    # Rsei = 1/s.aₖ/σ_sei # SEI resistivity

    η_sei = [ϕₛ.u[i] - ϕₑ.u[i] - Usei - ϕf[i] for i in 1:N]
    
    k_exp = [s.aₖ*k_sei*exp(-p.α*η_sei[i]*F/R/T.u) for i in 1:N]
    LoverD = [L_sei[i]/D_ec for i in 1:N]
    j_sei_f = [-F*c_0*k_exp[i]/(1+LoverD[i]*k_exp[i]) for i in 1:N]
    c_ec = [c_0/(1 + LoverD[i]*k_exp[i]) for i in 1:N]

    c_ec_av = sum(c_ec)/N

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        [Dt(L_sei[i]) ~ -j_sei[i]*M_sei/(n_sei*F*ρ_sei) for i in 1:N]...,
        [j_sei[i] ~ j_sei_f[i] for i in 1:N]...,
        [ϕf[i] ~ -J.u*L_sei[i]*Rsei/s.aₖ for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N
    ]

    System(eqns,t; name=name,systems=[J, T, ϕₑ, ϕₛ])
end