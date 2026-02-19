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
    c_0 = 4541 # Initial ec concentration for SEI growth in mol*m^-3
    
    @variables begin
        # I am adding two ghost nodes for the boundary conditions
        (c_sei(t))[1:N] = c_0
        c_sei_av(t)
        (L_sei(t))[1:N]
        L_sei_av(t)
        (j_sei(t))[1:N]
        j_sei_av(t)
        (ϕf(t))[1:N]
        ϕf_av(t)
    end
    
    k_sei = 1e-12 # SEI growth rate constant in m/s
    Usei = 0.4 # SEI open circuit potential
    D_ec = 2e-19 # Diffusivity of ec in SEI in m^2*s^-1
    alpha = 0.5 # SEI transfer coefficient
    n_sei = 2 # Number of electrons transferred in SEI reaction
    ρ_sei = 1690 # Density of SEI in kg*m^-3
    M_sei = 0.162 # Molar mass of SEI in kg*mol^-1
    σ_sei = 5e-6 # SEI conductivity in S*m^-1
    Rsei = 2e5
    V̄sei = 9.585e-05
    # Rsei = 1/s.aₖ/σ_sei # SEI resistivity

    η_sei = [ϕₛ.u[i] - ϕₑ.u[i] - Usei - ϕf[i] for i in 1:N]
    
    k_exp = [k_sei*exp(-p.α*η_sei[i]*F/R/T.u) for i in 1:N]
    LoverD = [L_sei[i]/D_ec for i in 1:N]
    c_ec = [c_0/(1 + LoverD[i]*k_exp[i]) for i in 1:N]

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        [Dt(c_sei[i]) ~ 0 for i in 1:N],#-j_sei[i]*s.aₖ/(n_sei*F) for i in 1:N]...,
        c_sei_av ~ sum(c_sei)/N,
        [L_sei[i] ~ c_sei[i]/V̄sei for i in 1:N]...,
        L_sei_av ~ sum(L_sei)/N,
        [j_sei[i] ~ -F*c_0*k_exp[i]/(1+LoverD[i]*k_exp[i]) for i in 1:N]...,
        j_sei_av ~ sum(j_sei)/N,
        [ϕf[i] ~ J.u*L_sei[i]*Rsei/s.aₖ for i in 1:N]...,
        ϕf_av ~ sum(ϕf)/N
    ]

    System(eqns,t; name=name,systems=[J, T, ϕₑ, ϕₛ])
end