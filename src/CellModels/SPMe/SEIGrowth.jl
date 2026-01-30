using ModelingToolkit

function SEIGrowth(; name, p::SideReactionParameters, g)
    
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
        (η_sei(t))[1:N]
        Lsei_x(t)
    end
    
    k_sei = 1
    Rsei = 200000 # SEI resistance factor
    Usei = 0.4 # SEI open circuit potential
    c_0 = 4541 # Initial ec concentration for SEI growth in mol*m^-3
    D_ec = 2e-18 # Diffusivity of ec in SEI in m^2*s^-1
    alpha = 0.5 # SEI transfer coefficient
    n_sei = 2 # Number of electrons transferred in SEI reaction
    ρ_sei = 1690 # Density of SEI in kg*m^-3
    M_sei = 0.162 # Molar mass of SEI in kg*mol^-1
    
    # Scott Marquis thesis (eq. 5.92)
    j_sei = -p.j_sei₀.*exp.(-p.α.*η_sei.*F/R/T.u)
    
    eqns = [
        [η_sei[i] ~ ϕₛ.u[i] - ϕₑ.u[i] - Usei - J.u.*L_sei[i].*Rsei for i in 1:N]...,
        [Dt(L_sei[i]) ~ -j_sei[i]* M_sei/(n_sei*F*ρ_sei) for i in 1:N]...,
        Lsei_x ~ sum([L_sei[i] for i in 1:N])/N
    ]

    System(eqns,t; name=name,systems=[J, T, ϕₑ, ϕₛ])
end