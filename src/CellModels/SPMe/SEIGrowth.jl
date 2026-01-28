using ModelingToolkit

function SEIGrowth(; name, p::SideReactionParameters)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant

    @named J = RealInput()
    @named T = RealInput()
    @named ϕₑ = RealInput()
    @named ϕₛ = RealInput()

    # Time derivative operator
    Dt = Differential(t)
    
    Lsei_0 = 5e-09 # Initial SEI thickness in m
    
    @variables begin
        # I am adding two ghost nodes for the boundary conditions
        L_sei(t) = Lsei_0
        j_sei(t)
        c_ec(t)
        η_sei(t)
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

    k_exp = k_sei * exp(-alpha * F/(R*T.u) * η_sei)
    
    eqns = [
        η_sei ~ ϕₛ.u - sum(ϕₑ.u)/length(ϕₑ.u) - Usei - J.u*L_sei*Rsei
        j_sei ~ -F*c_0*k_exp/(1 + k_exp*L_sei/D_ec)
        c_ec ~ c_0/(1 + k_exp*L_sei/D_ec)
        Dt(L_sei) ~ -j_sei* M_sei/(n_sei*F*ρ_sei)
    ]

    System(eqns,t; name=name,systems=[J, T, ϕₑ, ϕₛ])
end