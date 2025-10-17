# include("Base.jl")

n = SolidParticleParameters(
    Rₖ = 5.86e-6, # Radius of the electrode in m
    aₖ = 3.84e5, # Surface area density in m^-1
    Dₖ = c -> 3.3e-14, # electrode diffusivity in m^2*s^-1
    σₖ = 215, # Conductivity in S*m^-1
    c₀ = 29866, # Initial electrode concentration in mol*m^-3
    c₊ = 33133, # Maximum electrode concentration in mol*m^-3
    # Open-circuit potential in V
    Uₖ = z->1.9793*exp(-39.3631*z) + 0.2482-0.0909*tanh(29.8538*(z-0.1234)) - 0.04478*tanh(14.9159*(z-0.2769)) - 0.0205*tanh(30.4444*(z-0.6103)), 
    mₖ = 6.48e-7, # Reaction rate constant in A*m^-2*(mol*m^-3)^-1.5
    
    # SEI parameters
    k_sei = 1e-12, # SEI reaction rate constant in m/s
    c_sei₀ = 4541, # Initial SEI concentration in mol/m^3
    U_sei = 0.0, # SEI equilibrium potential in V
    D_sei = 2e-19, # SEI diffusivity in m^2/s
    M_sei = 0.162, # SEI molar mass in kg/mol
    ρ_sei = 1690, # SEI film density in kg/m^-3
    n_sei = 2, # SEI reaction order
    σ_sei = 5e-6, # SEI film conductivity in S/m
    L_sei₀ = 5e-9, # Initial SEI film thickness in
)

p = SolidParticleParameters(
    Rₖ = 5.22e-6, # Radius of the electrode in m
    aₖ = 3.82e5, # Surface area density in m^-1
    Dₖ = c -> 4.0e-15, # electrode diffusivity in m^2*s^-1
    σₖ = 0.18, # Conductivity in S*m^-1
    c₀ = 17038, # Initial electrode concentration in mol*m^-3
    c₊ = 63104, # Maximum electrode concentration in mol*m^-3
    # Open-circuit potential in V
    Uₖ = z->-0.8090*z + 4.4875 - 0.0428*tanh(18.5138*(z-0.5542)) - 17.7326*tanh(15.7890*(z-0.3117)) + 17.5842*tanh(15.9308*(z-0.3120)), 
    mₖ = 3.42e-6 # Reaction rate constant in A*m^-2*(mol*m^-3)^-1.5
)

e = ElectrolyteParameters(
    Lₚ= 75.6e-6, # Length of the electrode in m
    Lₛ= 12e-6, # Length of the separator in m
    Lₙ= 85.2e-6, # Length of the negative electrode in m
    Dₑ= c -> 8.794e−17 * c^2 − 3.972e−13 * c + 4.862e−10, # Diffusivity in m^2*s^-1
    σₑ= c -> 1.297e−10 * c^3 − 7.937e−5  * c^1.5 + 3.329e-3*c, # Conductivity in S*m^-1
    c₀= 1000.0, # Initial electrolyte concentration in mol*m^-3
    cₜ= 1000.0, # Typical Electrolyte concentration in mol*m^-3
    t₊= c->0.2594, # transer number
    ϵₚ= 0.335, # Porosity of the positive electrode
    ϵₛ= 0.47, # Porosity of the separator
    ϵₙ= 0.25, # Porosity of the negative electrode
    # bₚ = 2.43, # Bruggeman coefficient
    # bₛ = 2.57,
    # bₙ = 2.91,
    bₚ = 1.5, # Bruggeman coefficient
    bₛ = 1.5,
    bₙ = 1.5,
)

function Chen2020()
    return BatteryParameters(
        p = p,
        n = n,
        e = e,
        i₀ = 48.69, # Typical current density for 1C in A/m^-2
        Q₀ = 5, # Original battery capacity in Ah
        Vmin = 2.5,
        Vmax = 5
    )
end
