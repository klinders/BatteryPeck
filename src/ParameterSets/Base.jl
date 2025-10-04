
Base.@kwdef mutable struct ElectrolyteParameters
    Lₚ # Length of the electrode in m
    Lₛ # Length of the separator in m
    Lₙ # Length of the negative electrode in m
    Dₑ # Diffusivity in m^2*s^-1
    σₑ # Conductivity in S*m^-1
    c₀ # Initial electrolyte concentration in mol*m^-3
    cₜ # Typical Electrolyte concentration in mol*m^-3
    t₊ # transer number
    ϵₚ # Porosity of the positive electrode
    ϵₛ # Porosity of the separator
    ϵₙ # Porosity of the negative electrode
    bₙ # Negative electrode Bruggeman coefficient
    bₛ # Separator Bruggeman coefficient
    bₚ # Positive electrode Bruggeman coefficient
end

Base.@kwdef mutable struct SolidParticleParameters
    Rₖ # Radius of the electrode in m
    aₖ # Surface area density in m^-1
    Dₖ # electrode diffusivity in m^2*s^-1
    σₖ # Conductivity in S*m^-1
    c₀ # Initial electrode concentration in mol*m^-3
    c₊ # Maximum electrode concentration in mol*m^-3
    Uₖ # Open-circuit potential in V
    mₖ # Reaction rate constant in A*m^-2*(mol*m^-3)^-1.5

    # SEI parameters
    k_sei  = nothing # SEI reaction rate constant in m/s
    c_sei₀ = nothing # Initial SEI concentration in mol/m^3
    U_sei  = nothing # SEI equilibrium potential in V
    D_sei  = nothing # SEI diffusivity in m^2/s
    M_sei  = nothing # SEI molar mass in kg/mol
    ρ_sei  = nothing # SEI film density in kg/m^-3
    n_sei  = nothing # SEI reaction order
    σ_sei  = nothing # SEI film conductivity in S/m
    L_sei₀ = nothing # Initial SEI film thickness in
end

Base.@kwdef mutable struct BatteryParameters
    p::SolidParticleParameters # Parameters for the positive electrode
    n::SolidParticleParameters # Parameters for the negative electrode
    e::ElectrolyteParameters # Parameters for the electrolyte
    
    i₀ # Typical current density for 1C in A/m^-2
    Q₀ # Original battery capacity in Ah
end