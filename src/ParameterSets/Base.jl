
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

Base.@kwdef mutable struct SideReactionParameters
    name            # Name of the side reaction
    k               # Reaction rate
    α               # Side reaction transfer coefficient
    M               # Molar mass of SR product
    z               # Ratio of Li to SEI moles
    ρ               # Density of SR product
    σ               # Conductivity in the SEI layer
    U               # Open circuit potential of SR
    Lf₀             # Initial thickness of SR film
    V̄               # SEI Partial molar volume
    R               # SEI Resistivity
    j_sei₀          # Reaction exchange current
    c::Function     # Concentration dependence function
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
    L_sei₀ # Initial thickness of SEI film

    side_reactions::Vector{SideReactionParameters} = SideReactionParameters[]
end

Base.@kwdef mutable struct BatteryParameters
    p::SolidParticleParameters # Parameters for the positive electrode
    n::SolidParticleParameters # Parameters for the negative electrode
    e::ElectrolyteParameters # Parameters for the electrolyte
    
    Hcc # Current collector height
    Wcc # Current collector width
    n_el # Number of parallel electrodes
    Q₀ # Original battery capacity in Ah

    Vmin
    Vmax
end