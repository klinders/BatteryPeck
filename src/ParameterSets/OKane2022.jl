

sei_parameters = SideReactionParameters(
    name = :sei,
    k = 1.0e-10, # Reaction rate
    α = 0.5, # Side reaction transfer coefficient
    M = 0.162, # Molar mass of SR product
    z = 1, # Ratio of Li to SEI moles
    ρ = 1690, # Density of SR product
    σ = 5e-6, #8.95e-14 Conductivity in the SEI layer
    U = 0.4, # Open circuit potential of SR
    Lf₀ = 5e-9, # Initial thickness of SR film
    V̄ = 9.585e-05, # Partial molar volume
    R = 2e5, # SEI resistivity
    j_sei₀ = 1.5e-7, # SEI Reaction exchangcurrent [A.m-2]
    c = (el) -> nothing, # Concentration dependence function
    D_sol = 2.5e-22, # Solvent diffusivity in m^2/s
    c_sol = 2636.0, # Solvent concentration in mol/m^3
)

n = SolidParticleParameters(
    Rₖ = 5.86e-6, # Radius of the electrode in m
    aₖ = 383960, # Surface area density in m^-1
    Dₖ = c -> 3.3e-14, # electrode diffusivity in m^2*s^-1
    σₖ = 215, # Conductivity in S*m^-1
    c₀ = 29866, # Initial electrode concentration in mol*m^-3
    c₊ = 33133, # Maximum electrode concentration in mol*m^-3
    # Open-circuit potential in V
    Uₖ = z->1.9793*exp(-39.3631*z) + 0.2482-0.0909*tanh(29.8538*(z-0.1234)) - 0.04478*tanh(14.9159*(z-0.2769)) - 0.0205*tanh(30.4444*(z-0.6103)), 
    mₖ = 6.48e-7, # Reaction rate constant in A*m^-2*(mol*m^-3)^-1.5
    L_sei₀ = 5e-9,
    side_reactions = [sei_parameters]
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
    mₖ = 3.42e-6, # Reaction rate constant in A*m^-2*(mol*m^-3)^-1.5
    L_sei₀ = 0,

)

e = ElectrolyteParameters(
    Lₚ= 75.6e-6, # Length of the electrode in m
    Lₛ= 12e-6, # Length of the separator in m
    Lₙ= 85.2e-6, # Length of the negative electrode in m
    Dₑ= c -> 8.794e−11 * (c/1000)^2 − 3.972e−10 * (c/1000) + 4.862e−10, # Diffusivity in m^2*s^-1
    σₑ= c -> 0.1297 * (c/1000)^3 − 2.51*(c/1000)^1.5 + 3.329*(c/1000), # Conductivity in S*m^-1
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


"""
    OKane2022()

Create a battery parameter set based on O'Kane et al. 2022.

Returns a BatteryParameters object with parameters for a high-fidelity graphite/NMC
lithium-ion battery cell. These parameters are derived from extensive electrochemical
characterization and are suitable for detailed electrochemical simulations.

# Features
- Graphite negative electrode with SEI side reaction
- NMC positive electrode without side reactions
- High-precision electrolyte parameters
- Consistent with PyBaMM (Python Battery Mathematical Modelling) standard parameters
- Cell capacity: 5 Ah nominal

# Returns
- `BatteryParameters` object ready for use with `SPMe()` cell models

# Example
```julia
params = OKane2022()
sys = SPMe(params=params, N=10)  # 10 FVM nodes per domain
```

# References
See O'Kane et al. 2022 for detailed electrochemical characterization and model validation.
"""
function OKane2022()
    return BatteryParameters(
        p = p,
        n = n,
        e = e,
        Hcc = 0.065,
        Wcc = 1.58,
        n_el = 1,
        Q₀ = 5, # Original battery capacity in Ah
        Vmin = 2.0,
        Vmax = 4.2
    )
end
