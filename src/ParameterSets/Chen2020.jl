# Source [1]: https://doi.org/10.1149/1945-7111/ab9050        # Chen2020 dataset
# Source [2]: https://doi.org/10.1016/j.electacta.2022.140700 # Regan2022 stoichiometry vs. SoC

# Import base parameters
# include("Base.jl")

# Import packages
import NaNMath # Return NaN for log or sqrt of negative one

sei_parameters = SideReactionParameters(
    name = :sei,
    k = 1.0e-10, # Reaction rate
    α = 0.5, # Side reaction transfer coefficient
    M = 0.162, # Molar mass of side reaction product
    z = 2, # Ratio of lithium to side reaction moles
    ρ = 1690, # Density of side reaction product
    σ = 5e-6, # Conductivity in the side reaction layer
    U = 0.4, # Open circuit potential of side reaction
    Lf₀ = 5e-9, # Initial thickness of side reaction film
    V̄ = 9.585e-05, # Partial molar volume
    R = 2e5, # Side reaction resistivity
    j_sei₀ = 1.5e-7, # Side reaction exchange current
    c = (el) -> nothing # Concentration dependence function
)

n = SolidParticleParameters(
    Rₖ = 5.86e-6, # Radius of the electrode in m
    aₖ = 3.84e5, # Surface area density in m^-1
    Dₖ = c -> 3.3e-14, # Electrode diffusivity in m^2*s^-1
    σₖ = 215, # Conductivity in S*m^-1
    c₀ = 29866, # Initial electrode concentration in mol*m^-3
    c₊ = 33133, # Maximum electrode concentration in mol*m^-3
    # Open-circuit potential in V
    Uₖ = z->1.9793*exp(-39.3631*z) + 0.2482-0.0909*tanh(29.8538*(z-0.1234)) - 0.04478*tanh(14.9159*(z-0.2769)) - 0.0205*tanh(30.4444*(z-0.6103)), 
    mₖ = 6.48e-7, # Reaction rate constant in A*m^-2*(mol*m^-3)^-1.5
    L_sei₀ = 5e-9,
    side_reactions = [sei_parameters],
    # Stoichiometry versus state of charge
    z_0 = 0.0279,   # Stoichiometry at zero percent state of charge
    z_100 = 0.9014, # Stoichiometry at full state of charge
)

p = SolidParticleParameters(
    Rₖ = 5.22e-6, # Radius of the electrode in m
    aₖ = 3.82e5, # Surface area density in m^-1
    Dₖ = c -> 4.0e-15, # Electrode diffusivity in m^2*s^-1
    σₖ = 0.18, # Conductivity in S*m^-1
    c₀ = 17038, # Initial electrode concentration in mol*m^-3
    c₊ = 63104, # Maximum electrode concentration in mol*m^-3
    # Open-circuit potential in V
    Uₖ = z->-0.8090*z + 4.4875 - 0.0428*tanh(18.5138*(z-0.5542)) - 17.7326*tanh(15.7890*(z-0.3117)) + 17.5842*tanh(15.9308*(z-0.3120)), 
    mₖ = 3.42e-6, # Reaction rate constant in A*m^-2*(mol*m^-3)^-1.5
    L_sei₀ = 0,
    # Stoichiometry versus state of charge
    z_0 = 0.9072,   # Stoichiometry at zero percent state of charge
    z_100 = 0.2567, # Stoichiometry at full state of charge
)

e = ElectrolyteParameters(
    Lₚ= 75.6e-6, # Length of the electrode in m
    Lₛ= 12e-6, # Length of the separator in m
    Lₙ= 85.2e-6, # Length of the negative electrode in m
    Dₑ= c -> 8.794e−11 * (c/1000)^2 − 3.972e−10 * (c/1000) + 4.862e−10, # Diffusivity in m^2*s^-1
    # Safely handle negative concentrations using NaNMath to prevent solver crashing
    σₑ= c -> 0.1297 * (c/1000)^3 − 2.51 * NaNMath.pow(c/1000, 1.5) + 3.329*(c/1000), # Conductivity in S*m^-1
    c₀= 1000.0, # Initial electrolyte concentration in mol*m^-3
    cₜ= 1000.0, # Typical electrolyte concentration in mol*m^-3
    t₊= c->0.2594, # Transfer number
    ϵₚ= 0.335, # Porosity of the positive electrode
    ϵₛ= 0.47, # Porosity of the separator
    ϵₙ= 0.25, # Porosity of the negative electrode
    bₚ = 1.5, # Bruggeman coefficient
    bₛ = 1.5,
    bₙ = 1.5,
)

function Chen2020()
    return BatteryParameters(
        p = p,
        n = n,
        e = e,
        i₀ = 48.69, # Typical current density for single C rate
        Hcc = 0.065,
        Wcc = 1.58,
        n_el = 1,
        Q₀ = 5, # Original battery capacity in Ah
        Vmin = 2.5,
        Vmax = 4.2
    )
end