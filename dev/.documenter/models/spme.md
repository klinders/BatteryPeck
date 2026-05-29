
# SPMe Cell Model {#SPMe-Cell-Model}

Comprehensive guide to the Single Particle Model with Electrolyte (SPMe) implementation.

## Overview {#Overview}

The SPMe is a 1D electrochemical model capturing the essential physics of lithium-ion battery operation:
- **Solid-phase transport**: Lithium diffusion within electrode particles
  
- **Electrolyte transport**: Ionic migration and diffusion through porous media
  
- **Interfacial kinetics**: Charge transfer at electrode/electrolyte interfaces
  
- **Potential fields**: Electronic and ionic potential distributions
  
- **Side reactions**: SEI growth and lithium plating
  

## Model Structure {#Model-Structure}

The SPMe consists of several coupled components:

```
SPMe
├── SolidParticle (Negative)
├── SolidParticle (Positive)
├── Electrolyte
├── Potentials
├── SEI (optional)
└── LithiumPlating (optional)
```


## Components {#Components}

### SolidParticle (Electrode Model) {#SolidParticle-Electrode-Model}

Models radial lithium diffusion within spherical electrode particles using finite volume method.

**Key equations:**
- Diffusion: ∂c/∂t = ∇·(D∇c)
  
- Surface boundary condition: flux from interfacial current
  
- Center boundary condition: zero flux by symmetry
  

**Output variables:**
- `c`: Concentration profile (mol/m³)
  
- `c_surf`: Surface concentration
  
- `z`: Stoichiometry = c_surf/c_max
  
- `U₀`: Open-circuit potential (V)
  

**Configuration:**

```julia
sys = SPMe(params=params, N=10)  # N nodes per domain
```


### Electrolyte (Ion Transport) {#Electrolyte-Ion-Transport}

Models lithium-ion transport through porous electrode and separator domains.

**Key physics:**
- Ionic diffusion with concentration gradients
  
- Ionic migration driven by potential gradients
  
- Bruggeman correlation for tortuous paths in pores
  
- Three domains: negative electrode, separator, positive electrode
  

**Output variables:**
- `cₑ`: Electrolyte concentration profile (mol/m³)
  
- `ϕₑ`: Electrolyte potential (V)
  

### Potentials (Electrochemistry) {#Potentials-Electrochemistry}

Computes electrode potentials and overpotentials based on current and concentrations.

**Key calculations:**
- Butler-Volmer kinetics for interfacial current
  
- Reaction overpotential from concentration
  
- Ohmic potential drops in solids and electrolyte
  

### SEI (Solid Electrolyte Interface) {#SEI-Solid-Electrolyte-Interface}

Optional model for SEI film formation on negative electrode.

**Available implementations:**
1. **NoSEI**: No film formation (fastest)
  
2. **ReactionLimitedSEI**: Fast film growth limited by reaction kinetics
  
3. **SolventDiffusionLimitedSEI**: Slow film growth limited by solvent diffusion
  

**Physics:**
- SEI film grows during discharge
  
- Causes capacity fade and impedance rise
  
- Reduces cycling efficiency
  

**Selection:**

```julia
sys = SPMe(params=params, side_reactions=:NoSEI)
sys = SPMe(params=params, side_reactions=:ReactionLimitedSEI)
sys = SPMe(params=params, side_reactions=:SolventDiffusionLimitedSEI)
```


### LithiumPlating (Plating/Stripping) {#LithiumPlating-Plating/Stripping}

Optional model for lithium metal deposition on negative electrode during abuse conditions.

**Available implementations:**
1. **NoPlating**: No plating (fastest)
  
2. **IrreversiblePlating**: All plated Li becomes dead/inactive
  
3. **PartiallyReversiblePlating**: Some Li can be stripped (more realistic)
  

**Physics:**
- Plating occurs at low potentials (overdischarge)
  
- Represents failure mode in abuse scenarios
  
- Dead lithium reduces active inventory
  

**Selection:**

```julia
sys = SPMe(params=params, side_reactions=:NoPlating)
sys = SPMe(params=params, side_reactions=:IrreversiblePlating)
sys = SPMe(params=params, side_reactions=:PartiallyReversiblePlating)
```


## Creating an SPMe System {#Creating-an-SPMe-System}

### Basic Configuration {#Basic-Configuration}

```julia
using BatteryToolkit

params = Chen2020()
sys = SPMe(params=params)
```


### High-Fidelity Configuration {#High-Fidelity-Configuration}

```julia
params = OKane2022()
sys = SPMe(
    params=params,
    N=20,                              # More nodes for accuracy
    side_reactions=:SolventDiffusionLimitedSEI  # Include SEI
)
```


### With Lithium Plating {#With-Lithium-Plating}

```julia
sys = SPMe(
    params=params,
    N=15,
    side_reactions=:PartiallyReversiblePlating
)
```


## Model Resolution {#Model-Resolution}

The parameter `N` controls spatial discretization in each domain (negative, separator, positive):

|   N |  Accuracy |     Speed |              Use Case |
| ---:| ---------:| ---------:| ---------------------:|
|   5 |       Low | Very fast |   Feasibility studies |
|  10 |    Medium |      Fast | General use (default) |
|  15 |      Good |  Moderate |  Research simulations |
| 20+ | Excellent |      Slow |    High-fidelity work |


```julia
sys = SPMe(params=params, N=20)  # Use 20 nodes per domain
```


## Governing Equations {#Governing-Equations}

### Solid Phase {#Solid-Phase}

```
∂c_s/∂t = ∇·(D_s ∇c_s)  in solid particles
-D_s ∇c_s · n = j_Li    at electrode surface (current)
-D_s ∇c_s · n = 0       at particle center (symmetry)
```


### Electrolyte Phase {#Electrolyte-Phase}

```
ε ∂c_e/∂t = ∇·(D_e ∇c_e) + (1-t₊)/F ∇·(i_e)
i_e = -σ_e ∇ϕ_e - (σ_e RT/F) ∇ln(c_e)
```


### Electrochemistry {#Electrochemistry}

```
i = i₀ [exp(αF/RT η) - exp(-(1-α)F/RT η)]  Butler-Volmer
η = ϕ_s - ϕ_e - U₀
```


## Physical Parameters {#Physical-Parameters}

Key parameters controlling SPMe behavior:

|        Parameter | Symbol |               Range |                        Effect |
| ----------------:| ------:| -------------------:| -----------------------------:|
|      Diffusivity |     Dₖ | 1e-15 to 1e-13 m²/s |      Controls rate capability |
|     Conductivity |     σₖ |     0.1 to 1000 S/m |         Controls ohmic losses |
|     Surface area |     aₖ |      1e5 to 1e6 m⁻¹ | Controls interfacial kinetics |
| Exchange current |     i₀ |   1e-7 to 1e-5 A/m² |        Controls overpotential |
|     Transference |     t₊ |          0.2 to 0.4 |      Controls ionic transport |


## Output Variables {#Output-Variables}

Standard variables available in simulation results:

```julia
sol = simulate(sys, exp, Rodas4())

# Access variables
voltage = sol[sys.Pin]           # Battery voltage (V)
current = sol[sys.Iin]           # Battery current (A)
power = voltage .* current        # Instantaneous power (W)

# Cell internal states (if defined)
soc = sol[sys.cell.soc]          # State of charge
soh = sol[sys.cell.soh]          # State of health (if degradation model)
```


## Limitations and Assumptions {#Limitations-and-Assumptions}

**Assumptions:**
- Isothermal operation (no heat generation/dissipation)
  
- Perfect mixing in pores
  
- Spherical electrode particles
  
- No mechanical stress effects
  
- 1D spatial variation (no radial variation in electrodes)
  

**Limitations:**
- Does not capture local temperature effects
  
- Simplified SEI growth kinetics
  
- No active material cracking
  
- Limited to slow to moderate rates (&lt; 10C typical)
  

## Advanced Usage {#Advanced-Usage}

### Custom Parameter Sets {#Custom-Parameter-Sets}

```julia
custom_params = BatteryParameters(
    p = SolidParticleParameters(...),
    n = SolidParticleParameters(...),
    e = ElectrolyteParameters(...),
    Hcc=0.065, Wcc=1.58, n_el=1, Q₀=5,
    Vmin=2.0, Vmax=4.2
)
sys = SPMe(params=custom_params)
```


### Hybrid Pack Configuration {#Hybrid-Pack-Configuration}

Combine SPMe with pack geometry:

```julia
sys = SingleCellPack(params=params, config=(96,3), Qcell=5)
sys = MultiCellPack(params=params, config=(12,3), Qcell=5)
```


## Next Steps {#Next-Steps}
- [Pack Models](pack-models.md): Series/parallel battery packs
  
- [Side Reactions](side-reactions.md): SEI and plating models
  
- [Examples](../guide/examples.md): Practical simulation scenarios
  
- [API Reference](../api/cellmodels.md): Complete SPMe API
  
