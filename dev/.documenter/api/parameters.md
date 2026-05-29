
# Parameters API Reference {#Parameters-API-Reference}

Complete API documentation for battery parameter types and constructors.

## Parameter Constructors {#Parameter-Constructors}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.Chen2020' href='#BatteryToolkit.Chen2020'><span class="jlbinding">BatteryToolkit.Chen2020</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Chen2020()
```


Create a battery parameter set based on Chen et al. 2020.

Returns a BatteryParameters object with parameters for a pouch-type lithium-ion battery cell using graphite negative electrode and NMC positive electrode. These parameters are calibrated from experimental data and are suitable for general lithium-ion battery simulation.

**Features**
- Graphite negative electrode with SEI side reaction
  
- NMC positive electrode without side reactions
  
- Typical electrolyte composition parameters
  
- Cell capacity: 5 Ah nominal
  

**Returns**
- `BatteryParameters` object ready for use with `SPMe()` cell models
  

**Example**

```julia
params = Chen2020()
sys = SPMe(params=params)
```


**References**

See Chen et al. 2020 in battery literature for detailed parameter derivation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/ParameterSets/Chen2020.jl#L70-L96" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.OKane2022' href='#BatteryToolkit.OKane2022'><span class="jlbinding">BatteryToolkit.OKane2022</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
OKane2022()
```


Create a battery parameter set based on O'Kane et al. 2022.

Returns a BatteryParameters object with parameters for a high-fidelity graphite/NMC lithium-ion battery cell. These parameters are derived from extensive electrochemical characterization and are suitable for detailed electrochemical simulations.

**Features**
- Graphite negative electrode with SEI side reaction
  
- NMC positive electrode without side reactions
  
- High-precision electrolyte parameters
  
- Consistent with PyBaMM (Python Battery Mathematical Modelling) standard parameters
  
- Cell capacity: 5 Ah nominal
  

**Returns**
- `BatteryParameters` object ready for use with `SPMe()` cell models
  

**Example**

```julia
params = OKane2022()
sys = SPMe(params=params, N=10)  # 10 FVM nodes per domain
```


**References**

See O'Kane et al. 2022 for detailed electrochemical characterization and model validation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/ParameterSets/OKane2022.jl#L70-L97" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Parameter Types {#Parameter-Types}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.BatteryParameters' href='#BatteryToolkit.BatteryParameters'><span class="jlbinding">BatteryToolkit.BatteryParameters</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BatteryParameters
```


Complete parameter set for a lithium-ion battery cell.

Aggregates all electrochemical and geometric parameters needed for SPMe simulations, including positive/negative electrode, electrolyte, and pack configuration.

**Fields**
- `p::SolidParticleParameters`: Positive electrode (cathode) parameters
  
- `n::SolidParticleParameters`: Negative electrode (anode) parameters
  
- `e::ElectrolyteParameters`: Electrolyte parameters
  
- `Hcc::Float64`: Current collector height (m)
  
- `Wcc::Float64`: Current collector width (m)
  
- `n_el::Int`: Number of parallel electrode pairs
  
- `Q₀::Float64`: Nominal cell capacity (Ah)
  
- `Vmin::Float64`: Minimum safe voltage (V)
  
- `Vmax::Float64`: Maximum safe voltage (V)
  

**Example**

```julia
params = Chen2020()  # Pre-configured parameter set
sys = SPMe(params=params)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/ParameterSets/Base.jl#L111-L135" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.ElectrolyteParameters' href='#BatteryToolkit.ElectrolyteParameters'><span class="jlbinding">BatteryToolkit.ElectrolyteParameters</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ElectrolyteParameters
```


Parameters describing lithium-ion battery electrolyte properties.

**Fields**
- `Lₚ::Float64`: Positive electrode thickness (m)
  
- `Lₛ::Float64`: Separator thickness (m)
  
- `Lₙ::Float64`: Negative electrode thickness (m)
  
- `Dₑ::Function`: Electrolyte diffusivity as function of concentration (m²/s)
  
- `σₑ::Function`: Electrolyte conductivity as function of concentration (S/m)
  
- `c₀::Float64`: Initial electrolyte concentration (mol/m³)
  
- `cₜ::Float64`: Typical electrolyte concentration (mol/m³)
  
- `t₊::Float64`: Transference number (cation fraction in current)
  
- `ϵₚ::Float64`: Positive electrode porosity
  
- `ϵₛ::Float64`: Separator porosity
  
- `ϵₙ::Float64`: Negative electrode porosity
  
- `bₚ::Float64`: Positive electrode Bruggeman coefficient
  
- `bₛ::Float64`: Separator Bruggeman coefficient
  
- `bₙ::Float64`: Negative electrode Bruggeman coefficient
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/ParameterSets/Base.jl#L2-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SideReactionParameters' href='#BatteryToolkit.SideReactionParameters'><span class="jlbinding">BatteryToolkit.SideReactionParameters</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SideReactionParameters
```


Parameters for secondary reactions (SEI growth, lithium plating) on electrode surfaces.

**Fields**
- `name::Symbol`: Reaction identifier (e.g., `:sei`)
  
- `k::Float64`: Reaction rate coefficient (A/m²)
  
- `α::Float64`: Transfer coefficient (charge transfer kinetics parameter)
  
- `M::Float64`: Molar mass of reaction product (kg/mol)
  
- `z::Int`: Number of electrons transferred per product molecule
  
- `ρ::Float64`: Product film density (kg/m³)
  
- `σ::Float64`: Product film conductivity (S/m)
  
- `U::Float64`: Open circuit potential of reaction (V)
  
- `Lf₀::Float64`: Initial film thickness (m)
  
- `V̄::Float64`: Partial molar volume of product (m³/mol)
  
- `R::Float64`: Film resistivity (Ω·m)
  
- `j_sei₀::Float64`: Exchange current density (A/m²)
  
- `c::Function`: Concentration dependence function
  
- `D_sol::Float64`: Solvent diffusivity in film (m²/s)
  
- `c_sol::Float64`: Solvent concentration (mol/m³)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/ParameterSets/Base.jl#L40-L61" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SolidParticleParameters' href='#BatteryToolkit.SolidParticleParameters'><span class="jlbinding">BatteryToolkit.SolidParticleParameters</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SolidParticleParameters
```


Parameters describing lithium-ion electrode (active material particle) properties.

**Fields**
- `Rₖ::Float64`: Particle radius (m)
  
- `aₖ::Float64`: Interfacial area per unit volume (m⁻¹)
  
- `Dₖ::Function`: Solid-state diffusivity as function of concentration (m²/s)
  
- `σₖ::Float64`: Electronic conductivity (S/m)
  
- `c₀::Float64`: Initial lithium concentration (mol/m³)
  
- `c₊::Float64`: Maximum lithium concentration (mol/m³)
  
- `Uₖ::Function`: Open circuit potential as function of stoichiometry (V)
  
- `mₖ::Float64`: Reaction rate constant (A/m²·(mol/m³)⁻¹·⁵)
  
- `L_sei₀::Float64`: Initial SEI thickness (m)
  
- `side_reactions::Vector{SideReactionParameters}`: Secondary reactions on this electrode (default: [])
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/ParameterSets/Base.jl#L80-L96" target="_blank" rel="noreferrer">source</a></Badge>

</details>

