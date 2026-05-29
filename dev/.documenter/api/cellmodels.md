
# Cell Models API Reference {#Cell-Models-API-Reference}

API documentation for SPMe and related electrochemical models.

## Main Cell Model {#Main-Cell-Model}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SPMe' href='#BatteryToolkit.SPMe'><span class="jlbinding">BatteryToolkit.SPMe</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



The SPMe model implemented based on [Marquis _et al._ [1]](/references#MarquisEtAl2019) and [Brosa Planella and Widanage [2]](/references#BrosaPlanellaWidanage2023)

**Arguments**
- `name` (optional) defaults to SPMe
  
- `params ::BatteryParameters` Parameters from the given parameter set
  
- `Q ::Real` (optional) The capacity of the cell in Ah
  
- `N ::Dict(:Nₓ=>[::Int, ::Int, ::Int], :Nᵣ=>[::Int, ::Int])` (optional) number of mesh nodes in the particles and electrolyte
  
- `side_reactions ::Bool` (optional) enable side reactions
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/SPMe.jl#L20-L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## SPMe Components {#SPMe-Components}

### Solid Particle Model {#Solid-Particle-Model}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SolidParticle' href='#BatteryToolkit.SolidParticle'><span class="jlbinding">BatteryToolkit.SolidParticle</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
SolidParticle(; name, p::SolidParticleParameters, g)
```


Create a ModelingToolkit system for solid-state lithium diffusion in a battery electrode.

Models radial diffusion of lithium ions within spherical electrode particles using the finite volume method. Computes surface concentration, stoichiometry, and open-circuit potential.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::SolidParticleParameters`: Electrode material parameters
  
- `g`: FVM geometry object with node locations and volumes
  

**Input Ports**
- `J`: Surface current density (A/m²)
  
- `T`: Temperature (K)
  

**Output Variables**
- `c`: Concentration profile across particles (mol/m³)
  
- `c_surf`: Surface lithium concentration (mol/m³)
  
- `z`: Stoichiometry = c_surf/c_max
  
- `U₀`: Open-circuit potential (V)
  

**Notes**

Uses second-order accurate finite volume discretization with ghost nodes for boundary conditions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/SolidParticle.jl#L3-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Electrolyte Model {#Electrolyte-Model}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.Electrolyte' href='#BatteryToolkit.Electrolyte'><span class="jlbinding">BatteryToolkit.Electrolyte</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Electrolyte(; name, p::ElectrolyteParameters, g)
```


Create a ModelingToolkit system for electrolyte salt concentration and ionic transport.

Models lithium-ion transport through the electrolyte including diffusion in negative electrode, separator, and positive electrode. Computes concentration profiles and electrochemical potentials.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::ElectrolyteParameters`: Electrolyte material and transport parameters
  
- `g`: FVM geometry object defining domain and node locations
  

**Input Ports**
- `i_app`: Applied current density (A/m²)
  
- `T`: Temperature (K)
  
- `Δϕₙ`: Potential drop in negative electrode (V)
  
- `ϕₛn`: Solid potential in negative electrode (V)
  

**Output Variables**
- `cₑ`: Electrolyte concentration profile (mol/m³)
  
- `c̄ₑ`: Average electrolyte concentration (mol/m³)
  
- `ϕₑ`: Electrolyte potential profile (V)
  

**Notes**

Uses finite volume method with Bruggeman correlation for tortuosity in porous media. Automatically computes diffusion and migration based on concentration gradients.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/Electrolyte.jl#L3-L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Side Reactions {#Side-Reactions}

### SEI (Solid Electrolyte Interface) {#SEI-Solid-Electrolyte-Interface}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SEI.NoSEI' href='#BatteryToolkit.SEI.NoSEI'><span class="jlbinding">BatteryToolkit.SEI.NoSEI</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
NoSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
```


Create a zero SEI growth model (reaction disabled).

Returns a ModelingToolkit system where SEI film thickness remains zero and provides no ohmic resistance. Use this when SEI growth effects are negligible or you want to exclude them from the simulation.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::SideReactionParameters`: SEI reaction parameters (unused in this model)
  
- `s::SolidParticleParameters`: Electrode solid particle parameters
  
- `g`: FVM geometry object
  

**Output Variables**
- `L_sei`: SEI film thickness (always 0)
  
- `j_sei`: SEI current density (always 0)
  
- `ϕf`: Film potential (always 0)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/SEI.jl#L8-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SEI.ReactionLimitedSEI' href='#BatteryToolkit.SEI.ReactionLimitedSEI'><span class="jlbinding">BatteryToolkit.SEI.ReactionLimitedSEI</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
ReactionLimitedSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
```


Create a reaction-limited SEI growth model.

Models SEI film formation with reaction kinetics controlled by surface overpotential. The SEI current density follows Butler-Volmer kinetics. Use when SEI growth is fast (high overpotential) and film diffusion resistance is negligible.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::SideReactionParameters`: SEI reaction kinetic parameters
  
- `s::SolidParticleParameters`: Electrode solid particle parameters
  
- `g`: FVM geometry object
  

**Key Parameters Used**
- `p.j_sei₀`: Exchange current density (A/m²)
  
- `p.α`: Transfer coefficient (charge transfer kinetics)
  
- `p.U`: SEI formation potential (V vs Li/Li⁺)
  

**Output Variables**
- `L_sei`: SEI film thickness (grows over time)
  
- `j_sei`: SEI current density (determined by kinetics)
  
- `ϕf`: Film potential (typically small)
  

**Physical Assumption**

Reaction rate dominates over diffusion; film acts as perfect ionic conductor.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/SEI.jl#L81-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SEI.SolventDiffusionLimitedSEI' href='#BatteryToolkit.SEI.SolventDiffusionLimitedSEI'><span class="jlbinding">BatteryToolkit.SEI.SolventDiffusionLimitedSEI</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
SolventDiffusionLimitedSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
```


Create a solvent-diffusion-limited SEI growth model.

Models SEI film formation limited by solvent diffusion through the growing film. The SEI current density decreases as the film thickens due to increasing ionic resistance and decreasing solvent diffusion. Use when film resistance dominates over reaction kinetics.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::SideReactionParameters`: SEI reaction and film transport parameters
  
- `s::SolidParticleParameters`: Electrode solid particle parameters
  
- `g`: FVM geometry object
  

**Key Parameters Used**
- `p.D_sol`: Solvent diffusivity in film (m²/s)
  
- `p.c_sol`: Solvent concentration (mol/m³)
  
- `p.U`: SEI formation potential (V vs Li/Li⁺)
  
- `p.R`: Film resistivity (Ω·m)
  

**Output Variables**
- `L_sei`: SEI film thickness (grows over time, asymptotically)
  
- `j_sei`: SEI current density (decreases as L_sei increases)
  
- `ϕf`: Film potential drop (increases with thickness)
  

**Physical Assumption**

Film diffusion resistance and potential drop dominate; SEI growth self-limits via thickness.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/SEI.jl#L163-L191" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Lithium Plating {#Lithium-Plating}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.LithiumPlating.NoPlating' href='#BatteryToolkit.LithiumPlating.NoPlating'><span class="jlbinding">BatteryToolkit.LithiumPlating.NoPlating</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
NoPlating(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
```


Create a zero lithium plating model (reaction disabled).

Returns a ModelingToolkit system where lithium plating is disabled. Use this when lithium plating effects are negligible or you want to exclude them from the simulation.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::SideReactionParameters`: Lithium plating reaction parameters (unused)
  
- `s::SolidParticleParameters`: Electrode solid particle parameters
  
- `g`: FVM geometry object
  

**Output Variables**
- `c_plating`: Plated lithium concentration (always 0)
  
- `c_dead`: Dead lithium concentration (always 0)
  
- `j_stripping`: Lithium stripping current (always 0)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/LithiumPlating.jl#L8-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.LithiumPlating.IrreversiblePlating' href='#BatteryToolkit.LithiumPlating.IrreversiblePlating'><span class="jlbinding">BatteryToolkit.LithiumPlating.IrreversiblePlating</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
IrreversiblePlating(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
```


Create an irreversible lithium plating model.

Models lithium plating where all plated lithium becomes "dead" (electrochemically inactive). Once plated, lithium cannot be stripped (reversed). Use for conservative simulations where plating is the primary failure mechanism.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::SideReactionParameters`: Lithium plating reaction parameters
  
- `s::SolidParticleParameters`: Electrode solid particle parameters
  
- `g`: FVM geometry object
  

**Key Physics**
- Plating occurs when surface potential drops too low (overpotential driven)
  
- All plated lithium irreversibly becomes dead lithium
  
- Dead lithium consumes cyclable lithium inventory
  
- Stripping current is zero (irreversible assumption)
  

**Output Variables**
- `c_plating`: Active plated lithium (accumulates until fully dead)
  
- `c_dead`: Dead lithium (irreversibly lost from battery)
  
- `j_stripping`: Always zero (no reversibility)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/LithiumPlating.jl#L82-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.LithiumPlating.PartiallyReversiblePlating' href='#BatteryToolkit.LithiumPlating.PartiallyReversiblePlating'><span class="jlbinding">BatteryToolkit.LithiumPlating.PartiallyReversiblePlating</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
PartiallyReversiblePlating(; name, p::SideReactionParameters, s::SolidParticleParameters, g)
```


Create a partially reversible lithium plating model.

Models lithium plating where a fraction of plated lithium can be stripped (reversed) during charging, but some becomes dead due to structural damage or electrolyte reactions. Use for realistic simulations where partial reversibility and cycling damage occur.

**Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `p::SideReactionParameters`: Lithium plating reaction and reversibility parameters
  
- `s::SolidParticleParameters`: Electrode solid particle parameters
  
- `g`: FVM geometry object
  

**Key Physics**
- Plating occurs during overdischarge (high negative overpotential)
  
- Stripping (reversal) occurs during charging with kinetic limitations
  
- Some plated lithium irreversibly becomes dead due to cycling stress
  
- Dead lithium fraction accumulates with plating/stripping cycles
  

**Output Variables**
- `c_plating`: Active plated lithium (can grow or shrink with cycling)
  
- `c_dead`: Dead lithium (irreversibly lost, increases over time)
  
- `j_stripping`: Stripping current (depends on plating thickness and potential)
  

**Cycling Behavior**

Repeated plating/stripping cycles increase dead lithium fraction, modeling accelerated degradation under abuse conditions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/CellModels/SPMe/LithiumPlating.jl#L193-L222" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Potential Calculations {#Potential-Calculations}

Electrochemical potential and overpotential functions used internally in SPMe simulations.
