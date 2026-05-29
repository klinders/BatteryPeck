# Side Reactions: SEI and Lithium Plating

Model secondary electrochemical reactions that drive battery degradation.

## Overview

Side reactions represent failure mechanisms in lithium-ion batteries:

1. **SEI (Solid Electrolyte Interface)**: Passivating film on negative electrode
   - Grows during discharge
   - Increases impedance and reduces efficiency
   - Main source of capacity fade in early cycling

2. **Lithium Plating**: Metallic lithium deposition on negative electrode
   - Occurs during overdischarge or high-rate charging
   - Creates safety hazards
   - Main failure mode in abuse scenarios

## Solid Electrolyte Interface (SEI)

### Formation and Growth

SEI forms when electrolyte reacts with lithium metal on the negative electrode surface:

```
Li + organic electrolyte → SEI products (Li₂CO₃, LiF, etc.)
```

**Key characteristics:**
- Grows continuously during discharge (small rate)
- Ionic conductor (allows Li⁺ transport)
- Electronic insulator (prevents electron flow)
- Increases cell impedance over time

### SEI Models in BatteryToolkit

BatteryToolkit provides three SEI implementations with increasing realism:

#### 1. NoSEI

```julia
sys = SPMe(params=params, side_reactions=:NoSEI)
```

**Features:**
- No SEI film formation
- Zero resistance contribution
- Fastest simulation

**Use when:**
- SEI effects are negligible
- Simulating fresh cells on short timescales
- Exploring other phenomena
- Fast feasibility studies

**Output:**
```
L_sei = 0 (always)
j_sei = 0 (always)
ϕf = 0 (always)
```

#### 2. ReactionLimitedSEI

```julia
sys = SPMe(params=params, side_reactions=:ReactionLimitedSEI)
```

**Physics:**
- SEI current follows Butler-Volmer kinetics
- Growth rate controlled by surface overpotential
- Film acts as perfect ionic conductor

**Mathematical form:**
```
j_sei = j₀ exp(-α·F·η_sei/R/T)
L_sei = dL_sei/dt = j_sei·V̄/(F·aₖ)
```

**Characteristics:**
- Rapid initial SEI growth
- Growth slows as overpotential decreases
- Asymptotic limit as η_sei → 0

**Use when:**
- SEI growth is reaction-dominated
- Fast formation conditions (high overpotential)
- Typical discharge scenarios
- Most practical applications

**Parameters:**
- `j_sei₀`: Exchange current density (A/m²)
- `α`: Transfer coefficient (0-1)
- `U`: SEI formation potential (V)

**Output:**
```
L_sei: Growing thickness (m)
j_sei: Reaction current (A/m²)
ϕf: Usually small (< 0.01V)
```

#### 3. SolventDiffusionLimitedSEI

```julia
sys = SPMe(params=params, side_reactions=:SolventDiffusionLimitedSEI)
```

**Physics:**
- SEI current limited by solvent diffusion through film
- Reaction fast, diffusion slow
- Film provides significant resistance

**Mathematical form:**
```
j_sei ≈ D_sol·c_sol / (L_sei·V̄)
L_sei = dL_sei/dt = j_sei·V̄/(F·aₖ)
ϕf = j_app·L_sei·R
```

**Characteristics:**
- Rapid initial growth (diffusion not yet limiting)
- Gradual slowdown as thickness increases
- Growth asymptotically approaches zero
- Significant ohmic potential drop

**Use when:**
- SEI forms thick protective layer
- Long-term cycling/calendar aging
- Conservative capacity fade predictions
- Higher impedance important

**Parameters:**
- `D_sol`: Solvent diffusivity in film (m²/s)
- `c_sol`: Solvent concentration (mol/m³)
- `R`: Film resistivity (Ω·m)
- `V̄`: Partial molar volume (m³/mol)

**Output:**
```
L_sei: Growing thickness, self-limiting
j_sei: Decreasing current
ϕf: Increasing potential drop (0.01-0.1V)
```

## Lithium Plating

### Formation and Stripping

Lithium plating occurs when the negative electrode potential drops below 0V (vs. Li/Li⁺):

```
Li⁺ + e⁻ → Li (plating, potential too low)
Li → Li⁺ + e⁻ (stripping, during charging)
```

**Failure mechanisms:**
- Metallic lithium forms dendrites (safety hazard)
- Lost lithium ("dead lithium") reduces capacity
- Separator puncture from dendrites → short circuit

### Lithium Plating Models

BatteryToolkit provides three plating implementations:

#### 1. NoPlating

```julia
sys = SPMe(params=params, side_reactions=:NoPlating)
```

**Features:**
- No lithium plating
- Zero plating current

**Use when:**
- Operating within safe voltage window
- No abuse scenarios
- Fast simulations

**Output:**
```
c_plating = 0
c_dead = 0
j_stripping = 0
```

#### 2. IrreversiblePlating

```julia
sys = SPMe(params=params, side_reactions=:IrreversiblePlating)
```

**Physics:**
- All plated lithium becomes permanently inactive
- Represents worst-case capacity loss
- No stripping during charging

**Characteristics:**
- Rapid capacity fade if plating occurs
- Conservative estimate of degradation
- Useful for safety analysis

**Use when:**
- Simulating abuse conditions
- Conservative degradation predictions
- Analyzing overdischarge scenarios

**Output:**
```
c_plating: Grows when potential drops below 0V
c_dead: Plated Li → dead Li (monotonically increasing)
j_stripping = 0 (no reversibility)
```

#### 3. PartiallyReversiblePlating

```julia
sys = SPMe(params=params, side_reactions=:PartiallyReversiblePlating)
```

**Physics:**
- Some plated lithium can be stripped (reversed)
- Fraction becomes irreversibly dead
- Stripping limited by kinetics

**Characteristics:**
- More realistic cycling behavior
- Dead lithium accumulates over cycles
- Reversible plating enables recovery

**Use when:**
- Realistic abuse + recovery scenarios
- Multi-cycle degradation studies
- Detailed capacity fade modeling
- Realistic cycle life prediction

**Output:**
```
c_plating: Reversible plated lithium
c_dead: Dead lithium (irreversible loss)
j_stripping: Stripping current (kinetically limited)
```

## Configuration Examples

### Conservative Simulation (Fast)

```julia
sys = SPMe(
    params=Chen2020(),
    N=10,
    side_reactions=:NoSEI  # Ignore SEI
)
# Fastest, ignores degradation
```

### Realistic Simulation (Balanced)

```julia
sys = SPMe(
    params=OKane2022(),
    N=15,
    side_reactions=:ReactionLimitedSEI  # Include SEI growth
)
# Good speed/fidelity balance
```

### High-Fidelity Simulation (Slow)

```julia
sys = SPMe(
    params=OKane2022(),
    N=20,
    side_reactions=:SolventDiffusionLimitedSEI  # Include detailed SEI
)
# Best accuracy, slowest
```

### Abuse Scenario (Plating Risk)

```julia
sys = SPMe(
    params=Chen2020(),
    N=10,
    side_reactions=:PartiallyReversiblePlating  # Realistic plating
)
# For overdischarge/overcharge studies
```

## Impact on Battery Performance

### SEI Effects

**ReactionLimitedSEI:**
```
Discharge 1: Coulombic efficiency ≈ 99.5%
Discharge 2: Coulombic efficiency ≈ 99.7%
Discharge 10: Coulombic efficiency ≈ 99.8%
```

**SolventDiffusionLimitedSEI:**
```
Discharge 1: Coulombic efficiency ≈ 97%
Discharge 2: Coulombic efficiency ≈ 98%
Discharge 10: Coulombic efficiency ≈ 98.5%
(Higher impedance, lower efficiency)
```

### Plating Effects

**Without plating:**
```
Overdischarge: Stops at Vmin (2.0V)
Capacity: Unchanged
```

**With IrreversiblePlating:**
```
Overdischarge: Capacity loss
Stripping impossible: Lost lithium gone
```

**With PartiallyReversiblePlating:**
```
Overdischarge: Partial capacity loss
Charge: Lithium stripped back (partially recoverable)
```

## Performance Tips

### Speed Optimization

```julia
# Fastest: No side reactions
sys = SPMe(params=params, N=10, side_reactions=:NoSEI)
# ~10-30 seconds per 1-hour simulation

# Medium: Simple SEI
sys = SPMe(params=params, N=10, side_reactions=:ReactionLimitedSEI)
# ~15-40 seconds per 1-hour simulation

# Slow: Complex SEI
sys = SPMe(params=params, N=10, side_reactions=:SolventDiffusionLimitedSEI)
# ~20-60 seconds per 1-hour simulation
```

### Accuracy vs Speed Trade-off

```julia
# Fast, less accurate
fast_config = SPMe(params=params, N=5, side_reactions=:ReactionLimitedSEI)

# Balanced (default)
balanced_config = SPMe(params=params, N=10, side_reactions=:ReactionLimitedSEI)

# Slow, accurate
accurate_config = SPMe(params=params, N=20, side_reactions=:SolventDiffusionLimitedSEI)
```

## Next Steps

- [SPMe Model](spme.md): Core electrochemistry
- [Examples](../guide/examples.md): Degradation in practice
- [API Reference](../api/cellmodels.md): Complete model API
