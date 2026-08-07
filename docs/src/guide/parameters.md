# Battery Parameter Sets

BatteryToolkit provides pre-configured parameter sets for lithium-ion battery simulations.

## Overview

Two main parameter sets are available:

| Aspect | Chen2020 | OKane2022 |
|--------|----------|-----------|
| **Cell Type** | Graphite/NMC pouch cell | Graphite/NMC pouch cell |
| **Capacity** | 5 Ah nominal | 5 Ah nominal |
| **Use Case** | General-purpose simulations | High-fidelity research |
| **Fit Quality** | Good for typical conditions | Excellent across wide range |
| **Negative Electrode** | Graphite with SEI | Graphite with SEI |
| **Positive Electrode** | NMC (no side reactions) | NMC (no side reactions) |
| **Calibration Data** | Single test point | Multi-point characterization |
| **Experimental Basis** | Standard lithium-ion | Detailed electrochemistry |

## Creating Parameter Sets

### Chen2020 Parameters

```julia
params = Chen2020()
```

**Suitable for:**
- General battery simulations
- Initial feasibility studies
- Teaching and learning
- Rapid prototyping

### OKane2022 Parameters

```julia
params = OKane2022()
```

**Suitable for:**
- High-fidelity research
- Detailed performance prediction
- Physics-based design
- Validation studies

## Using Parameter Sets

Both parameter sets are used identically:

```julia
params = Chen2020()
sys = SPMe(params=params, Q=5)
exp = Experiment([PowerStep(1000, 3600)])
sol = simulate(sys, exp, Rodas4())
```

## Parameter Structure

Each parameter set returns a `BatteryParameters` object containing:

### Electrode Parameters

**Negative Electrode (graphite anode)**
- Particle radius and surface area
- Lithium diffusivity and conductivity
- Initial and max lithium concentration
- Open-circuit potential curve
- SEI (Solid Electrolyte Interface) film

**Positive Electrode (NMC cathode)**
- Particle radius and surface area
- Lithium diffusivity and conductivity
- Initial and max lithium concentration
- Open-circuit potential curve

### Electrolyte Parameters

- Electrolyte thickness in each domain (negative, separator, positive)
- Ion diffusivity (concentration-dependent)
- Ionic conductivity (concentration-dependent)
- Transference number
- Porosity in each domain
- Bruggeman tortuosity coefficient

### Cell Geometry

- Current collector dimensions: 0.065 m height, 1.58 m width
- Number of parallel electrode pairs
- Voltage limits: 2.0 V min, 4.2 V max

## Customizing Parameters

To create custom parameters, construct a `BatteryParameters` object:

```julia
# Define electrode parameters
n = SolidParticleParameters(
    Rₖ = 5e-6,              # Particle radius (m)
    aₖ = 200000,            # Surface area density (m⁻¹)
    Dₖ = c -> 1e-14,        # Diffusivity (m²/s)
    σₖ = 100,               # Conductivity (S/m)
    c₀ = 30000,             # Initial concentration (mol/m³)
    c₊ = 33000,             # Max concentration (mol/m³)
    Uₖ = z -> -0.8*z + 4.5, # OCP as function of stoichiometry
    mₖ = 1e-6,              # Reaction rate constant
    L_sei₀ = 1e-8,
)

# Define electrolyte parameters
e = ElectrolyteParameters(
    Lₚ = 75e-6,
    Lₛ = 12e-6,
    Lₙ = 85e-6,
    Dₑ = c -> 1e-10,
    σₑ = c -> 0.5,
    c₀ = 1000.0,
    t₊ = 0.26,
    ϵₚ = 0.3, ϵₛ = 0.5, ϵₙ = 0.3,
    bₚ = 1.5, bₛ = 1.5, bₙ = 1.5,
)

# Assemble into BatteryParameters
params = BatteryParameters(
    p = positive_electrode,
    n = negative_electrode,
    e = electrolyte,
    Hcc = 0.065, Wcc = 1.58,
    n_el = 1, Q₀ = 5,
    Vmin = 2.0, Vmax = 4.2
)
```

## Next Steps

- [Experiments](experiments.md): Define complex test profiles
- [Examples](examples.md): See parameter sets in action
- [SPMe Model](../models/spme.md): Understand the physics
