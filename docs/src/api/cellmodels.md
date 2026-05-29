# Cell Models API Reference

API documentation for SPMe and related electrochemical models.

## Main Cell Model

```@docs
BatteryToolkit.SPMe
```

## SPMe Components

### Solid Particle Model

```@docs
BatteryToolkit.SolidParticle
```

### Electrolyte Model

```@docs
BatteryToolkit.Electrolyte
```

## Side Reactions

### SEI (Solid Electrolyte Interface)

```@docs
BatteryToolkit.SEI.NoSEI
BatteryToolkit.SEI.ReactionLimitedSEI
BatteryToolkit.SEI.SolventDiffusionLimitedSEI
```

### Lithium Plating

```@docs
BatteryToolkit.LithiumPlating.NoPlating
BatteryToolkit.LithiumPlating.IrreversiblePlating
BatteryToolkit.LithiumPlating.PartiallyReversiblePlating
```

## Potential Calculations

Electrochemical potential and overpotential functions used internally in SPMe simulations.
