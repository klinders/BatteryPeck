
# BatteryToolkit.jl

A Julia package for electrochemical battery modeling and simulation using the Single Particle Model with Electrolyte (SPMe).

## Features

- **SPMe Physics**: Single particle model with electrolyte transport, charge transfer kinetics, and potential calculations
- **Degradation Models**: SEI film growth and lithium plating representations
- **Multiple Parameter Sets**: Chen2020 and O'Kane2022 calibrated for graphite/NMC batteries
- **Pack Modeling**: Series and parallel battery pack configurations
- **Flexible Experiments**: Power, current, rest, charging, and driving cycle profiles
- **High Performance**: Efficient finite volume discretization with ModelingToolkit integration

## Quick Start

```julia
using BatteryToolkit, OrdinaryDiffEq

# Create battery system
params = Chen2020()
sys = SPMe(params=params)

# Define experiment
exp = Experiment([PowerStep(1000, 3600)])  # 1000W for 1 hour

# Run simulation
sol = simulate(sys, exp, Rodas4())

# Extract results
voltage = sol[sys.Pin]
current = sol[sys.Iin]
```

See the [Quick Start](guide/quickstart.md) for more details.

## Documentation

- **[Getting Started](guide/quickstart.md)**: Installation and first steps
- **[Parameter Sets](guide/parameters.md)**: Available battery parameter sets
- **[Experiments](guide/experiments.md)**: Define test profiles and protocols
- **[Examples](guide/examples.md)**: Real-world simulation scenarios
- **[SPMe Model](models/spme.md)**: Electrochemistry and model structure
- **[API Reference](api/parameters.md)**: Complete function documentation

## Installation

```julia
using Pkg
Pkg.add("BatteryToolkit")
```

## Supported Julia Versions

Julia 1.8 and later

## Citation

If you use BatteryToolkit in your research, please cite this package and the underlying parameter sets.

## License

See LICENSE file for details.

```@bibliography
```
