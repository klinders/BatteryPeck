# Quick Start

Get up and running with BatteryToolkit in 5 minutes.

## Installation

Add the package to your Julia environment:

```julia
using Pkg
Pkg.add("BatteryToolkit")
```

## Your First Simulation

Here's a minimal working example that creates a battery cell, defines an experiment, and runs a simulation:

```julia
using BatteryToolkit

# 1. Load battery parameters (Chen et al. 2020)
params = Chen2020()

# 2. Create the battery cell model (SPMe = Single Particle Model with Electrolyte)
sys = SPMe(params=params, Q=5)

# 3. Define an experiment: discharge at 1000W for 1 hour
experiment = Experiment([PowerStep(1000, 3600)])

# 4. Run the simulation
using OrdinaryDiffEq
sol = simulate(sys, experiment, Rodas4())

# 5. Extract and plot results
time = sol.t
voltage = sol[sys.Pin]
current = sol[sys.Iin]
```

That's it! You now have voltage and current time series data for your battery simulation.

## What's Next?

- **Learn about parameter sets**: See [Parameter Sets](parameters.md) to understand Chen2020 vs OKane2022
- **Explore experiments**: See [Experiments](experiments.md) for multi-step profiles, drive cycles, and charging
- **View worked examples**: See [Examples](../guide/examples.md) for common simulation scenarios
- **Dive into the models**: See [SPMe Model](../models/spme.md) for electrochemistry details
- **API reference**: See [API Reference](../api/parameters.md) for complete function documentation

## Common Tasks

### Run a multi-step experiment

```julia
steps = [
    PowerStep(1000, 1800),    # 1000W for 30 minutes
    RestStep(300),             # 5 minute rest
    PowerStep(-500, 3600)      # Discharge at 500W for 1 hour
]
exp = Experiment(steps)
sol = simulate(sys, exp, Rodas4())
```

### Use different parameter sets

```julia
# Chen et al. 2020 (typical graphite/NMC)
params_chen = Chen2020()

# O'Kane et al. 2022 (high-fidelity graphite/NMC)
params_okane = OKane2022()

sys1 = SPMe(params=params_chen)
sys2 = SPMe(params=params_okane)
```

### Increase model resolution

```julia
# Use more finite volume nodes for finer spatial resolution
sys = SPMe(params=params, N=20)  # Default is N=10
```

## Need Help?

- Check [Experiments](experiments.md) for step types and experiment composition
- See [Examples](../guide/examples.md) for real-world simulation scenarios
- Read the full [API Reference](../api/parameters.md)
