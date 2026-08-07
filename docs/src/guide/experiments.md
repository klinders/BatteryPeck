# Experiments and Step Types

Design complex battery test profiles by composing step types into experiments.

## Overview

An `Experiment` combines multiple `Step` objects into a single test profile. Each step defines
a distinct phase of operation (charging, discharging, rest, or custom cycle).

## Step Types

### PowerStep

Apply constant power for a specified duration.

```julia
PowerStep(power::Real, period::Real)
```

**Arguments:**
- `power`: Power in watts (W). Negative for discharge.
- `period`: Duration in seconds

**Example:**
```julia
PowerStep(1000, 3600)   # 1000W for 1 hour
PowerStep(-500, 1800)   # Discharge at 500W for 30 minutes
```

**Use case:** Constant power charging/discharging

### CurrentStep

Apply constant current for a specified duration.

```julia
CurrentStep(current::Real, period::Real)
```

**Arguments:**
- `current`: Current in amperes (A). Negative for discharge.
- `period`: Duration in seconds

**Example:**
```julia
CurrentStep(10, 3600)   # 10A for 1 hour
CurrentStep(-10, 1800)  # Discharge at 10A for 30 minutes
```

**Use case:** Constant current (CC) charging or discharging

### ChargeStep

Charge to a target state-of-charge (SoC) using constant power.

```julia
ChargeStep(soc::Real, period::Real=0, power::Real=11000)
```

**Arguments:**
- `soc`: Target state-of-charge (0 to 1)
- `period`: Maximum duration in seconds (0 = unlimited)
- `power`: Charging power in watts (default: 11000W)

**Example:**
```julia
ChargeStep(0.8, 3600, 5000)  # Charge to 80% SoC within 1 hour at 5000W
ChargeStep(1.0)              # Charge to 100% SoC with default 11000W
```

**Use case:** Constant-power charging to target SoC

### RestStep

Hold the battery at open-circuit (no charge/discharge) for a rest period.

```julia
RestStep(period::Real)
```

**Arguments:**
- `period`: Rest duration in seconds

**Example:**
```julia
RestStep(300)  # 5-minute rest
RestStep(3600) # 1-hour rest
```

**Use case:** Relaxation, thermal stabilization, OCV measurement

### DriveStep

Apply a driving cycle from a CSV file with time and power columns.

```julia
DriveStep(file::String, period::Real=nothing)
```

**Arguments:**
- `file`: Path to CSV file with columns [time, power]
- `period`: Maximum duration in seconds (default: full duration from CSV)

**Example:**
```julia
DriveStep("drive_cycle.csv", 3600)  # Apply cycle for 1 hour maximum
```

**CSV Format:**
```
time,power
0,0
10,500
20,1000
...
```

**Use case:** Real-world driving cycles, custom power profiles

## Composing Experiments

### Single Step

```julia
exp = Experiment([PowerStep(1000, 3600)])
sol = simulate(sys, exp, Rodas4())
```

### Multiple Steps (sequential)

```julia
steps = [
    PowerStep(1000, 1800),   # 1000W for 30 min
    RestStep(300),            # 5 min rest
    PowerStep(-500, 3600),    # 500W discharge for 1 hour
]
exp = Experiment(steps)
sol = simulate(sys, exp, Rodas4())
```

### Repeat Steps

Use `*` operator to repeat steps:

```julia
base_step = [PowerStep(1000, 1800), RestStep(300)]
steps = base_step * 3  # Repeat 3 times
exp = Experiment(steps)
```

### Concatenate Step Sequences

Use `+` operator to combine sequences:

```julia
charge_steps = [PowerStep(1000, 3600)]
discharge_steps = [PowerStep(-500, 3600)]
steps = charge_steps + discharge_steps
exp = Experiment(steps)
```

## Experiment Properties

Once created, an `Experiment` object contains:

| Property | Type | Description |
|----------|------|-------------|
| `steps` | Vector{AbstractStep} | All steps in order |
| `tstops` | Vector{Float64} | Time at end of each step except last (s) |
| `tend` | Float64 | Total simulation duration (s) |
| `step_count` | Int64 | Number of steps |
| `p0` | Float64 | Initial power/current (W or A) |
| `start_time` | DateTime | Experiment start timestamp |

## Example: Complete Test Protocol

```julia
using BatteryToolkit
using OrdinaryDiffEq

# Create system
params = OKane2022()
sys = SPMe(params=params, N=15)

# Define multi-stage test
steps = [
    # Initialization
    RestStep(60),
    
    # Charging
    PowerStep(5000, 7200),      # 5000W for 2 hours
    RestStep(300),               # Cool down 5 min
    
    # Discharging
    PowerStep(-3000, 10800),     # 3000W discharge for 3 hours
    RestStep(600),               # Cool down 10 min
    
    # Repeat 3 times
] * 3

# Run simulation
exp = Experiment(steps)
sol = simulate(sys, exp, Rodas4())

# Extract results
time = sol.t
voltage = sol[sys.Pin]
```

## Tips

- **Timestamps**: Use `Experiment(steps, DateTime(2024,1,15))` to set experiment start time
- **Long simulations**: Use larger tolerances (`abstol=1e-8, reltol=1e-6`) for faster computation
- **Safety limits**: Battery stops if voltage hits Vmin or Vmax
- **High fidelity**: Increase `N` parameter in SPMe for better spatial resolution

## Next Steps

- [Examples](examples.md): See real-world experiment scenarios
- [Quick Start](quickstart.md): Return to basics
- [SPMe Model](../models/spme.md): Understand underlying physics
