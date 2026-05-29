# Examples and Worked Scenarios

Practical examples showing common battery simulation use cases.

## Example 1: Simple Constant-Power Discharge

The most basic scenario: discharge a fully charged battery at constant power.

```julia
using BatteryToolkit
using OrdinaryDiffEq
using Plots

# Setup
params = Chen2020()
sys = SPMe(params=params)

# Discharge at 1000W for 1 hour
exp = Experiment([PowerStep(1000, 3600)])

# Run simulation
sol = simulate(sys, exp, Rodas4())

# Extract and plot
time = sol.t
voltage = sol[sys.Pin]
current = sol[sys.Iin]

plot(time, voltage, label="Voltage", ylabel="Voltage (V)")
```

**Expected Results:**
- Voltage starts at ~3.8V and drops to ~2.0V
- Simulation stops when reaching Vmin (safety limit)

## Example 2: Multi-Step Charging and Discharging

A realistic test protocol with charging, resting, and discharge phases.

```julia
using BatteryToolkit
using OrdinaryDiffEq

params = OKane2022()  # Higher fidelity
sys = SPMe(params=params, N=15)  # Use 15 nodes per domain

steps = [
    RestStep(300),                    # 5 min initialization rest
    PowerStep(5000, 3600),            # Charge at 5000W for 1 hour
    RestStep(600),                    # 10 min cool-down
    PowerStep(-3000, 7200),           # Discharge at 3000W for 2 hours
    RestStep(300),                    # 5 min final rest
]

exp = Experiment(steps)
sol = simulate(sys, exp, Rodas4())

# Access results
time = sol.t
voltage = sol[sys.Pin]
energy_out = sol[sys.cell.energy_out]  # If available

println("Simulation duration: $(time[end]/3600) hours")
println("Final voltage: $(voltage[end]) V")
```

## Example 3: Comparing Parameter Sets

Evaluate performance differences between Chen2020 and OKane2022 on the same scenario.

```julia
using BatteryToolkit
using OrdinaryDiffEq
using Plots

# Define test profile
steps = [
    PowerStep(2000, 1800),      # 2000W for 30 min
    RestStep(300),               # 5 min rest
    PowerStep(-1000, 3600),      # 1000W discharge for 1 hour
]

# Run with Chen2020
params_c = Chen2020()
sys_c = SPMe(params=params_c)
exp = Experiment(steps)
sol_c = simulate(sys_c, exp, Rodas4())

# Run with OKane2022
params_o = OKane2022()
sys_o = SPMe(params=params_o)
sol_o = simulate(sys_o, exp, Rodas4())

# Compare voltages
p = plot(
    sol_c.t, sol_c[sys_c.Pin], label="Chen2020",
    ylabel="Voltage (V)", xlabel="Time (s)"
)
plot!(p, sol_o.t, sol_o[sys_o.Pin], label="OKane2022")
display(p)

# Print differences
println("Max voltage difference: $(maximum(abs.(sol_c[sys_c.Pin] .- sol_o[sys_o.Pin]))) V")
```

## Example 4: Drive Cycle Simulation

Apply a real-world driving cycle from CSV data.

**Example CSV file (drive_cycle.csv):**
```
time,power
0,0
60,1000
120,2000
180,3000
240,2000
300,1000
360,0
420,-500
480,-1000
540,-1500
600,0
```

**Julia code:**
```julia
using BatteryToolkit
using OrdinaryDiffEq

params = Chen2020()
sys = SPMe(params=params)

# Apply drive cycle with 1 hour maximum
exp = Experiment([DriveStep("drive_cycle.csv", 3600)])
sol = simulate(sys, exp, Rodas4())

# Analyze
time = sol.t
voltage = sol[sys.Pin]
avg_voltage = mean(voltage)

println("Drive cycle completed in $(time[end]) seconds")
println("Average voltage during cycle: $avg_voltage V")
```

## Example 5: Temperature Effects

Study battery performance at different temperatures (if available in model).

```julia
using BatteryToolkit
using OrdinaryDiffEq

params = OKane2022()
sys = SPMe(params=params)

# Standard discharge
exp = Experiment([PowerStep(1000, 3600)])

# Note: Temperature effects may require custom system setup
# This example shows the structure
sol = simulate(sys, exp, Rodas4())

time = sol.t
voltage = sol[sys.Pin]
```

## Example 6: Long-Duration Cycling

Simulate multiple charge-discharge cycles to study degradation (SEI growth).

```julia
using BatteryToolkit
using OrdinaryDiffEq

params = Chen2020()
sys = SPMe(params=params)

# Single cycle: charge and discharge
single_cycle = [
    PowerStep(2000, 1800),       # Charge 30 min
    RestStep(300),               # Rest 5 min
    PowerStep(-1000, 3600),      # Discharge 1 hour
    RestStep(300),               # Rest 5 min
]

# Repeat 10 cycles
cycles = single_cycle * 10

exp = Experiment(cycles)
sol = simulate(sys, exp, Rodas4())

# SEI thickness growth (if enabled)
# time_end_of_cycle_1 = cumsum([s.period for s in single_cycle])[end]
# time_end_of_cycle_2 = 2 * time_end_of_cycle_1
# etc.

println("Total simulation time: $(sol.t[end]/3600) hours")
```

## Tips and Tricks

### Adjust Solver Tolerances for Speed

```julia
# Loose tolerances for quick estimates (faster)
sol = simulate(sys, exp, Rodas4(), abstol=1e-6, reltol=1e-4)

# Tight tolerances for accurate results (slower but more precise)
sol = simulate(sys, exp, Rodas4(), abstol=1e-10, reltol=1e-8)
```

### Increase Model Resolution for Better Accuracy

```julia
# Default resolution (10 nodes per domain)
sys = SPMe(params=params)

# Higher resolution (20 nodes) - slower but more accurate
sys = SPMe(params=params, N=20)
```

### Extract Specific Variables

```julia
# Access state variables from solution
voltage = sol[sys.Pin]          # Battery voltage
current = sol[sys.Iin]          # Battery current
soc = sol[sys.cell.soc]         # State of charge (if defined)
```

### Save Results for Analysis

```julia
using JLD2

sol = simulate(sys, exp, Rodas4())

# Save for later use
@save "my_simulation.jld2" sol time voltage current

# Load later
@load "my_simulation.jld2" sol time voltage current
```

## Next Steps

- [Experiments](experiments.md): Learn more about step types
- [Parameters](parameters.md): Customize parameter sets
- [SPMe Model](../models/spme.md): Understand the physics
- [API Reference](../api/parameters.md): Browse complete function documentation
