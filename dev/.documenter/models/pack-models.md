
# Pack Models {#Pack-Models}

Extend single-cell simulations to series and parallel battery pack configurations.

## Overview {#Overview}

BatteryToolkit provides two pack models for simulating multi-cell batteries:
1. **SingleCellPack**: Single cell with series/parallel electrode configuration
  
2. **MultiCellPack**: Multiple cells connected in series and/or parallel
  

## SingleCellPack {#SingleCellPack}

Models a single cell with configurable internal series/parallel electrode arrangement.

### Configuration {#Configuration}

```julia
SingleCellPack(;
    name,
    params=Chen2020(),
    config=(96, 3),  # (n_series, n_parallel)
    Qcell=5
)
```


### Arguments {#Arguments}
- `name`: System name (required, use `Symbol` like `:pack`)
  
- `params`: Battery parameters (default: Chen2020)
  
- `config`: Tuple of (n_series, n_parallel) electrodes
  
- `Qcell`: Cell capacity in Ah (default: 5)
  

### Example {#Example}

```julia
using BatteryToolkit

# Single cell with 96 series, 3 parallel electrodes
pack = SingleCellPack(
    name=:battery,
    params=Chen2020(),
    config=(96, 3),
    Qcell=5
)

# Run simulation
exp = Experiment([PowerStep(1000, 3600)])
sol = simulate(pack, exp, Rodas4())

# Extract results
voltage = sol[pack.V]
current = sol[pack.I]
```


### Electrical Topology {#Electrical-Topology}

For `config=(96, 3)`:
- 96 electrodes in series → voltage multiplied by 96
  
- 3 electrode pairs in parallel → current divided by 3
  
- Overall: high voltage, moderate current
  

**Voltage scaling:**

```
V_total = config[1] * V_single_electrode
```


**Current scaling:**

```
I_total = config[2] * I_single_electrode
```


### Input/Output Signals {#Input/Output-Signals}

**Input ports:**
- `Pin`: Input power (W)
  
- `Iin`: Input current (A)
  
- `Tin`: Input temperature (K)
  

**Output variables:**
- `V`: Pack voltage (V)
  
- `I`: Pack current (A)
  

## MultiCellPack {#MultiCellPack}

Models multiple battery cells connected in series and/or parallel.

### Configuration {#Configuration-2}

```julia
MultiCellPack(;
    name,
    params=Chen2020(),
    config=(12, 3),  # (n_series, n_parallel)
    Qcell=5
)
```


### Arguments {#Arguments-2}
- `name`: System name (required)
  
- `params`: Battery parameters
  
- `config`: Tuple of (n_series, n_parallel) cells
  
- `Qcell`: Single cell capacity (Ah)
  

### Example {#Example-2}

```julia
using BatteryToolkit

# 12 cells in series, 3 modules in parallel
# Total: 36 cells arranged as (12s × 3p)
pack = MultiCellPack(
    name=:battery_pack,
    params=OKane2022(),
    config=(12, 3),
    Qcell=5
)

exp = Experiment([PowerStep(5000, 3600)])
sol = simulate(pack, exp, Rodas4())

voltage = sol[pack.V]  # ~48V for 12s (4V per cell)
current = sol[pack.I]  # Distributed across 3p
```


### Electrical Topology {#Electrical-Topology-2}

For `config=(12, 3)`:
- 12 cells in series → 12 × single-cell voltage
  
- 3 parallel branches → current split equally
  
- Total cells: 12 × 3 = 36 cells
  

**Total voltage:**

```
V_total = config[1] * V_single_cell
```


**Total current:**

```
I_total = config[2] * I_single_cell
```


**Total capacity:**

```
Q_total = config[2] * Q_single_cell  # Parallel adds capacity
```


### Input/Output Signals {#Input/Output-Signals-2}

**Input ports:**
- `P`: Input power (W)
  
- `T`: Input temperature (K)
  

**Output variables:**
- `V`: Total pack voltage (V)
  
- `I`: Total pack current (A)
  

## Comparison: SingleCellPack vs MultiCellPack {#Comparison:-SingleCellPack-vs-MultiCellPack}

|               Aspect |                        SingleCellPack |                    MultiCellPack |
| --------------------:| -------------------------------------:| --------------------------------:|
|          **Physics** | Single cell with internal arrangement |       Multiple independent cells |
|            **Cells** |             One (internal electrodes) | Multiple (config[1] × config[2]) |
|         **Use Case** |                  Single-cell research |          Multi-cell battery pack |
|          **Voltage** |               config[1] × V_electrode |               config[1] × V_cell |
|          **Current** |               config[2] × I_electrode |               config[2] × I_cell |
|       **Complexity** |                                 Lower |                           Higher |
| **Simulation Speed** |                                Faster |                           Slower |


## Pack Configuration Examples {#Pack-Configuration-Examples}

### Example 1: Small Pack (12V nominal) {#Example-1:-Small-Pack-12V-nominal}

```julia
# 4 cells in series, 1 parallel = 12V, 5Ah
pack = MultiCellPack(
    name=:small_pack,
    params=Chen2020(),
    config=(4, 1),
    Qcell=5
)
# Nominal voltage: ~16.8V (4.2V × 4)
# Capacity: 5Ah
```


### Example 2: Medium Pack (48V nominal) {#Example-2:-Medium-Pack-48V-nominal}

```julia
# 12 cells in series, 3 parallel = 48V, 15Ah
pack = MultiCellPack(
    name=:medium_pack,
    params=OKane2022(),
    config=(12, 3),
    Qcell=5
)
# Nominal voltage: ~50.4V (4.2V × 12)
# Capacity: 15Ah (5Ah × 3)
```


### Example 3: Large Pack (400V nominal) {#Example-3:-Large-Pack-400V-nominal}

```julia
# 100 cells in series, 5 parallel = 400V, 25Ah
pack = MultiCellPack(
    name=:large_pack,
    params=OKane2022(),
    config=(100, 5),
    Qcell=5
)
# Nominal voltage: ~420V (4.2V × 100)
# Capacity: 25Ah (5Ah × 5)
```


## Cell Balancing Considerations {#Cell-Balancing-Considerations}

**Current model assumptions:**
- All cells receive equal current (ideal balancing)
  
- Temperature identical across all cells
  
- No cell-to-cell variations
  

**In reality:**
- Cell balancing requirements increase with pack size
  
- Voltage mismatches develop over cycling
  
- Temperature gradients occur in large packs
  

## Degradation in Packs {#Degradation-in-Packs}

When using parameter sets with SEI growth:

```julia
pack = MultiCellPack(
    params=Chen2020(),
    config=(12, 3),
    Qcell=5
)
# All cells experience identical SEI growth
# Capacity fade is uniform across pack
```


## Power Scaling {#Power-Scaling}

### Example: How Current Affects Simulation {#Example:-How-Current-Affects-Simulation}

```julia
pack = MultiCellPack(config=(12, 3), Qcell=5)

# 48V pack, 15Ah capacity = 720Wh energy

# Conservative discharge: 100W for ~7 hours
exp1 = Experiment([PowerStep(100, 25000)])

# High-rate discharge: 5000W for ~0.14 hours
exp2 = Experiment([PowerStep(5000, 500)])

# Both drain the same energy but at different rates
# High-rate discharge produces higher heat and degradation
```


## Next Steps {#Next-Steps}
- [SPMe Model](spme.md): Single-cell electrochemistry
  
- [Side Reactions](side-reactions.md): Degradation mechanisms
  
- [Examples](../guide/examples.md): Real scenarios
  
- [API Reference](../api/packmodels.md): Complete pack API
  
