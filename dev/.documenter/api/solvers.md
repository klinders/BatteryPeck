
# Solvers API Reference {#Solvers-API-Reference}

API documentation for simulation solvers and utilities.

## Simulation {#Simulation}

### Main Solver {#Main-Solver}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.simulate' href='#BatteryToolkit.simulate'><span class="jlbinding">BatteryToolkit.simulate</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
simulate(sys::AbstractSystem, experiment::Experiment, args...; kwargs...)
```


Run a battery simulation using the provided ModelingToolkit system and experiment profile.

Simulates the battery system through all steps defined in the experiment, updating the solver after each step and checking safety limits. Returns the solution object with time series data.

**Arguments**
- `sys::AbstractSystem`: ModelingToolkit system (typically from `SPMe()` or pack models)
  
- `experiment::Experiment`: Experiment profile containing steps and parameters
  
- `args...`: Positional arguments passed to ODE solver (e.g., solver algorithm)
  
- `kwargs...`: Keyword arguments passed to ODE solver (e.g., `abstol`, `reltol`)
  

**Returns**
- Solution object with fields `t` (time) and `u` (state variables) for plotting/analysis
  

**Example**

```julia
params = OKane2022()
sys = SPMe(params=params)
exp = Experiment([PowerStep(1000, 3600)])  # 1000W for 1 hour
sol = simulate(sys, exp, Rodas4())
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/solvers.jl#L4-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>

