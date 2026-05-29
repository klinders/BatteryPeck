
# Experiments API Reference {#Experiments-API-Reference}

API documentation for experiment types and control functions.

## Experiment Type {#Experiment-Type}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.Experiment' href='#BatteryToolkit.Experiment'><span class="jlbinding">BatteryToolkit.Experiment</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Experiment(steps::Vector{<:AbstractStep}, start_time::DateTime=DateTime(2020, 1, 1))
```


Create an experimental profile composing multiple battery operation steps.

Combines a sequence of operation steps (power, current, rest, charge, drive cycle) into a single experiment. Automatically calculates step timing and prepares parameters for simulation with the `simulate()` function.

**Arguments**
- `steps::Vector{<:AbstractStep}`: Vector of step objects (RestStep, PowerStep, CurrentStep, ChargeStep, DriveStep)
  
- `start_time::DateTime`: Real-world timestamp for first step (default: 2020-01-01)
  

**Fields (automatically calculated)**
- `steps::Vector`: Original step vector
  
- `tstops::Vector{Float64}`: Cumulative time at end of each step except the last (s)
  
- `tend::Float64}`: Total simulation duration (s)
  
- `step_count::Int64`: Number of steps
  
- `p0::Float64`: Initial power/current value (W or A)
  
- `start_time::DateTime`: Experiment start timestamp
  

**Example**

```julia
steps = [
    PowerStep(1000, 1800),      # 1000W for 30 min
    RestStep(300),               # 5 min rest
    PowerStep(-500, 3600)        # -500W (discharge) for 1 hour
]
exp = Experiment(steps)
sol = simulate(sys, exp, Rodas4())
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L101-L132" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Step Types {#Step-Types}

### AbstractStep Hierarchy {#AbstractStep-Hierarchy}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.AbstractStep' href='#BatteryToolkit.AbstractStep'><span class="jlbinding">BatteryToolkit.AbstractStep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for all step types in the experiment.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L6-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Step Implementations {#Step-Implementations}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.RestStep' href='#BatteryToolkit.RestStep'><span class="jlbinding">BatteryToolkit.RestStep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Rest for a period


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L11-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.ChargeStep' href='#BatteryToolkit.ChargeStep'><span class="jlbinding">BatteryToolkit.ChargeStep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Charge up to specified SoC using the given power for the given period.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L18-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.PowerStep' href='#BatteryToolkit.PowerStep'><span class="jlbinding">BatteryToolkit.PowerStep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Apply a given  `power` for a given `period`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L28-L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.CurrentStep' href='#BatteryToolkit.CurrentStep'><span class="jlbinding">BatteryToolkit.CurrentStep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Apply a given  `current` for a given `period`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L36-L38" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.DriveStep' href='#BatteryToolkit.DriveStep'><span class="jlbinding">BatteryToolkit.DriveStep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Apply a drivecycle from the given csv

**Arguments**
- `csv ::Vector{Any}` Path to the csv driving cycle
  
- `period ::Real` (optional) time to apply the cycle in seconds
  

At the moment, the period can be up to the lenght of the csv. 

TODO: Repeat the cycle when period is longer then csv


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L44-L54" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Control Functions {#Control-Functions}

### Step Execution {#Step-Execution}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.step!' href='#BatteryToolkit.step!'><span class="jlbinding">BatteryToolkit.step!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
step!(integrator::SciMLBase.DEIntegrator, sys::ModelingToolkit.AbstractSystem, step::AbstractStep)
```


Apply a single step of the experiment to the DE integrator. This function modifies the integrator's input parameters according to the step type and advances the simulation by the step's period. The behavior depends on the specific step type (PowerStep, RestStep, ChargeStep, CurrentStep, DriveStep).

**Arguments**
- `integrator::SciMLBase.DEIntegrator`: The DE integrator to modify and step
  
- `sys::ModelingToolkit.AbstractSystem`: The system being simulated, used to access input variables
  
- `step::AbstractStep`: The step to apply, which determines how the integrator's inputs are modified and how long to step the simulation
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L159-L170" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Helper Functions {#Helper-Functions}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.get_p0' href='#BatteryToolkit.get_p0'><span class="jlbinding">BatteryToolkit.get_p0</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Get the initial value for the first step of the experiment. This is used to set the initial conditions for the simulation. The behavior depends on the specific step type (PowerStep, RestStep, ChargeStep, CurrentStep, DriveStep).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/experiment.jl#L72-L74" target="_blank" rel="noreferrer">source</a></Badge>

</details>

