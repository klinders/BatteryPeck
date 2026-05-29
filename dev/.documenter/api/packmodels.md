
# Pack Models API Reference {#Pack-Models-API-Reference}

API documentation for single-cell and multi-cell battery packs.

## Pack Constructors {#Pack-Constructors}

### Single Cell Pack {#Single-Cell-Pack}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.SingleCellPack' href='#BatteryToolkit.SingleCellPack'><span class="jlbinding">BatteryToolkit.SingleCellPack</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
SingleCellPack(; name, params=Chen2020(), config=(96,3), Qcell=5)
```


Create a single battery cell pack with configurable series/parallel arrangement.

Builds a ModelingToolkit system representing a battery pack of series-parallel connected cells. Each cell is modeled using the SPMe (Single Particle Model with Electrolyte) framework.

**Keyword Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `params::BatteryParameters`: Parameter set for the cell model (default: Chen2020 parameters)
  
- `config::Tuple{Int,Int}`: Series and parallel configuration, `(n_series, n_parallel)` (default: 96 series, 3 parallel)
  
- `Qcell::Float64`: Nominal cell capacity in Ah (default: 5 Ah)
  

**Returns**
- ModelingToolkit system with inputs (Pin, Iin, Tin) and outputs (V, I)
  

**Input/Output Signals**
- Inputs: `Pin` (power), `Iin` (current), `Tin` (temperature)
  
- Outputs: `V` (voltage), `I` (current)
  

**Example**

```julia
params = OKane2022()
pack = SingleCellPack(name=:pack, params=params, config=(96,3), Qcell=5)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/PackModels/SingleCellPack.jl#L5-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Multi-Cell Pack {#Multi-Cell-Pack}
<details class='jldocstring custom-block' open>
<summary><a id='BatteryToolkit.MultiCellPack' href='#BatteryToolkit.MultiCellPack'><span class="jlbinding">BatteryToolkit.MultiCellPack</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
MultiCellPack(; name, params=Chen2020(), config=(12,1), Qcell=5)
```


Create a multi-cell battery pack with flexible series/parallel configuration.

Builds a ModelingToolkit system representing a large battery pack with multiple cells in series and parallel. Each cell is modeled using SPMe, and they are electrically connected according to the specified configuration.

**Keyword Arguments**
- `name`: System name for ModelingToolkit (required)
  
- `params::BatteryParameters`: Parameter set for all cells (default: Chen2020 parameters)
  
- `config::Tuple{Int,Int}`: Series and parallel configuration, `(n_series, n_parallel)` (default: 12 series, 1 parallel)
  
- `Qcell::Float64`: Nominal cell capacity in Ah (default: 5 Ah)
  

**Returns**
- ModelingToolkit system with inputs (P, T) and outputs (V, I) for the entire pack
  

**Input/Output Signals**
- Inputs: `P` (power), `T` (temperature)
  
- Outputs: `V` (pack voltage = sum of series cell voltages), `I` (pack current)
  

**Notes**
- All cells share the same temperature input
  
- Total pack voltage is the sum of series-connected cell voltages
  
- Parallel cells share equal current
  

**Example**

```julia
params = Chen2020()
pack = MultiCellPack(name=:battery, params=params, config=(96,3), Qcell=5)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/klinders/BatteryToolkit/blob/b4d229b7acb19e134fa7d4efec0646d85d1a6085/src/PackModels/MultiCellPack.jl#L5-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

