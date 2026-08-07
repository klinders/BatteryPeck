# Parameters API Reference

Complete API documentation for battery parameter types and constructors.

## Parameter Constructors

```@docs
Chen2020
OKane2022
```

## Parameter Types

```@autodocs
Modules = [BatteryToolkit]
Pages = ["ParameterSets/Base.jl"]
Order = [:type, :function]
Filter = t -> typeof(t) <: Type || typeof(t) <: Function
```
