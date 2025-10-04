using BatteryPeck

using ModelingToolkit, DifferentialEquations, ModelingToolkitStandardLibrary.Blocks, ModelingToolkitStandardLibrary.Electrical
using ModelingToolkit: t_nounits as t

params = Chen2020()

Ncell = 10

function PackF(;name)
    @named current = Constant(k=-5)
    @named source = Current() # Current source
    
    battery = [SPMe(name=Symbol("battery_$i"), params=params) for i in 1:Ncell] # Battery model
    
    eqns = [
        connect(current.output, source.I)
        connect(source.p, battery[1].p)
        [connect(battery[i].n, battery[i+1].p) for i in 1:Ncell-1]...
        connect(source.n, battery[end].n)
    ]

    return System(eqns,t; name=name, systems=[current,source,battery...])
end

function build()
    print("Start building the system...\n")

    @mtkbuild pde_sys = PackF()

    print("Start building the problem...\n")

    prob = ODEProblem(pde_sys, [], (0,500); sparse=true)
   return prob
end

@time prob = build()

print("Start solving the problem...\n")

@time sol = solve(prob, Rodas5(), abstol=1e-6, reltol=1e-6);

print("Done\n")