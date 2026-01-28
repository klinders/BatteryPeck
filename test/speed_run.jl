using Pkg
Pkg.activate(".")

using BatteryPeck

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks, ModelingToolkitStandardLibrary.Electrical

params = Chen2020()

Ns = [1, 10, 100]

@mtkmodel Pack begin

    @structural_parameters begin
        Ncell = 1
    end

    @components begin
        current = Constant(k=-5)
        temperature = Constant(k=298)
        source = Current() # Current source
        ground = Ground()
        cell = [SPMe(name=Symbol("battery_$i"), params=p) for i in 1:Ncell] # Battery model
    end

    @equations begin
        connect(current.output, source.I)
        connect(source.p, cell[1].p)
        if Ncell > 1
            [connect(cell[i].n, cell[i+1].p) for i in 1:Ncell-1]...
        end
        connect(source.n, cell[end].n)
        connect(ground.g, source.n)
        [connect(temperature.output, cell[i].T) for i in 1:Ncell]...
    end
end

for Ncell in Ns
    print("Running for $Ncell cells...\n")
    

    print("Start building the system...\n")

    @time @mtkbuild pde_sys = Pack(Ncell=Ncell)

    print("Start building the problem...\n")

    @time prob = ODEProblem(pde_sys, [], (0,500); sparse=true, jac=false)

    print("Start solving the problem...\n")

    @time sol = solve(prob, TRBDF2(), abstol=1e-6, reltol=1e-6);

    print("Start solving second time...\n")

    @time sol = solve(prob, TRBDF2(), abstol=1e-6, reltol=1e-6);

    print("\n\n")
end


print("Done\n")