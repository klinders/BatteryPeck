# =====================================================================================================================
# speed_run.jl
#
# Benchmarking file measuring time to:
# - Build; generate equations
# - First solve; compile
# - Second solve; execute
# =====================================================================================================================

# Import package manager
using Pkg
# Activate Julia environment in current directory
Pkg.activate(".")

# Import module
using BatteryToolkit

# Import packages
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks, ModelingToolkitStandardLibrary.Electrical

# Select battery parameter set
params = Chen2020()

# Number of cells to be simulated
Ns = [1, 10, 100]

# Create new component model
@mtkmodel Pack begin

    # Initialise loop parameter for number of cells
    @structural_parameters begin
        Ncell = 1
    end

    # Create component blocks
    @components begin
        # Constant signal for current source
        current = Constant(k=-5)
        # Ideal current source
        source = Current()
        # Reference node
        ground = Ground()
        # Constant signal for ambient temperature
        temp_source = Constant(k=298.15)
        # Create array of battery models
        battery = [SPMe(name=Symbol("battery_$i"), params=params) for i in 1:Ncell]
    end

    # Define component connections
    @equations begin
        # Signal constant value > current source
        connect(current.output, source.I)
        # Source+ > battery1+
        connect(source.p, battery[1].p)
        # Complete series loop of battery connections
        if Ncell > 1
            [connect(battery[i].n, battery[i+1].p) for i in 1:Ncell-1]...
        end
        # Source- > batteryN-
        connect(source.n, battery[end].n)
        # Source- > ground
        connect(source.n, ground.g)
        # Initialise all battery temperatures at ambient
        [connect(temp_source.output, battery[i].T) for i in 1:Ncell]...
    end
end

# Print time to [build; first solve; second solve] for every defined cell configuration
for Ncell in Ns
    print("Running for $Ncell cells...\n")

    print("Building system...\n")
    @time @mtkbuild pde_sys = Pack(Ncell=Ncell)

    print("Building problem...\n")
    @time prob = ODEProblem(pde_sys, [], (0,500); sparse=true, jac=false)

    print("First solve...\n")
    @time sol = solve(prob, TRBDF2(), abstol=1e-6, reltol=1e-6);

    print("Second solve...\n")
    @time sol = solve(prob, TRBDF2(), abstol=1e-6, reltol=1e-6);

    print("\n\n")
end

# End of test
print("Done!\n")