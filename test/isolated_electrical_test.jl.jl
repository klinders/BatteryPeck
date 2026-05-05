# ==============================================================================
# isolated_electrical_test.jl
# Diagnostic script to isolate SPMe series electrical wiring
# ==============================================================================

using ModelingToolkit
using OrdinaryDiffEq
using Plots
using Revise
using Logging
using BatteryToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks

Revise.revise()

"""
    run_test()

Execute purely electrical multi-cell series diagnostic.
"""
function run_test()
    println("Initialising pure electrical diagnostic...")

    # Load base parameters and loosen constraints
    p = Chen2020()
    p.Vmin = 1.0 
    p.Vmax = 5.5

    @parameters t
    
    # Instantiate bare SPMe cells
    @named cell1 = SPMe(params=p)
    @named cell2 = SPMe(params=p)
    @named cell3 = SPMe(params=p)
    @named cell4 = SPMe(params=p)
    
    # Instantiate electrical boundaries
    @named source = Current()
    @named ground = Ground()
    
    @variables V(t)
    
    # Assemble pure series circuit equations
    eqs = [
        # Define pack voltage
        V ~ cell1.p.v - cell4.n.v
        
        # Apply constant 1.25A discharge
        source.I.u ~ 1.25
        
        # Provide standard 25C temperature to the bare cells
        cell1.T.u ~ 298.15
        cell2.T.u ~ 298.15
        cell3.T.u ~ 298.15
        cell4.T.u ~ 298.15
        
        # Wire components in series
        connect(source.p, cell1.p)
        connect(cell1.n, cell2.p)
        connect(cell2.n, cell3.p)
        connect(cell3.n, cell4.p)
        connect(cell4.n, source.n)
        
        # Ground the negative terminal of the entire string
        connect(source.n, ground.g)
    ]
    
    # Construct uncoupled ODE system
    @named sys = ODESystem(eqs, t, [V], []; systems=[cell1, cell2, cell3, cell4, source, ground])
    
    println("Simplifying system...")
    sys_simp = structural_simplify(sys)
    
    # Provide the staggered voltage guesses we derived earlier to help the series initialisation
    my_guesses = Dict()
    my_guesses[sys_simp.V] = 16.76
    
    for (i, c) in enumerate([sys_simp.cell1, sys_simp.cell2, sys_simp.cell3, sys_simp.cell4])
        offset = (4 - i) * 4.19
        my_guesses[c.v] = 4.19
        my_guesses[c.n.v] = offset
        my_guesses[c.p.v] = offset + 4.19
    end
    
    println("Compiling ODE problem...")
    prob = ODEProblem(sys_simp, [], (0.0, 3600.0); guesses=my_guesses, sparse=true, jac=false)
    
    println("Running simulation...")
    # Solve directly using standard QNDF
    sol = solve(prob, QNDF(autodiff=false); saveat=1.0)
    
    # Plot pack voltage
    p_volt = plot(sol.t, sol[sys_simp.V], label="Pack voltage", xlabel="Time [s]", ylabel="Voltage [V]", title="Pure electrical series check", lw=2)
    display(p_volt)
    
    println("Pure electrical test complete!")
end

Base.invokelatest(run_test)