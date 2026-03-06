# ==============================================================================
# benchmark_solver.jl
# Compares ODE algorithms and tolerances for simulation speed and accuracy
# Evaluates both base SPMe and coupled core-shell thermal SPMe
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots
using Plots.Measures
using Printf

using Revise
using BatteryPeck

Revise.revise()

function optimise_solvers()
    # Initialise plot
    m = 7mm
    p_volt = plot(xlabel="Time (s)", ylabel="Voltage (V)", legend=:bottomleft, margin=m)

    # Setup parameters for 1C discharge from full
    p = Chen2020()
    soc_init = 1.0
    p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
    p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
    
    # Set voltage limits to stop naturally at empty
    p.Vmin = 2.5
    p.Vmax = 4.3

    # Build systems once to avoid recompilation overhead
    println("Building systems...")
    @mtkbuild sys_spme = SingleCellPack(params=p, config=(1,1))
    @mtkbuild sys_cs = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=30.0, T_ambient=298.15)

    # Define 20W continuous power discharge
    exp = Experiment([PowerStep(20.0, 3600.0)])

    # Define algorithms and tolerances to test
    algorithms = [
        ("QNDF", QNDF()),
        ("FBDF", FBDF()),
        ("Rodas4", Rodas4()),
        ("Rodas5P", Rodas5P()),
        ("Rosenbrock23", Rosenbrock23())
    ]

    tolerances = [
        (1e-3, 1e-6),
        (1e-4, 1e-7),
        (1e-5, 1e-8)
    ]

    # Initialise metrics array
    results = []
    colors = palette(:tab10)
    plot_idx = 1

    println("\nStarting benchmark sequence...")

    for (alg_name, alg) in algorithms
        
        # Warm up compiler for both systems
        println("\nWarming up $alg_name...")
        simulate(sys_spme, exp, alg, reltol=1e-3, abstol=1e-6, saveat=10.0)
        simulate(sys_cs, exp, alg, reltol=1e-3, abstol=1e-6, saveat=10.0)

        for (rtol, atol) in tolerances
            label_cs = "$(alg_name) CS (rt: $(rtol))"
            println("Testing $alg_name at reltol=$rtol...")

            # Time SPMe execution
            t_start = time()
            sol_spme = simulate(sys_spme, exp, alg, reltol=rtol, abstol=atol, saveat=10.0)
            time_spme = time() - t_start
            v_spme = sol_spme[sys_spme.cell.v][end]

            # Time Core-Shell execution
            t_start = time()
            sol_cs = simulate(sys_cs, exp, alg, reltol=rtol, abstol=atol, saveat=10.0)
            time_cs = time() - t_start
            v_cs = sol_cs[sys_cs.cell.v][end]

            # Store metrics
            push!(results, (alg_name, rtol, atol, time_spme, v_spme, time_cs, v_cs))

            # Add core-shell curve to plot
            plot!(p_volt, sol_cs.t, sol_cs[sys_cs.cell.v], label=label_cs, lw=2, color=colors[plot_idx % 10 + 1], alpha=0.7)
            plot_idx += 1
        end
    end

    # Output text table
    println("\n" * "="^92)
    println("Solver Performance Comparison")
    println("-"^92)
    @printf("%-15s | %-8s | %-8s | %-12s | %-8s | %-14s | %-8s\n", 
            "Algorithm", "reltol", "abstol", "Time SPMe(s)", "V SPMe", "Time TSPMe(s)", "V TSPMe")
    println("-"^92)

    for res in results
        @printf("%-15s | %-8.1e | %-8.1e | %-12.3f | %-8.4f | %-14.3f | %-8.4f\n", 
                res[1], res[2], res[3], res[4], res[5], res[6], res[7])
    end
    
    println("-"^92)
    println("="^92 * "\n")

    # Display plot
    display(p_volt)
end

optimise_solvers()

# Results:
# ==========================================================================================
# Solver Performance Comparison
# ------------------------------------------------------------------------------------------
# Algorithm       | reltol   | abstol   | Time SPMe(s) | V SPMe   | Time TSPMe(s)| V TSPMe    
# ------------------------------------------------------------------------------------------
# QNDF            | 1.0e-03  | 1.0e-06  | 2.080        | 2.5000   | 2.215        | 2.5000  
# QNDF            | 1.0e-04  | 1.0e-07  | 2.089        | 2.5000   | 2.134        | 2.5000  
# QNDF            | 1.0e-05  | 1.0e-08  | 2.154        | 2.5000   | 2.243        | 2.5000
# FBDF            | 1.0e-03  | 1.0e-06  | 2.227        | 2.5000   | 2.297        | 2.5000
# FBDF            | 1.0e-04  | 1.0e-07  | 2.173        | 2.5000   | 2.231        | 2.5000
# FBDF            | 1.0e-05  | 1.0e-08  | 2.242        | 2.5000   | 2.274        | 2.5000
# Rodas4          | 1.0e-03  | 1.0e-06  | 2.140        | 2.5000   | 2.361        | 2.5000
# Rodas4          | 1.0e-04  | 1.0e-07  | 3.039        | 2.5000   | 2.268        | 2.5000
# Rodas4          | 1.0e-05  | 1.0e-08  | 2.358        | 2.5000   | 2.391        | 2.5000
# Rodas5P         | 1.0e-03  | 1.0e-06  | 2.200        | 2.5000   | 2.251        | 2.5000
# Rodas5P         | 1.0e-04  | 1.0e-07  | 2.185        | 2.5000   | 2.329        | 2.5000
# Rodas5P         | 1.0e-05  | 1.0e-08  | 2.250        | 2.5000   | 2.399        | 2.5000
# Rosenbrock23    | 1.0e-03  | 1.0e-06  | 2.948        | 2.5000   | 2.815        | 2.5000
# Rosenbrock23    | 1.0e-04  | 1.0e-07  | 4.666        | 2.5000   | 3.435        | 2.5000
# Rosenbrock23    | 1.0e-05  | 1.0e-08  | 10.573       | 2.5000   | 5.423        | 2.5000
# ------------------------------------------------------------------------------------------
# ==========================================================================================