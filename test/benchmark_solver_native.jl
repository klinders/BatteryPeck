# ==============================================================================
# benchmark_solver_native.jl
# Compares ODE algorithms and tolerances using pure mathematical integration
# Bypasses experimental control loop overhead
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

function optimise_solvers_native()
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

    println("\nStarting native benchmark sequence...")

    for (alg_name, alg) in algorithms
        
        # Build ODE problems with constant -20W load and sparse AutoDiff
        prob_spme = ODEProblem(sys_spme, [sys_spme.P => -20.0], (0.0, 3600.0), sparse=true)
        prob_cs = ODEProblem(sys_cs, [sys_cs.P => -20.0], (0.0, 3600.0), sparse=true)

        # Warm up compiler for both systems
        println("\nWarming up $alg_name...")
        solve(prob_spme, alg, reltol=1e-3, abstol=1e-6, saveat=10.0)
        solve(prob_cs, alg, reltol=1e-3, abstol=1e-6, saveat=10.0)

        for (rtol, atol) in tolerances
            label_cs = "$(alg_name) CS (rt: $(rtol))"
            println("Testing $alg_name at reltol=$rtol...")

            # Time SPMe execution
            t_start = time()
            sol_spme = solve(prob_spme, alg, reltol=rtol, abstol=atol, saveat=10.0)
            time_spme = time() - t_start
            v_spme = sol_spme[sys_spme.cell.v][end]

            # Time Core-Shell execution
            t_start = time()
            sol_cs = solve(prob_cs, alg, reltol=rtol, abstol=atol, saveat=10.0)
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
    println("Native Solver Performance Comparison (No Overhead)")
    println("-"^92)
    
    # Increased the sixth column width from 12 to 14 to fit "Time TSPMe(s)"
    @printf("%-15s | %-8s | %-8s | %-12s | %-8s | %-14s | %-8s\n", 
            "Algorithm", "reltol", "abstol", "Time SPMe(s)", "V SPMe", "Time TSPMe(s)", "V TSPMe")
    println("-"^92)

    for res in results
        # Matched the sixth column width here (14.3f)
        @printf("%-15s | %-8.1e | %-8.1e | %-12.3f | %-8.4f | %-14.3f | %-8.4f\n", 
                res[1], res[2], res[3], res[4], res[5], res[6], res[7])
    end
    
    println("-"^92)
    println("="^92 * "\n")

    # Display plot
    display(p_volt)
end

optimise_solvers_native()

# Results:
# ==========================================================================================
# Native Solver Performance Comparison (No Overhead)
# ------------------------------------------------------------------------------------------
# Algorithm       | reltol   | abstol   | Time SPMe(s) | V SPMe   | Time TSPMe(s)| V TSPMe
# ------------------------------------------------------------------------------------------
# QNDF            | 1.0e-03  | 1.0e-06  | 0.013        | 2.5000   | 0.014        | 2.5000
# QNDF            | 1.0e-04  | 1.0e-07  | 0.018        | 2.5000   | 0.021        | 2.5000
# QNDF            | 1.0e-05  | 1.0e-08  | 0.026        | 2.5000   | 0.107        | 2.5000
# FBDF            | 1.0e-03  | 1.0e-06  | 0.017        | 2.5000   | 0.018        | 2.5000
# FBDF            | 1.0e-04  | 1.0e-07  | 0.022        | 2.5000   | 0.024        | 2.5000
# FBDF            | 1.0e-05  | 1.0e-08  | 0.031        | 2.5000   | 0.081        | 2.5000
# Rodas4          | 1.0e-03  | 1.0e-06  | 0.098        | 2.5000   | 0.075        | 2.5000
# Rodas4          | 1.0e-04  | 1.0e-07  | 0.055        | 2.5000   | 0.086        | 2.5000
# Rodas4          | 1.0e-05  | 1.0e-08  | 0.112        | 2.5000   | 0.119        | 2.5000
# Rodas5P         | 1.0e-03  | 1.0e-06  | 0.217        | 2.5000   | 0.132        | 2.5000
# Rodas5P         | 1.0e-04  | 1.0e-07  | 0.981        | 2.5000   | 5.605        | 2.5000
# Rodas5P         | 1.0e-05  | 1.0e-08  | 0.276        | 2.5000   | 0.292        | 2.5000
# Rosenbrock23    | 1.0e-03  | 1.0e-06  | 0.327        | 2.5000   | 0.362        | 2.5000
# Rosenbrock23    | 1.0e-04  | 1.0e-07  | 0.881        | 2.5000   | 0.954        | 2.5000
# Rosenbrock23    | 1.0e-05  | 1.0e-08  | 2.782        | 2.5000   | 2.916        | 2.5000
# ------------------------------------------------------------------------------------------
# ==========================================================================================