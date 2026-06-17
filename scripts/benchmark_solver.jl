# ==============================================================================
# benchmark_solver.jl
# Compares ODE algorithms and tolerances for simulation speed and accuracy
# Evaluates base SPMe and coupled thermal SPMe
# Modes: :native, :precompiled, :overhead
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots
using Plots.Measures
using Printf
using Logging

using Revise
using BatteryToolkit

Revise.revise()

function optimise_solvers(mode=:precompiled, include_rosenbrock=false)
    # Check for valid execution mode
    if mode ∉ [:native, :precompiled, :overhead]
        error("Invalid mode. Choose :native, :precompiled, or :overhead")
    end

    # Initialise plot
    m = 7mm
    p_volt = plot(xlabel="Time (s)", ylabel="Voltage (V)", legend=:bottomleft, margin=m)

    # Setup parameters for single C discharge from full
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
    @mtkbuild sys_tspme = SingleCellCoreShellPack(params=p, config=(1,1), h_conv=30.0, T_ambient=298.15)

    # Define twenty watt continuous power discharge
    exp = Experiment([PowerStep(20.0, 3600.0)])

    # Compile problems ahead of time if bypassing overhead
    if mode == :native || mode == :precompiled
        println("Compiling ODE problems...")
        prob_spme = ODEProblem(sys_spme, [sys_spme.Pin => exp.p0, sys_spme.Iin => 0.0], (0.0, exp.tend))
        prob_tspme = ODEProblem(sys_tspme, [sys_tspme.Pin => exp.p0, sys_tspme.Iin => 0.0], (0.0, exp.tend))
    end

    # Define algorithms and tolerances to test
    algorithms = Any[
        ("QNDF", QNDF()),
        ("FBDF", FBDF()),
        ("Rodas4", Rodas4()),
        ("Rodas5P", Rodas5P())
    ]
    
    # Toggle slow solver
    if include_rosenbrock
        push!(algorithms, ("Rosenbrock23", Rosenbrock23()))
    end

    tolerances = [
        (1e-3, 1e-6),
        (1e-4, 1e-7),
        (1e-5, 1e-8)
    ]

    # Initialise metrics array
    results = []
    colors = palette(:tab10)
    plot_idx = 1

    println("\nStarting benchmark sequence using mode: $mode")

    for (alg_name, alg) in algorithms
        
        # Warm up compiler based on selected mode
        println("\nWarming up $alg_name...")
        with_logger(NullLogger()) do
            if mode == :native
                solve(prob_spme, alg, reltol=1e-3, abstol=1e-6, saveat=10.0, maxiters=1e7, verbose=false)
                solve(prob_tspme, alg, reltol=1e-3, abstol=1e-6, saveat=10.0, maxiters=1e7, verbose=false)
            elseif mode == :precompiled
                simulate(sys_spme, prob_spme, exp, alg, reltol=1e-3, abstol=1e-6, saveat=10.0, maxiters=1e7, verbose=false)
                simulate(sys_tspme, prob_tspme, exp, alg, reltol=1e-3, abstol=1e-6, saveat=10.0, maxiters=1e7, verbose=false)
            elseif mode == :overhead
                simulate(sys_spme, exp, alg, reltol=1e-3, abstol=1e-6, saveat=10.0, maxiters=1e7, verbose=false)
                simulate(sys_tspme, exp, alg, reltol=1e-3, abstol=1e-6, saveat=10.0, maxiters=1e7, verbose=false)
            end
        end

        for (rtol, atol) in tolerances
            label_tspme = "$(alg_name) TSPMe (rt: $(rtol))"
            println("Testing $alg_name at reltol=$rtol...")

            # Time SPMe execution while silencing internal package warnings
            t_start = time()
            sol_spme = with_logger(NullLogger()) do
                if mode == :native
                    solve(prob_spme, alg, reltol=rtol, abstol=atol, saveat=10.0, maxiters=1e7, verbose=false)
                elseif mode == :precompiled
                    simulate(sys_spme, prob_spme, exp, alg, reltol=rtol, abstol=atol, saveat=10.0, maxiters=1e7, verbose=false)
                elseif mode == :overhead
                    simulate(sys_spme, exp, alg, reltol=rtol, abstol=atol, saveat=10.0, maxiters=1e7, verbose=false)
                end
            end
            time_spme = time() - t_start
            v_spme = sol_spme[sys_spme.cell.v][end]

            # Time thermal SPMe execution while silencing internal package warnings
            t_start = time()
            sol_tspme = with_logger(NullLogger()) do
                if mode == :native
                    solve(prob_tspme, alg, reltol=rtol, abstol=atol, saveat=10.0, maxiters=1e7, verbose=false)
                elseif mode == :precompiled
                    simulate(sys_tspme, prob_tspme, exp, alg, reltol=rtol, abstol=atol, saveat=10.0, maxiters=1e7, verbose=false)
                elseif mode == :overhead
                    simulate(sys_tspme, exp, alg, reltol=rtol, abstol=atol, saveat=10.0, maxiters=1e7, verbose=false)
                end
            end
            time_tspme = time() - t_start
            v_tspme = sol_tspme[sys_tspme.cell.v][end]

            # Print single clean status message per condition
            println("  -> Terminated with status: $(sol_tspme.retcode)")

            # Store metrics
            push!(results, (alg_name, rtol, atol, time_spme, v_spme, time_tspme, v_tspme))

            # Add thermal SPMe curve to plot
            plot!(p_volt, sol_tspme.t, sol_tspme[sys_tspme.cell.v], label=label_tspme, lw=2, color=colors[plot_idx % 10 + 1], alpha=0.7)
            plot_idx += 1
        end
    end

    # Output text table
    println("\n" * "="^92)
    println("Solver performance comparison ($mode mode)")
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

# Safely invoke function to avoid world age errors
# Pass mode and boolean to include slow solver
Base.invokelatest(optimise_solvers, :precompiled, false)



# Overhead - Old results (with symbolics):
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

# Overhead - New results:
# ============================================================================================
# Solver performance comparison
# --------------------------------------------------------------------------------------------
# Algorithm       | reltol   | abstol   | Time SPMe(s) | V SPMe   | Time TSPMe(s)  | V TSPMe 
# --------------------------------------------------------------------------------------------
# QNDF            | 1.0e-03  | 1.0e-06  | 12.383       | 2.5000   | 11.574         | 2.5000  
# QNDF            | 1.0e-04  | 1.0e-07  | 10.909       | 2.5000   | 10.211         | 2.5000  
# QNDF            | 1.0e-05  | 1.0e-08  | 10.168       | 2.5000   | 10.553         | 2.5000
# FBDF            | 1.0e-03  | 1.0e-06  | 9.728        | 2.5000   | 10.091         | 2.5000
# FBDF            | 1.0e-04  | 1.0e-07  | 9.592        | 2.5000   | 9.948          | 2.5000
# FBDF            | 1.0e-05  | 1.0e-08  | 10.643       | 2.5000   | 9.803          | 2.5000
# Rodas4          | 1.0e-03  | 1.0e-06  | 9.507        | 2.5000   | 10.104         | 2.5000
# Rodas4          | 1.0e-04  | 1.0e-07  | 9.700        | 2.5000   | 9.981          | 2.5000
# Rodas4          | 1.0e-05  | 1.0e-08  | 9.628        | 2.5000   | 10.065         | 2.5000
# Rodas5P         | 1.0e-03  | 1.0e-06  | 9.635        | 2.5000   | 10.287         | 2.5000
# Rodas5P         | 1.0e-04  | 1.0e-07  | 9.879        | 2.5000   | 10.246         | 2.5000
# Rodas5P         | 1.0e-05  | 1.0e-08  | 10.607       | 2.5000   | 10.278         | 2.5000
# --------------------------------------------------------------------------------------------
# ============================================================================================




# Native - Old results (with symbolics):
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

# Native - New results:
# ============================================================================================
# Native solver performance comparison (no overhead)
# --------------------------------------------------------------------------------------------
# Algorithm       | reltol   | abstol   | Time SPMe(s) | V SPMe   | Time TSPMe(s)  | V TSPMe 
# --------------------------------------------------------------------------------------------
# QNDF            | 1.0e-03  | 1.0e-06  | 0.005        | 2.5000   | 0.006          | 2.5000  
# QNDF            | 1.0e-04  | 1.0e-07  | 0.006        | 2.5000   | 0.008          | 2.5000  
# QNDF            | 1.0e-05  | 1.0e-08  | 0.009        | 2.5000   | 0.011          | 2.5000
# FBDF            | 1.0e-03  | 1.0e-06  | 0.008        | 2.5000   | 0.099          | 2.5000
# FBDF            | 1.0e-04  | 1.0e-07  | 0.010        | 2.5000   | 0.012          | 2.5000
# FBDF            | 1.0e-05  | 1.0e-08  | 0.013        | 2.5000   | 0.015          | 2.5000
# Rodas4          | 1.0e-03  | 1.0e-06  | 0.013        | 2.5000   | 0.017          | 2.5000
# Rodas4          | 1.0e-04  | 1.0e-07  | 0.019        | 2.5000   | 0.024          | 2.5000
# Rodas4          | 1.0e-05  | 1.0e-08  | 0.030        | 2.5000   | 0.034          | 2.5000
# Rodas5P         | 1.0e-03  | 1.0e-06  | 0.013        | 2.5000   | 0.015          | 2.5000
# Rodas5P         | 1.0e-04  | 1.0e-07  | 0.018        | 2.5000   | 0.021          | 2.5000
# Rodas5P         | 1.0e-05  | 1.0e-08  | 0.022        | 2.5000   | 0.085          | 2.5000
# Rosenbrock23    | 1.0e-03  | 1.0e-06  | 54.532       | 2.5000   | 63.252         | 2.5000
# Rosenbrock23    | 1.0e-04  | 1.0e-07  | 180.059      | 2.5000   | 200.781        | 2.5000
# Rosenbrock23    | 1.0e-05  | 1.0e-08  | 280.586      | 3.2621   | 331.222        | 3.2426
# --------------------------------------------------------------------------------------------
# ============================================================================================



# ============================================================================================
# Solver performance comparison (precompiled mode)
# --------------------------------------------------------------------------------------------
# Algorithm       | reltol   | abstol   | Time SPMe(s) | V SPMe   | Time TSPMe(s)  | V TSPMe
# --------------------------------------------------------------------------------------------
# QNDF            | 1.0e-03  | 1.0e-06  | 0.031        | 2.5000   | 0.030          | 2.5000
# QNDF            | 1.0e-04  | 1.0e-07  | 0.007        | 2.5000   | 0.009          | 2.5000
# QNDF            | 1.0e-05  | 1.0e-08  | 0.009        | 2.5000   | 0.011          | 2.5000
# FBDF            | 1.0e-03  | 1.0e-06  | 0.034        | 2.5000   | 0.033          | 2.5000
# FBDF            | 1.0e-04  | 1.0e-07  | 0.009        | 2.5000   | 0.013          | 2.5000
# FBDF            | 1.0e-05  | 1.0e-08  | 0.014        | 2.5000   | 0.014          | 2.5000
# Rodas4          | 1.0e-03  | 1.0e-06  | 0.038        | 2.5000   | 0.037          | 2.5000
# Rodas4          | 1.0e-04  | 1.0e-07  | 0.020        | 2.5000   | 0.025          | 2.5000
# Rodas4          | 1.0e-05  | 1.0e-08  | 0.031        | 2.5000   | 0.037          | 2.5000
# Rodas5P         | 1.0e-03  | 1.0e-06  | 0.036        | 2.5000   | 0.036          | 2.5000  
# Rodas5P         | 1.0e-04  | 1.0e-07  | 0.018        | 2.5000   | 0.021          | 2.5000
# Rodas5P         | 1.0e-05  | 1.0e-08  | 0.023        | 2.5000   | 0.028          | 2.5000
# --------------------------------------------------------------------------------------------
# ============================================================================================