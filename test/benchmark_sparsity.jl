# ==============================================================================
# benchmark_sparsity.jl
# Evaluates impact of sparsity patterns on solver speed using AutoDiff
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

function benchmark_jacobian()
    # Setup parameters for 1C continuous discharge
    p = Chen2020()
    soc_init = 1.0
    p.n.c₀ = (p.n.z_0 + soc_init * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
    p.p.c₀ = (p.p.z_0 + soc_init * (p.p.z_100 - p.p.z_0)) * p.p.c₊
    
    # Set voltage limits
    p.Vmin = 2.5
    p.Vmax = 4.3

    println("Building system...")
    @mtkbuild sys = SingleCellPack(params=p, config=(1,1))

    # Define algorithm and tolerances based on previous optimisation
    alg = QNDF()
    rtol = 1e-4
    atol = 1e-7

    # Define Jacobian test configurations
    # jac=true is disabled due to array symbolics limitations
    configs = [
        ("Dense AutoDiff (Default)", false, false),
        ("Sparse AutoDiff (MTK pattern)", false, true)
    ]

    results = []
    
    println("\nStarting Jacobian benchmark sequence...")

    for (label, use_jac, use_sparse) in configs
        println("\nEvaluating: $label...")
        
        # Build ODE problem manually to pass jacobian arguments directly
        prob = ODEProblem(sys, [sys.P => -20.0], (0.0, 3600.0), jac=use_jac, sparse=use_sparse)
        
        # Warm up compiler for specific configuration
        println("  Warming up...")
        solve(prob, alg, reltol=1e-3, abstol=1e-6, saveat=10.0)

        # Time execution
        println("  Running timed execution...")
        t_start = time()
        sol = solve(prob, alg, reltol=rtol, abstol=atol, saveat=10.0)
        solve_time = time() - t_start

        final_v = sol[sys.cell.v][end]
        push!(results, (label, solve_time, final_v))
    end

    # Output text table
    println("\n" * "="^65)
    println("Sparsity performance comparison")
    println("-"^65)
    @printf("%-30s | %-15s | %-10s\n", "Configuration", "Time (s)", "Final V")
    println("-"^65)

    for res in results
        @printf("%-30s | %-15.3f | %-10.4f\n", res[1], res[2], res[3])
    end
    
    println("-"^65)
    println("="^65 * "\n")
end

benchmark_jacobian()

# Results:
# =================================================================
# Sparsity performance comparison
# -----------------------------------------------------------------
# Configuration                  | Time (s)        | Final V
# -----------------------------------------------------------------
# Dense AutoDiff (Default)       | 0.146           | 2.5000
# Sparse AutoDiff (MTK pattern)  | 0.063           | 2.5000
# -----------------------------------------------------------------
# =================================================================