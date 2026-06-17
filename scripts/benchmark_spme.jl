# ==============================================================================
# benchmark_spme.jl
# Radial nodes & solver tolerances Pareto sweep
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using BatteryToolkit
using Plots
using DataFrames
using CSV
using LinearAlgebra

function build_benchmark_cell(elec_params, nodes::Int)
    N_dict = Dict(:Nₓ => [10, 10, 10], :Nᵣ => [nodes, nodes])
    @named cell = SPMe(params=elec_params, N=N_dict, side_reactions=true)
    @named load = Current(); @named ground = Ground()
    @parameters t I_app=0.0 T_ext=298.15
    
    eqs = [
        cell.T.u ~ T_ext,
        load.I.u ~ I_app,
        connect(load.p, cell.p),
        connect(cell.n, load.n, ground.g)
    ]
    
    sys = ODESystem(eqs, t, [], [I_app, T_ext]; systems=[cell, load, ground], name=:cell_template)
    return structural_simplify(sys)
end

function simulate_single_cell(sys, total_t, get_load_func, tol, is_truth=false)
    prob = ODEProblem(sys, [sys.I_app => 0.0, sys.T_ext => 298.15], (0.0, total_t); jac=true, sparse=true)
    
    integrator = init(prob, QNDF(autodiff=true); reltol=tol, abstol=tol/10.0, save_everystep=false)
    
    hist_t = Float64[]
    hist_v = Float64[]
    
    current_t = 0.0
    dt_sync = 1.0
    
    start_wall = time()
    last_print = time()
    
    while current_t < total_t
        I_val, is_done = get_load_func(current_t, 4.0, 0.5)
        if is_done break end
        
        integrator.ps[sys.I_app] = I_val
        
        SciMLBase.u_modified!(integrator, true)
        
        t_target = current_t + dt_sync
        SciMLBase.add_tstop!(integrator, t_target)
        
        substeps = 0
        while integrator.t < t_target
            SciMLBase.step!(integrator)
            substeps += 1
            if substeps > 1000 || integrator.sol.retcode == SciMLBase.ReturnCode.Terminated 
                break 
            end
        end
        
        if integrator.sol.retcode == SciMLBase.ReturnCode.Terminated break end
        
        current_t = integrator.t
        push!(hist_t, current_t)
        push!(hist_v, integrator[sys.cell.v])
        
        if time() - last_print >= 0.5
            pct = clamp((current_t / total_t) * 100.0, 0.0, 100.0)
            filled = round(Int, 30 * (pct / 100.0))
            bar = "[" * repeat("=", filled) * repeat(" ", 30 - filled) * "]"
            v_val = integrator[sys.cell.v]
            prefix = is_truth ? "[Truth] " : "[Test]  "
            print("\r$(prefix)$(bar) $(round(pct, digits=1))% | Time: $(round(current_t, digits=0))s | V: $(round(v_val, digits=4))V   ")
            last_print = time()
        end
    end
    
    wall_time = time() - start_wall
    pct = clamp((current_t / total_t) * 100.0, 0.0, 100.0)
    filled = round(Int, 30 * (pct / 100.0))
    bar = "[" * repeat("=", filled) * repeat(" ", 30 - filled) * "]"
    prefix = is_truth ? "[Truth] " : "[Test]  "
    print("\r$(prefix)$(bar) 100.0% | Time: $(round(current_t, digits=0))s | Wall: $(round(wall_time, digits=2))s\n")
    
    return wall_time, hist_t, hist_v
end

function run_spme_benchmark()
    println("=======================================================")
    println("--- SPMe Micro-Optimization & Pareto Sweep ---")
    println("=======================================================")

    elec_params = Chen2020()
    elec_params.Vmin = 1.0; elec_params.Vmax = 5.0
    soc_init = 0.5
    elec_params.n.c₀ = (elec_params.n.z_0 + soc_init * (elec_params.n.z_100 - elec_params.n.z_0)) * elec_params.n.c₊
    elec_params.p.c₀ = (elec_params.p.z_0 + soc_init * (elec_params.p.z_100 - elec_params.p.z_0)) * elec_params.p.c₊
    
    pack_1C_amps = 5.0 # Single cell 1C
    
    stress_cycle = [
        RestStep(100.0),
        CurrentStep(-1.0 * pack_1C_amps, 400.0),  
        RestStep(100.0),
        CurrentStep(1.75 * pack_1C_amps, 150.0),  
        CurrentStep(-1.0 * pack_1C_amps, 400.0), 
        RestStep(300.0)
    ]
    exp = Experiment(stress_cycle)
    total_t = sum([s.period for s in stress_cycle])
    get_load_func = compile_experiment(exp)

    println("\n[!] Compiling Ground Truth (20 Nodes)...")
    sys_ref = build_benchmark_cell(elec_params, 20)
    
    time_ref, t_ref, V_ref = simulate_single_cell(sys_ref, total_t, get_load_func, 1e-8, true)

    nodes_list = [10]
    reltol_list = [1e-5]
    
    matrix_sweep = Iterators.product(nodes_list, reltol_list)
    total_runs = length(nodes_list) * length(reltol_list)
    results = []
    run_count = 1

    println("\n=======================================================")
    println(">>> INITIATING 2D PARETO SWEEP ($total_runs Configurations)")
    println("=======================================================")

    for (n_val, tol_val) in matrix_sweep
        println("\n[Run $run_count/$total_runs] Testing: Nodes=$(n_val) | reltol=$(tol_val)")
        
        sys_test = build_benchmark_cell(elec_params, n_val)
        wall_time, t_test, V_test = simulate_single_cell(sys_test, total_t, get_load_func, tol_val, false)
        
        min_len = min(length(V_ref), length(V_test))
        max_v_err = maximum(abs.(V_ref[1:min_len] .- V_test[1:min_len]))
        
        speedup = time_ref / wall_time
        
        push!(results, (
            nodes = n_val,
            reltol = tol_val,
            wall_time = wall_time,
            speedup = speedup,
            max_v_err = max_v_err,
            raw_time = copy(t_test[1:min_len]),
            raw_v = copy(V_test[1:min_len])
        ))
        
        sys_test = nothing
        GC.gc()
        
        run_count += 1
    end

    sweep_df = DataFrame(results)
    
    csv_out = select(sweep_df, Not([:raw_time, :raw_v]))
    CSV.write("benchmark_spme_results.csv", csv_out)
    println("\n[!] Sweep Complete. Metrics saved to 'benchmark_spme_results.csv'")

    valid_runs = filter(row -> row.max_v_err <= 0.010, sweep_df)
    best_run = nothing
    
    if !isempty(valid_runs)
        s_min, s_max = minimum(valid_runs.speedup), maximum(valid_runs.speedup)
        s_range = s_max > s_min ? (s_max - s_min) : 1.0
        
        e_min, e_max = minimum(valid_runs.max_v_err), maximum(valid_runs.max_v_err)
        e_range = e_max > e_min ? (e_max - e_min) : 1.0
        
        max_kneedle_dist = -Inf
        
        for row in eachrow(valid_runs)
            x_norm = (row.speedup - s_min) / s_range
            y_norm = (row.max_v_err - e_min) / e_range
            dist = x_norm - y_norm 
            
            if dist > max_kneedle_dist
                max_kneedle_dist = dist
                best_run = row
            end
        end
        println("\nOPTIMAL SPMe CONFIGURATION FOUND (< 10 mV Peak Error):")
    else
        println("\n[!] No configurations met the 10 mV criteria. Defaulting to lowest error run.")
        sort!(sweep_df, :max_v_err)
        best_run = sweep_df[1, :]
    end

    println("   Radial Nodes: $(best_run.nodes) | reltol: $(best_run.reltol)")
    println("   Execution Time: $(round(best_run.wall_time, digits=4))s (Speedup: $(round(best_run.speedup, digits=2))x)")
    println("   Peak Voltage Error: $(round(best_run.max_v_err * 1000.0, digits=2)) mV")

    unique_tols = sort(unique(sweep_df.reltol))
    tol_shapes = Dict(unique_tols[i] => [:circle, :square, :utriangle, :diamond][i] for i in 1:length(unique_tols))
    pt_shapes = [tol_shapes[t] for t in sweep_df.reltol]

    p1 = Plots.scatter(sweep_df.speedup, sweep_df.max_v_err .* 1000.0, 
        zcolor=sweep_df.nodes, markershape=pt_shapes, markersize=8,
        cmap=:viridis, colorbar_title="Nodes",
        xlabel="Speedup Factor [x]", ylabel="Peak Error [mV]",
        title="Pareto: Execution vs Voltage Error", label="")
    Plots.scatter!(p1, [best_run.speedup], [best_run.max_v_err .* 1000.0], markershape=:star5, markersize=14, color=:red, label="Optimal")
    Plots.hline!(p1, [10.0], color=:red, linestyle=:dash, label="10 mV Limit")

    p2 = Plots.plot(t_ref, V_ref, 
        label="Truth (20 nodes)", color=:black, linewidth=2,
        xlabel="Time [s]", ylabel="Voltage [V]", 
        title="Voltage Equivalence Trace")
    Plots.plot!(p2, best_run.raw_time, best_run.raw_v, 
        label="Optimal", color=:red, linestyle=:dash, linewidth=2)

    p3 = Plots.plot(t_ref[1:length(best_run.raw_v)], abs.(best_run.raw_v .- V_ref[1:length(best_run.raw_v)]) .* 1000.0, 
        label="Abs Error", color=:purple, linewidth=1, fill=(0, 0.2, :purple),
        xlabel="Time [s]", ylabel="Abs Error [mV]",
        title="Transient Voltage Error")
    Plots.hline!(p3, [10.0], color=:red, linestyle=:dash, label="10 mV Limit")

    final_dashboard = Plots.plot(p1, p2, p3, layout=(3, 1), size=(800, 1000), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
    Plots.display(final_dashboard)
end

run_spme_benchmark()