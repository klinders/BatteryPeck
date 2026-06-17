# ==============================================================================
# equivalence_proof.jl
# Mathematical equivalence proof: monolithic DAE vs. explicit co-simulation
# Contains:
# 1. run_equivalence_proof: Execute comparative study between monolithic and explicit co-simulation solvers
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using BatteryToolkit
using Plots

"""
    run_equivalence_proof()

Execute comparative study between monolithic and explicit co-simulation solvers.

Builds native MTK monolithic system and separate explicit co-simulation loop. Solves both systems against identical load profiles and calculates voltage variance to demonstrate mathematical equivalence.

# Returns
- Nothing (displays equivalence dashboard)
"""
function run_equivalence_proof()
    println("=====================================================")
    println("--- THE EXPLICIT CO-SIMULATION EQUIVALENCE TEST ---")
    println("=====================================================")
    
    elec_params = Chen2020()
    elec_params.Vmin = 1.0; elec_params.Vmax = 5.5

    # Define dynamic pulse profile for comparative testing
    exp = Experiment([
        RestStep(10.0),
        CurrentStep(5.0, 800.0),  # 5A Discharge
        RestStep(200.0),
        CurrentStep(-5.0, 400.0), # 5A Charge
        RestStep(200.0),
    ])
    
    get_load_func = compile_experiment(exp)

    # Build native monolithic system and structural simplification
    println("\n[1/4] Building Native Monolithic System...")
    @named cell_mono = SPMe(params=elec_params, side_reactions=true)
    @named load_mono = Current(); @named ground_mono = Ground()
    @parameters t Pin=0.0 Iin=0.0
    @variables V(t) I(t)
    
    eqs_mono = [
        V ~ load_mono.p.v - load_mono.n.v,
        I ~ Iin + Pin/V,
        load_mono.I.u ~ I,
        cell_mono.T.u ~ 298.15,
        connect(load_mono.p, cell_mono.p),
        connect(cell_mono.n, load_mono.n, ground_mono.g)
    ]
    @named sys_mono = ODESystem(eqs_mono, t, [V, I], [Pin, Iin]; systems=[cell_mono, load_mono, ground_mono])
    sys_simp_mono = structural_simplify(sys_mono)
    
    prob_mono = ODEProblem(sys_simp_mono, [sys_simp_mono.Pin => 0.0, sys_simp_mono.Iin => 0.0], (0.0, 2000.0); sparse=true, jac=false)

    # Solve monolithic reference system
    println("[2/4] Solving Monolithic System (High Precision)...")
    @time sol_mono = simulate(sys_simp_mono, prob_mono, exp, QNDF(autodiff=true); saveat=0.1, dtmax=1.0, reltol=1e-6, abstol=1e-6)
    
    # Initialise explicit co-simulation system and integrators
    println("\n[3/4] Building Explicit Co-Simulation System...")
    @named cell_exp = SPMe(params=elec_params, side_reactions=true)
    @named load_exp = Current(); @named ground_exp = Ground()
    @parameters t I_app=0.0
    
    eqs_exp = [
        cell_exp.T.u ~ 298.15,
        load_exp.I.u ~ I_app,
        connect(load_exp.p, cell_exp.p),
        connect(cell_exp.n, load_exp.n, ground_exp.g)
    ]
    @named sys_exp = ODESystem(eqs_exp, t, [], [I_app]; systems=[cell_exp, load_exp, ground_exp])
    sys_simp_exp = structural_simplify(sys_exp)

    local I_app_sym = nothing
    for p in parameters(sys_simp_exp)
        if contains(string(p), "I_app") I_app_sym = p; break end
    end

    prob_exp = ODEProblem(sys_simp_exp, [I_app_sym => 0.0], (0.0, 2000.0); sparse=true, jac=false)
    int_exp = init(prob_exp, QNDF(autodiff=true); reltol=1e-4, abstol=1e-4)

    # Iterate through explicit time-stepping loop
    println("[4/4] Running Explicit Step Loop (dt_sync = 0.1s)...")
    dt_sync = 0.1
    last_I = 0.0

    history_t = Float64[]; history_v_exp = Float64[]; history_i = Float64[]

    @time while true
        current_t = int_exp.t
        
        target_current, is_done = get_load_func(current_t, 4.0, 0.5)
        if is_done break end

        # Adjust time step size upon current signal transients
        if abs(target_current - last_I) > 1e-3
            set_proposed_dt!(int_exp, 1e-4)
            last_I = target_current
        end

        int_exp.ps[I_app_sym] = target_current
        SciMLBase.step!(int_exp, dt_sync, true)
        
        # Halt simulation if voltage limits are exceeded
        if int_exp.sol.retcode == SciMLBase.ReturnCode.Terminated
            println("  [!] Explicit Co-Sim hit voltage limits at t = $(current_t)s")
            break
        end

        push!(history_t, current_t)
        push!(history_v_exp, int_exp[sys_simp_exp.cell_exp.v])
        push!(history_i, target_current)
    end
    
    # Calculate error relative to monolithic truth and plot comparison results
    println("\nSimulation Complete! Calculating Error...")
    
    max_valid_t = min(sol_mono.t[end], history_t[end])
    valid_indices = history_t .<= max_valid_t
    
    t_valid = history_t[valid_indices]
    v_exp_valid = history_v_exp[valid_indices]
    i_valid = history_i[valid_indices]
    
    v_mono_truth = [sol_mono(t, idxs=sys_simp_mono.cell_mono.v) for t in t_valid]
    v_error_mv = abs.(v_exp_valid .- v_mono_truth) .* 1000.0

    p1 = plot(t_valid, v_mono_truth, label="Monolithic (Truth)", lw=4, color=:black, xlabel="Time (s)", ylabel="Voltage (V)", legend=:bottomright)
    plot!(p1, t_valid, v_exp_valid, label="Explicit Co-Sim", lw=2, color=:cyan, linestyle=:dash)
    
    p2 = plot(t_valid, v_error_mv, label="Absolute Error (mV)", lw=2, color=:red, xlabel="Time (s)", ylabel="Error (mV)")
    
    p3 = plot(t_valid, i_valid, label="Load Profile (A)", fill=true, fillalpha=0.2, color=:orange, xlabel="Time (s)", ylabel="Current (A)")
    
    display(plot(p1, p2, p3, layout=(3,1), size=(800,900), title="Equivalence Proof: Monolithic vs Explicit"))
    
    max_error = maximum(v_error_mv)
    println("-----------------------------------------------------")
    println("MAXIMUM ABSOLUTE ERROR: $(round(max_error, digits=4)) mV")
    println("-----------------------------------------------------")
end

Base.invokelatest(run_equivalence_proof)