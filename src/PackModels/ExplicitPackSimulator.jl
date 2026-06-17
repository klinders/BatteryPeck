# ==============================================================================
# ExplicitPackSimulator.jl
# Automated Master Orchestrator for 3D Pack Explicit Co-Simulation
# Contains:
# 1. get_sym: Retrieve symbol from system by matching parameter name
# 2. find_state: Locate state variable within system by string name
# 3. get_observed: Find observed equation left-hand side variable by string name
# 4. build_isolated_cell_template: Construct template ODE system for isolated SPMe cell
# 5. build_pack_simulator: Generate explicit pack simulator orchestrator struct
# 6. solve_pack_currents!: Calculate branch currents and total pack voltage
# 7. simulate_pack!: Execute full explicit co-simulation for battery pack
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal
using OrdinaryDiffEq
using BatteryToolkit
using CSV, DataFrames, DataInterpolations
using Logging
using Random

export ExplicitPackSimulator, build_pack_simulator, simulate_pack!

# Set random seed to ensure repeatable capacity distributions across experiments
const CAPACITY_SEED = 1234 

"""
    get_sym(sys, param_name)

Retrieve symbol from system by matching parameter name.

Searches through system parameters and returns first matching symbol. Throws error if parameter is missing.

# Arguments
- `sys`: ModelingToolkit system to search
- `param_name`: String containing parameter name to match

# Returns
- Matched parameter symbol
"""
function get_sym(sys, param_name)
    # Scan system parameters for string match and return symbol
    for p in parameters(sys)
        if contains(string(p), param_name) return p end
    end
    error("Parameter '$param_name' not found in system")
end

"""
    find_state(sys, name::String)

Locate state variable within system by string name.

Attempts to retrieve unknowns or falls back to states. Iterates through retrieved array and returns first match.

# Arguments
- `sys`: ModelingToolkit system to search
- `name::String`: Substring to locate within state names

# Returns
- Matched state symbol, or `nothing` if absent
"""
function find_state(sys, name::String)
    # Search state arrays for target substring and return match
    sts = try unknowns(sys) catch; states(sys) end
    for s in sts
        if occursin(name, string(s))
            return s
        end
    end
    return nothing
end

"""
    get_observed(sys, name::String)

Find observed equation left-hand side variable by string name.

Searches observed equations array and matches left-hand side string representation against target name.

# Arguments
- `sys`: ModelingToolkit system to evaluate
- `name::String`: Substring identifying target observed variable

# Returns
- Matched variable symbol, or `nothing` if absent
"""
function get_observed(sys, name::String)
    # Search observed equations for target left hand side variable
    obs_eqs = try observed(sys) catch; [] end
    for eq in obs_eqs
        if occursin(name, string(eq.lhs))
            return eq.lhs
        end
    end
    return nothing
end

"""
    build_isolated_cell_template(elec_params)

Construct template ODE system for isolated SPMe cell.

Initialises single cell with current load and ground components. Applies external temperature and applied current parameters to connections.

# Arguments
- `elec_params`: Electrical parameter set for SPMe model

# Output Variables
- `sys`: Structurally simplified ODE system representing cell
"""
function build_isolated_cell_template(elec_params)
    # Initialise SPMe cell with electrical load and ground components
    @named cell = SPMe(params=elec_params, side_reactions=true)
    @named load = Current(); @named ground = Ground()
    @parameters t I_app=0.0 T_ext=298.15
    
    # Formulate connections and constraint equations for ODE system
    eqs = [
        cell.T.u ~ T_ext,
        load.I.u ~ I_app,
        connect(load.p, cell.p),
        connect(cell.n, load.n, ground.g)
    ]
    
    sys = ODESystem(eqs, t, [], [I_app, T_ext]; systems=[cell, load, ground], name=:cell_template)
    return structural_simplify(sys)
end

struct ExplicitPackSimulator
    num_cells::Int
    rows_series::Int
    cols_parallel::Int
    cell_integrators::Vector
    therm_integrator
    sys_cell_template
    sys_therm
    I_sym
    T_sym
    Q_syms::Vector
    m_flow_sym 
    T_inlet_sym 
end

"""
    build_pack_simulator(geom::PackGeometry, elec_params, therm_params, rows_series::Int, cols_parallel::Int; verbose::Bool=true, geom_sigma::Float64=0.0)

Generate explicit pack simulator orchestrator struct.

Compiles single SPMe cell template and full 3D thermal LPTN system. Clones integrators for each cell in pack architecture to enable parallel state tracking.

# Arguments
- `geom::PackGeometry`: Pack spatial geometry definitions
- `elec_params`: Electrical parameter set
- `therm_params`: Thermal parameter set
- `rows_series::Int`: Number of cell rows connected in series
- `cols_parallel::Int`: Number of cell columns connected in parallel

# Keyword Arguments
- `verbose::Bool`: Flag to print compilation progress
- `geom_sigma::Float64`: Geometric standard deviation for capacity distribution

# Returns
- `ExplicitPackSimulator`: Fully initialised simulator structure
"""
function build_pack_simulator(geom::PackGeometry, elec_params, therm_params, rows_series::Int, cols_parallel::Int; verbose::Bool=true, geom_sigma::Float64=0.0)
    
    # Calculate pack dimensions and build single cell template
    num_cells = rows_series * cols_parallel
    if verbose println("[1/3] Compiling Single SPMe Template (Fast Mode)...") end
    sys_cell_template = build_isolated_cell_template(elec_params)
    @parameters t
    
    # Construct full thermal network and map volumetric heat inputs to respective cells
    if verbose println("[2/3] Compiling Full 3D Thermal & Fluid LPTN...") end
    sys_therm = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        sys_therm_base = build_pack_system(:therm_pack_base, geom, therm_params)
        
        @parameters Q_in[1:num_cells] = zeros(num_cells)
        cls = [getproperty(sys_therm_base, Symbol("cell_$i")) for i in 1:num_cells]
        forcing_eqs = [cls[i].Q_volumetric_in.u ~ Q_in[i] for i in 1:num_cells]
        
        sys_therm_wrapped = ODESystem(forcing_eqs, t, [], collect(Q_in); systems=[sys_therm_base], name=:therm_pack)
        structural_simplify(sys_therm_wrapped)
    end

    # Extract system symbols for external inputs and heat generation
    if verbose println("[3/3] Cloning 28 Integrators...") end
    I_sym = get_sym(sys_cell_template, "I_app")
    T_sym = get_sym(sys_cell_template, "T_ext")
    Q_syms = [get_sym(sys_therm, "Q_in[$i]") for i in 1:num_cells]
    
    m_flow_sym = get_sym(sys_therm, "m_flow_in")
    T_inlet_sym = get_sym(sys_therm, "T_inlet")

    # Duplicate initialised SPMe integrators for each cell in pack
    prob_template = ODEProblem(sys_cell_template, [I_sym => 0.0, T_sym => therm_params.ambient_temperature], (0.0, 31536000.0); sparse=true, jac=false)
    cell_integrators = []
    for i in 1:num_cells
        prob_clone = remake(prob_template, u0=copy(prob_template.u0), p=copy(prob_template.p))
        push!(cell_integrators, init(prob_clone, QNDF(autodiff=true); reltol=1e-5, abstol=1e-6, save_everystep=false))
    end
    
    # Initialise thermal integrator with starting guesses
    therm_p_guess = [Q_syms[i] => 0.0 for i in 1:num_cells]
    prob_therm = ODEProblem(sys_therm, therm_p_guess, (0.0, 31536000.0); sparse=true, jac=true)
    int_therm = init(prob_therm, QNDF(autodiff=false); reltol=1e-2, abstol=1e-3, save_everystep=false)

    return ExplicitPackSimulator(num_cells, rows_series, cols_parallel, cell_integrators, int_therm, sys_cell_template, sys_therm, I_sym, T_sym, Q_syms, m_flow_sym, T_inlet_sym)
end

"""
    solve_pack_currents!(I_branches::Vector{Float64}, V_cores::Vector{Float64}, G_branches::Vector{Float64}, rows::Int, cols::Int, I_total::Float64)

Calculate branch currents and total pack voltage for parallel-series architecture.

Iterates through series tiers to compute parallel voltage drops. Resolves individual cell currents based on local conductances.

# Arguments
- `I_branches::Vector{Float64}`: In-place array for output cell currents
- `V_cores::Vector{Float64}`: Array of cell core voltages
- `G_branches::Vector{Float64}`: Array of cell branch conductances
- `rows::Int`: Number of series-connected rows
- `cols::Int`: Number of parallel-connected columns
- `I_total::Float64`: Total applied pack current

# Returns
- Total calculated pack voltage
"""
function solve_pack_currents!(I_branches::Vector{Float64}, V_cores::Vector{Float64}, G_branches::Vector{Float64}, rows::Int, cols::Int, I_total::Float64)
    pack_voltage = 0.0
    
    # Iterate through series rows to compute parallel tier voltage drops
    for r in 1:rows
        start_idx = (r - 1) * cols + 1
        end_idx = r * cols
        
        V_tier = @view V_cores[start_idx:end_idx]
        G_tier = @view G_branches[start_idx:end_idx]
        
        V_drop_tier = (sum(V_tier .* G_tier) + I_total) / sum(G_tier)
        pack_voltage += V_drop_tier
        
        # Resolve individual branch currents based on local conductances
        for local_idx in 1:cols
            global_idx = start_idx + local_idx - 1
            I_branches[global_idx] = (V_drop_tier - V_cores[global_idx]) * G_branches[global_idx]
        end
    end
    return pack_voltage
end

"""
    simulate_pack!(cosim::ExplicitPackSimulator, exp::Experiment, therm_params::PackParameters, tms_strategy::TMSStrategy; total_time::Float64=0.0, dt_max::Float64=25.0, alpha::Float64=1.0, save_csv::Bool=true, force_even_current::Bool=false, is_isothermal::Bool=false, geom_sigma::Float64=0.0, m_active::Float64=0.03, m_passive::Float64=0.01)

Execute full explicit co-simulation for battery pack.

Orchestrates parallel cell integration alongside thermal system progression. Applies smart BMS tracker for dynamic load adjustment and manages active/passive TMS cooling logic.

# Arguments
- `cosim::ExplicitPackSimulator`: Initialised pack simulator struct
- `exp::Experiment`: Load profile and termination conditions
- `therm_params::PackParameters`: Thermal parameter definitions
- `tms_strategy::TMSStrategy`: Thermal management operational logic

# Keyword Arguments
- `total_time::Float64`: Expected simulation duration for progress tracking
- `dt_max::Float64`: Maximum permitted time step
- `alpha::Float64`: Step size scaling factor
- `save_csv::Bool`: Flag to write output data to disk
- `force_even_current::Bool`: Bypass resistance solver and split current equally
- `is_isothermal::Bool`: Clamp cell temperatures to ambient
- `geom_sigma::Float64`: Standard deviation for capacity variance proxy
- `m_active::Float64`: Active cooling mass flow rate
- `m_passive::Float64`: Passive resting mass flow rate

# Returns
- `DataFrame`: Compiled simulation results history
"""
function simulate_pack!(cosim::ExplicitPackSimulator, exp::Experiment, therm_params::PackParameters, tms_strategy::TMSStrategy; 
                        total_time::Float64=0.0, dt_max::Float64=25.0, alpha::Float64=1.0, save_csv::Bool=true, 
                        force_even_current::Bool=false, is_isothermal::Bool=false, geom_sigma::Float64=0.0,
                        m_active::Float64=0.03, m_passive::Float64=0.01)
    
    start_wall_time = time()
    
    # Configure output directory if saving results
    run_dir = ""
    if save_csv
        timestamp = round(Int, time())
        run_dir = joinpath(pwd(), "results_run_$(timestamp)")
        mkpath(run_dir)
    end
    
    # Compile load profiles and advance future lookahead window if required
    get_load_func = compile_experiment(exp)
    get_load_func_future = compile_experiment(exp)
    window = get_lookahead(tms_strategy)
    if window > 0.0
        for t_dummy in 0.0 : 1.0 : window get_load_func_future(t_dummy, 4.0 * cosim.rows_series, 0.5) end
    end
    
    # Generate linear interpolation model for DCIR lookup table
    lut_path = joinpath(@__DIR__, "..", "..", "data", "Chen2020", "soc_dcir_lut.csv")
    dcir_df = CSV.read(lut_path, DataFrame)
    sort!(dcir_df, :SoC) 
    dcir_interp = LinearInterpolation(dcir_df.DCIR_Ohms, dcir_df.SoC)
    
    # Initialise core tracking arrays and operational variables
    A_c = therm_params.tms_geometry.channel_width * therm_params.tms_geometry.channel_height * therm_params.tms_geometry.number_of_channels
    rho = therm_params.fluid.density
    
    num_cells = cosim.num_cells
    V_cores = zeros(num_cells); Q_gens = zeros(num_cells); T_cores = zeros(num_cells)
    I_branches = zeros(num_cells); G_branches = zeros(num_cells)
    last_I_total = 0.0
    
    # Apply capacity variance distribution using geometric standard deviation
    cap_multipliers = ones(num_cells)
    if geom_sigma > 0.0
        rng = Random.Xoshiro(CAPACITY_SEED)
        for i in 1:num_cells
            cap_multipliers[i] = 1.0 / max(0.5, 1.0 + randn(rng) * geom_sigma)
        end
    end
    
    # Preallocate history buffers for data chunking
    chunk_limit = 5000
    chunk_idx = 1
    
    history_t = Float64[]; sizehint!(history_t, chunk_limit)
    history_pack_v = Float64[]; sizehint!(history_pack_v, chunk_limit)
    history_I = Float64[]; sizehint!(history_I, chunk_limit)
    history_T_max = Float64[]; sizehint!(history_T_max, chunk_limit)
    history_vel = Float64[]; sizehint!(history_vel, chunk_limit)
    history_dt = Float64[]; sizehint!(history_dt, chunk_limit) 
    history_SoH = zeros(num_cells, chunk_limit)
    
    # Configure initial simulation state and tracking flags
    buffer_idx = 0
    in_memory_master_df = DataFrame() 
    
    last_save_t = -999.0
    last_save_I = 0.0
    last_save_V = 0.0
    last_save_T = 0.0
    last_save_flow = 0.0
    
    current_dt_sync = 0.1
    
    initial_T_C = therm_params.ambient_temperature - 273.15
    prev_step_T_max = initial_T_C
    T_max_C = initial_T_C
    current_t = 0.0
    is_done = false
    
    R_sei_sym = find_state(cosim.sys_cell_template, "R_sei")
    L_sei_sym = find_state(cosim.sys_cell_template, "L_sei")
    SoH_sym = get_observed(cosim.sys_cell_template, "SoH")
    if isnothing(SoH_sym) SoH_sym = try getproperty(cosim.sys_cell_template, :SoH) catch; nothing end end
    
    R_sei_update_interval = 3600.0
    last_R_sei_update_t = -3600.0 
    current_R_sei = zeros(num_cells)

    function flush_chunk!(b_idx)
        # Construct and flush history buffer to disk or memory
        if b_idx == 0 return end
        df = DataFrame(
            Time_s = history_t[1:b_idx],
            Pack_Voltage_V = history_pack_v[1:b_idx],
            Pack_Current_A = history_I[1:b_idx],
            Max_Temp_C = history_T_max[1:b_idx],
            Velocity_ms = history_vel[1:b_idx],
            dt_sync_s = history_dt[1:b_idx]
        )
        for i in 1:num_cells df[!, Symbol("SoH_Cell_$i")] = history_SoH[i, 1:b_idx] end
        
        if save_csv
            CSV.write(joinpath(run_dir, "chunk_$(lpad(chunk_idx, 4, '0')).csv"), df)
        else
            append!(in_memory_master_df, df)
        end
        
        empty!(history_t); empty!(history_pack_v); empty!(history_I); empty!(history_T_max); empty!(history_vel); empty!(history_dt)
        chunk_idx += 1
        return 0 
    end

    last_print_wall_time = time()

    while true
        # Retrieve current time and estimate pack voltage
        current_t = cosim.cell_integrators[1].t
        approx_pack_v = length(history_pack_v) > 0 ? history_pack_v[end] : (4.0 * cosim.rows_series)
        
        # Track minimum and maximum cell voltage and state of charge for BMS
        min_cell_v = Inf; max_cell_v = -Inf
        min_soc = Inf; max_soc = -Inf
        
        for i in 1:num_cells
            v_c = cosim.cell_integrators[i][cosim.sys_cell_template.cell.v]
            soc_c = cosim.cell_integrators[i][cosim.sys_cell_template.cell.soc]
            
            if v_c < min_cell_v min_cell_v = v_c end
            if v_c > max_cell_v max_cell_v = v_c end
            if soc_c < min_soc min_soc = soc_c end
            if soc_c > max_soc max_soc = soc_c end
        end
        
        # Determine worst case operational bounds based on intended current direction
        intended_I, _ = get_load_func(current_t, (4.0 * cosim.rows_series), 0.5)
        
        if intended_I > 0
            worst_v = min_cell_v * cosim.rows_series
            worst_soc = min_soc
        elseif intended_I < 0
            worst_v = max_cell_v * cosim.rows_series
            worst_soc = max_soc
        else
            worst_v = approx_pack_v
            worst_soc = max_soc
        end
        
        # Fetch actual total current and check termination condition
        I_total, is_done = get_load_func(current_t, worst_v, worst_soc)
        if is_done break end
        
        # Periodically update SEI resistance values across all cells
        if current_t - last_R_sei_update_t >= R_sei_update_interval
            for i in 1:num_cells
                if !isnothing(R_sei_sym) current_R_sei[i] = cosim.cell_integrators[i][R_sei_sym]
                elseif !isnothing(L_sei_sym) current_R_sei[i] = cosim.cell_integrators[i][L_sei_sym] / 5e-6
                else current_R_sei[i] = 0.002 end
            end
            last_R_sei_update_t = current_t
        end
        
        # Calculate local virtual resistance and core temperatures
        for i in 1:num_cells
            V_term = cosim.cell_integrators[i][cosim.sys_cell_template.cell.v]
            I_prev = cosim.cell_integrators[i].ps[cosim.I_sym]
            cell_soc = clamp(cosim.cell_integrators[i][cosim.sys_cell_template.cell.soc], 0.0, 1.0)
            
            R_virt = dcir_interp(cell_soc) + current_R_sei[i]
            
            G_branches[i] = 1.0 / R_virt
            V_cores[i] = V_term - (I_prev * R_virt)
            Q_gens[i] = cosim.cell_integrators[i][cosim.sys_cell_template.cell.Q_total]
            
            core_T_var = getproperty(getproperty(cosim.sys_therm.therm_pack_base, Symbol("cell_$i")), :core_cap).T
            
            if is_isothermal
                T_cores[i] = therm_params.ambient_temperature 
            else
                T_cores[i] = cosim.therm_integrator[core_T_var]
            end
        end

        # Resolve current distribution based on solver mode
        if force_even_current
            pack_voltage = 0.0
            I_cell = I_total / cosim.cols_parallel
            for r in 1:cosim.rows_series
                tier_v = 0.0
                for c in 1:cosim.cols_parallel
                    idx = (r - 1) * cosim.cols_parallel + c
                    predicted_v = V_cores[idx] + (I_cell / G_branches[idx])
                    tier_v += predicted_v
                    I_branches[idx] = I_cell
                end
                pack_voltage += tier_v / cosim.cols_parallel
            end
        else
            pack_voltage = solve_pack_currents!(I_branches, V_cores, G_branches, cosim.rows_series, cosim.cols_parallel, I_total)
        end
        
        # Evaluate future load for predictive thermal management
        T_max_C = maximum(T_cores) - 273.15
        
        if window > 0.0
            future_I, _ = get_load_func_future(current_t + window, worst_v, worst_soc)
            cell_future_load = abs(future_I) / cosim.cols_parallel
        else
            cell_future_load = 0.0
        end
        
        # Override TMS flow and temperature based on active cooling strategy
        if m_active == 0.0 && m_passive == 0.0
            target_flow = 0.0
            target_T = therm_params.inlet_temperature
        else
            current_flow = cosim.therm_integrator.ps[cosim.m_flow_sym]
            target_flow, dynamic_T = evaluate_tms_state(tms_strategy, T_max_C, cell_future_load, current_flow, m_passive, m_active)
            
            if target_flow >= m_active && m_active > 0.0
                target_T = dynamic_T
            else
                target_T = therm_params.inlet_temperature
            end
        end
        
        # Apply updated parameters to individual cell and thermal integrators
        cosim.therm_integrator.ps[cosim.m_flow_sym] = target_flow
        cosim.therm_integrator.ps[cosim.T_inlet_sym] = target_T
        
        for i in 1:num_cells
            cosim.cell_integrators[i].ps[cosim.I_sym] = I_branches[i] * cap_multipliers[i]
            cosim.cell_integrators[i].ps[cosim.T_sym] = T_cores[i]
            cosim.therm_integrator.ps[cosim.Q_syms[i]] = Q_gens[i]
        end
        
        # Trigger SciML modification flag if load current shifts significantly
        if abs(I_total - last_I_total) > 1e-3
            for i in 1:num_cells
                SciMLBase.u_modified!(cosim.cell_integrators[i], true)
            end
        end
        
        # Evaluate transient state to determine data logging eligibility
        time_since_save = current_t - last_save_t
        is_transient   = abs(I_total - last_save_I) > 1e-3 || abs(target_flow - last_save_flow) > 1e-6
        is_fast_moving = abs(pack_voltage - last_save_V) > 0.01 || abs(T_max_C - last_save_T) > 0.05
        is_active      = abs(I_total) > 1e-3 
        
        if is_transient || (is_fast_moving && time_since_save >= 0.2) || (is_active && time_since_save >= 1.0) || (!is_active && time_since_save >= 10.0)
            # Append current state to history buffers and print progress
            buffer_idx += 1
            push!(history_t, current_t)
            push!(history_pack_v, pack_voltage)
            push!(history_T_max, T_max_C)
            push!(history_vel, target_flow / (rho * A_c))
            push!(history_I, I_total)
            push!(history_dt, current_dt_sync)
            
            if !isnothing(SoH_sym)
                for i in 1:num_cells history_SoH[i, buffer_idx] = cosim.cell_integrators[i][SoH_sym] end
            else
                for i in 1:num_cells history_SoH[i, buffer_idx] = 100.0 end
            end
            
            last_save_t = current_t
            last_save_I = I_total
            last_save_V = pack_voltage
            last_save_T = T_max_C
            last_save_flow = target_flow
            
            current_wall_time = time()
            if current_wall_time - last_print_wall_time >= 0.5
                if total_time > 0.0
                    pct = clamp((current_t / total_time) * 100.0, 0.0, 100.0)
                    filled = round(Int, 30 * (pct / 100.0))
                    bar = "[" * repeat("=", filled) * repeat(" ", 30 - filled) * "]"
                    print("\r$(bar) $(round(pct, digits=1))% | Time: $(round(current_t, digits=0))s | dt: $(round(current_dt_sync, digits=2))s | Max T: $(round(T_max_C, digits=2))°C        ")
                else
                    print("\r[Running] Time: $(round(current_t, digits=0))s | dt: $(round(current_dt_sync, digits=2))s | Max T: $(round(T_max_C, digits=2))°C        ")
                end
                last_print_wall_time = current_wall_time
            end
            
            if buffer_idx >= chunk_limit
                buffer_idx = flush_chunk!(buffer_idx)
            end
        end

        # Propose adaptive synchronisation time step based on temperature derivative
        dT_dt = abs(T_max_C - prev_step_T_max) / current_dt_sync
        epsilon = 1e-4 
        
        if is_transient
            current_dt_sync = 0.1
        else
            proposed_dt = alpha / (dT_dt + epsilon)
            current_dt_sync = clamp(proposed_dt, 0.1, dt_max)
        end
        
        prev_step_T_max = T_max_C
        last_I_total = I_total

        # Step individual cell ODE solvers forward in parallel
        t_target = current_t + current_dt_sync
        stalled = Threads.Atomic{Bool}(false) 
        
        Threads.@threads for i in 1:num_cells 
            int_c = cosim.cell_integrators[i]
            SciMLBase.add_tstop!(int_c, t_target)
            substeps = 0
            
            while int_c.t < t_target && !stalled[]
                SciMLBase.step!(int_c)
                substeps += 1
                
                if substeps > 1000
                    Threads.atomic_xchg!(stalled, true)
                    break
                end
                
                if stalled[] || any(int_c -> int_c.sol.retcode == SciMLBase.ReturnCode.Terminated, cosim.cell_integrators)
                    if stalled[]
                        println("\n[!] CRITICAL: ODE Solver Stalled! (Likely Surface Lithium Depletion)")
                    end
                    break
                end
            end
        end
        
        # Step thermal system forward to match cell time target
        if !is_isothermal
            SciMLBase.add_tstop!(cosim.therm_integrator, t_target)
            while cosim.therm_integrator.t < t_target && !stalled[]
                SciMLBase.step!(cosim.therm_integrator)
            end
        end
        
        # Halt simulation if solvers stall or terminate
        if stalled[] || any(int_c -> int_c.sol.retcode == SciMLBase.ReturnCode.Terminated, cosim.cell_integrators)
            break
        end
    end 
    
    # Calculate final elapsed time and print completion metrics
    elapsed_total = round(time() - start_wall_time, digits=2)
    if total_time > 0.0
        pct = is_done ? 100.0 : clamp((current_t / total_time) * 100.0, 0.0, 100.0)
        filled = round(Int, 30 * (pct / 100.0))
        bar = "[" * repeat("=", filled) * repeat(" ", 30 - filled) * "]"
        print("\r$(bar) $(round(pct, digits=1))% | Time: $(round(current_t, digits=0))s | Max T: $(round(T_max_C, digits=2))°C | Wall: $(elapsed_total)s\n")
    else
        print("\r[Done] Time: $(round(current_t, digits=0))s | Max T: $(round(T_max_C, digits=2))°C | Wall: $(elapsed_total)s\n")
    end
    
    # Consolidate intermediate chunk files into master results dataset
    flush_chunk!(buffer_idx)
    
    master_df = DataFrame()
    if save_csv
        chunk_files = sort(filter(f -> endswith(f, ".csv") && startswith(f, "chunk_"), readdir(run_dir)))
        for file in chunk_files
            path = joinpath(run_dir, file)
            append!(master_df, CSV.read(path, DataFrame))
            rm(path) 
        end
        master_path = joinpath(run_dir, "master_results.csv")
        CSV.write(master_path, master_df)
    else
        master_df = in_memory_master_df
    end
    
    return master_df
end