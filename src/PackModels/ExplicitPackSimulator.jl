# ==============================================================================
# ExplicitPackSimulator.jl
# Automated Master Orchestrator for 3D Pack Explicit Co-Simulation
# Contains:
# 1. get_sym: Retrieve symbol from system by matching parameter name
# 2. find_state: Locate state variable within system unknowns or states
# 3. get_observed: Retrieve observed equation left-hand side matching string
# 4. build_isolated_cell_template: Generate isolated cell template system
# 5. build_pack_simulator: Construct explicit pack simulator with coupled physics
# 6. reset_simulator!: Reinitialise cell and thermal integrators
# 7. solve_pack_currents!: Compute parallel branch currents and pack voltage
# 8. simulate_pack!: Execute transient multi-physics simulation of battery pack
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal
using OrdinaryDiffEq
using BatteryToolkit
using CSV, DataFrames, DataInterpolations
using Logging
using Random

export ExplicitPackSimulator, build_pack_simulator, reset_simulator!, simulate_pack!

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
    # Loop through system parameters and return parameter if string matches name
    for p in parameters(sys)
        if contains(string(p), param_name) return p end
    end
    # Throw error if parameter not found
    error("Parameter '$param_name' not found in system")
end

"""
    find_state(sys, name)

Locate state variable within system unknowns or states.

Iterates through available system states and returns first element matching target name string.

# Arguments
- `sys`: ModelingToolkit system to search
- `name`: String containing state name to match

# Returns
- Matched state variable or nothing if absent
"""
function find_state(sys, name::String)
    # Attempt to retrieve system unknowns or states and iterate to find match
    sts = try unknowns(sys) catch; states(sys) end
    for s in sts
        if occursin(name, string(s)) return s end
    end
    return nothing
end

"""
    get_observed(sys, name)

Retrieve observed equation left-hand side matching string.

Searches system observed equations and returns left-hand side if variable name matches.

# Arguments
- `sys`: ModelingToolkit system to search
- `name`: String containing observed variable name

# Returns
- Left-hand side of observed equation or nothing if absent
"""
function get_observed(sys, name::String)
    # Attempt to retrieve observed equations and iterate to find matching variable
    obs_eqs = try observed(sys) catch; [] end
    for eq in obs_eqs
        if occursin(name, string(eq.lhs)) return eq.lhs end
    end
    return nothing
end

"""
    build_isolated_cell_template(elec_params)

Generate isolated cell template system.

Constructs base electrochemical system with applied current and external temperature boundaries.

# Arguments
- `elec_params`: Struct containing cell electrochemical parameters

# Returns
- Structurally simplified ODESystem representing isolated cell
"""
function build_isolated_cell_template(elec_params)
    # Initialise base components and define external boundary parameters
    @named cell = SPMe(params=elec_params, side_reactions=true)
    @named load = Current()
    @named ground = Ground()
    @parameters t I_app=0.0 dI_dt_app=0.0 t_start=0.0 T_ext=298.15
    
    # Formulate governing equations and return simplified system
    eqs = [
        cell.T.u ~ T_ext,
        load.I.u ~ I_app + dI_dt_app * (t - t_start), 
        connect(load.p, cell.p),
        connect(cell.n, load.n, ground.g)
    ]
    sys = ODESystem(eqs, t, [], [I_app, dI_dt_app, t_start, T_ext]; systems=[cell, load, ground], name=:cell_template)
    return structural_simplify(sys)
end

mutable struct ExplicitPackSimulator
    num_cells::Int
    rows_series::Int
    cols_parallel::Int
    cell_integrators::Vector
    therm_integrator
    sys_cell_template
    sys_therm
    I_sym
    dI_dt_sym    
    t_start_sym  
    T_sym
    Q_syms::Vector
    dQ_dt_syms::Vector   
    t_start_therm_sym    
    m_flow_sym 
    T_inlet_sym 
    prob_template
    prob_therm
end

"""
    build_pack_simulator(geom, elec_params, therm_params, rows_series, cols_parallel; verbose, geom_sigma, R_axial_val, R_radial_val, R_contact_val)

Construct explicit pack simulator with coupled physics.

Compiles separate electrochemical cell templates and holistic 3D fluid thermal network.

# Arguments
- `geom`: PackGeometry layout struct
- `elec_params`: Cell electrochemical parameters
- `therm_params`: Pack thermal parameters
- `rows_series`: Integer cell count in series
- `cols_parallel`: Integer cell count in parallel
- `verbose`: Boolean flag to print initialisation progress
- `geom_sigma`: Float64 standard deviation for capacity distribution
- `R_axial_val`: Float64 axial thermal resistance
- `R_radial_val`: Float64 radial thermal resistance
- `R_contact_val`: Float64 contact thermal resistance

# Returns
- Initialised ExplicitPackSimulator instance
"""
function build_pack_simulator(geom::PackGeometry, elec_params, therm_params, rows_series::Int, cols_parallel::Int; verbose::Bool=true, geom_sigma::Float64=0.0, R_axial_val::Float64=0.25, R_radial_val::Float64=0.2, R_contact_val::Float64=4.0)
    # Determine total cell count and build base isolated cell template
    num_cells = rows_series * cols_parallel
    if verbose println("[1/3] Compiling Single SPMe Template (Fast Mode)...") end
    sys_cell_template = build_isolated_cell_template(elec_params)
    @parameters t
    
    # Construct and simplify full thermal network with external heat generation forcing equations
    if verbose println("[2/3] Compiling Full 3D Thermal & Fluid LPTN (FOH Enabled)...") end
    sys_therm = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        sys_therm_base = build_pack_system(:therm_pack_base, geom, therm_params; R_axial_val=R_axial_val, R_radial_val=R_radial_val, R_contact_val=R_contact_val)
        @parameters Q_in[1:num_cells] = zeros(num_cells)
        @parameters dQ_dt_in[1:num_cells] = zeros(num_cells)
        @parameters t_start_therm = 0.0
        
        cls = [getproperty(sys_therm_base, Symbol("cell_$i")) for i in 1:num_cells]
        forcing_eqs = [cls[i].Q_volumetric_in.u ~ Q_in[i] + dQ_dt_in[i] * (t - t_start_therm) for i in 1:num_cells]
        
        sys_therm_wrapped = ODESystem(forcing_eqs, t, [], vcat(collect(Q_in), collect(dQ_dt_in), [t_start_therm]); systems=[sys_therm_base], name=:therm_pack)
        structural_simplify(sys_therm_wrapped)
    end

    # Retrieve boundary condition symbols and clone distinct electrochemical integrators for each cell
    if verbose println("[3/3] Cloning $(num_cells) Integrators...") end
    I_sym = get_sym(sys_cell_template, "I_app")
    dI_dt_sym = get_sym(sys_cell_template, "dI_dt_app")
    t_start_sym = get_sym(sys_cell_template, "t_start")
    T_sym = get_sym(sys_cell_template, "T_ext")
    
    Q_syms = [get_sym(sys_therm, "Q_in[$i]") for i in 1:num_cells]
    dQ_dt_syms = [get_sym(sys_therm, "dQ_dt_in[$i]") for i in 1:num_cells]
    t_start_therm_sym = get_sym(sys_therm, "t_start_therm")
    m_flow_sym = get_sym(sys_therm, "m_flow_in")
    T_inlet_sym = get_sym(sys_therm, "T_inlet")

    prob_template = ODEProblem(sys_cell_template, [I_sym => 0.0, dI_dt_sym => 0.0, t_start_sym => 0.0, T_sym => therm_params.ambient_temperature], (0.0, 315360000.0); sparse=true, jac=false, warn_initialize_determined = false)
    cell_integrators = []
    for i in 1:num_cells
        prob_clone = remake(prob_template, u0=copy(prob_template.u0), p=copy(prob_template.p))
        push!(cell_integrators, init(prob_clone, QNDF(autodiff=true); reltol=1e-5, abstol=1e-6, save_everystep=false, maxiters=Int(1e7), verbose=false))
    end
    
    # Initialise thermal integrator and return final explicit simulator instance
    therm_p_guess = vcat([Q_syms[i] => 0.0 for i in 1:num_cells], [dQ_dt_syms[i] => 0.0 for i in 1:num_cells], [t_start_therm_sym => 0.0])
    prob_therm = ODEProblem(sys_therm, therm_p_guess, (0.0, 315360000.0); sparse=true, jac=true, warn_initialize_determined = false)
    int_therm = init(prob_therm, QNDF(autodiff=false); reltol=1e-2, abstol=1e-3, save_everystep=false, maxiters=Int(1e7), verbose=false)

    return ExplicitPackSimulator(num_cells, rows_series, cols_parallel, cell_integrators, int_therm, sys_cell_template, sys_therm, I_sym, dI_dt_sym, t_start_sym, T_sym, Q_syms, dQ_dt_syms, t_start_therm_sym, m_flow_sym, T_inlet_sym, prob_template, prob_therm)
end

"""
    reset_simulator!(cosim)

Reinitialise cell and thermal integrators.

Resets existing OrdinaryDiffEq integrators to original problem states.

# Arguments
- `cosim`: Target ExplicitPackSimulator instance
"""
function reset_simulator!(cosim::ExplicitPackSimulator)
    # Loop through and reset each electrochemical cell integrator alongside thermal integrator
    for i in 1:cosim.num_cells
        cosim.cell_integrators[i] = init(remake(cosim.prob_template), QNDF(autodiff=true); reltol=1e-5, abstol=1e-6, save_everystep=false, maxiters=Int(1e7), verbose=false)
    end
    cosim.therm_integrator = init(remake(cosim.prob_therm), QNDF(autodiff=false); reltol=1e-2, abstol=1e-3, save_everystep=false, maxiters=Int(1e7), verbose=false)
end

"""
    solve_pack_currents!(I_branches, V_cores, G_branches, rows, cols, I_total)

Compute parallel branch currents and pack voltage.

Divides current based solely on resistance to guarantee numerical stability.

# Arguments
- `I_branches`: Vector to store calculated branch currents
- `V_cores`: Vector of cell core voltages
- `G_branches`: Vector of cell conductances
- `rows`: Integer count of series rows
- `cols`: Integer count of parallel columns
- `I_total`: Float64 total pack current constraint

# Returns
- Float64 computed total pack voltage
"""
function solve_pack_currents!(I_branches::Vector{Float64}, V_cores::Vector{Float64}, G_branches::Vector{Float64}, rows::Int, cols::Int, I_total::Float64)
    # Initialise pack voltage and iterate through parallel tiers
    pack_voltage = 0.0
    for r in 1:rows
        start_idx = (r - 1) * cols + 1
        end_idx = r * cols
        
        V_tier = @view V_cores[start_idx:end_idx]
        G_tier = @view G_branches[start_idx:end_idx]
        
        # Compute tier voltage drop and append to total pack voltage
        sum_G_tier = sum(G_tier)
        V_drop_tier = (sum(V_tier .* G_tier) + I_total) / sum_G_tier
        pack_voltage += V_drop_tier
        
        # Distribute current stably omitting open circuit voltage balancing
        for local_idx in 1:cols
            global_idx = start_idx + local_idx - 1
            I_branches[global_idx] = I_total * (G_branches[global_idx] / sum_G_tier)
        end
    end
    return pack_voltage
end

"""
    simulate_pack!(cosim, exp, therm_params, tms_strategy; ...)

Execute transient multi-physics simulation of battery pack.

Steps through experimental protocol managing thermal, electrical, and fluid states.

# Arguments
- `cosim`: ExplicitPackSimulator instance
- `exp`: Experiment struct defining schedule
- `therm_params`: PackParameters detailing thermal properties
- `tms_strategy`: TMSStrategy detailing cooling logic

# Returns
- DataFrame containing logged simulation history
"""
function simulate_pack!(cosim::ExplicitPackSimulator, exp::Experiment, therm_params::PackParameters, tms_strategy::TMSStrategy;
                        total_time::Float64=0.0, save_csv::Bool=true, 
                        derate_exponent::Float64=1.0, print_prefix::String="",
                        force_even_current::Bool=false, is_isothermal::Bool=false, geom_sigma::Float64=0.0,
                        m_active::Float64=0.03, m_passive::Float64=1e-6, heat_multiplier::Float64=1.0,
                        
                        opt_w_ratio::Dict{Symbol, Float64}=Dict{Symbol, Float64}(:smooth => 1.75, :drive => 4.25, :grid => 2.5, :aging => 2.0, :rest => 1.0), 
                        opt_alpha_elec::Dict{Symbol, Float64}=Dict{Symbol, Float64}(:smooth => 30.0, :drive => 70.0, :grid => 7.75, :aging => 20.0, :rest => 10.0), 
                        opt_dt_max_elec::Dict{Symbol, Float64}=Dict{Symbol, Float64}(:smooth => 1000.0, :drive => 1.25, :grid => 12.5, :aging => 400.0, :rest => 3600.0),
                        
                        opt_alpha_therm::Dict{Symbol, Float64}=Dict{Symbol, Float64}(:smooth => 0.1, :drive => 12.0, :grid => 12.0, :aging => 0.5, :rest => 0.1),
                        opt_dt_max_therm::Dict{Symbol, Float64}=Dict{Symbol, Float64}(:smooth => 558.3, :drive => 150.0, :grid => 150.0, :aging => 2062.5, :rest => 3600.0),
                        
                        verbose::Bool=true, dense_logging::Bool=false, sparse_logging::Bool=false, ultra_sparse_logging::Bool=false)
    
    # Initialise runtime variables and structural constraints
    start_wall_time = time()
    DT_MIN_ELEC = 0.01
    V_MAX = 5.0 * cosim.rows_series
    I_MAX = 7.5 * cosim.cols_parallel
    EPSILON_ELEC = 0.01
    
    start_t = cosim.therm_integrator.t
    abs_total_time = start_t + total_time
    
    v_nom_min = 2.5 
    v_nom_max = 4.2 
    try v_nom_min = cosim.sys_cell_template.cell.p.Vmin catch end
    try v_nom_max = cosim.sys_cell_template.cell.p.Vmax catch end
    V_MIN_BMS = v_nom_min * cosim.rows_series
    V_MAX_BMS = v_nom_max * cosim.rows_series
    
    I_1C_PACK = 5.0 * cosim.cols_parallel
    MAX_CHG_BASE = 0.7 * I_1C_PACK  
    MAX_DSG_BASE = -1.5 * I_1C_PACK 
    
    dyn_chg_limit = MAX_CHG_BASE
    dyn_dsg_limit = MAX_DSG_BASE

    # Create output directory if saving enabled
    run_dir = ""
    if save_csv
        timestamp = round(Int, time())
        run_dir = joinpath(pwd(), "results_run_$(timestamp)")
        mkpath(run_dir)
    end
    
    # Compile dynamic schedule closures and load internal resistance maps
    get_load_func = compile_experiment(exp)
    get_load_func_future = compile_experiment(exp)
    boundaries = vcat([0.0], exp.tstops) .+ start_t
    
    lut_path = joinpath(@__DIR__, "..", "..", "data", "Chen2020", "soc_dcir_lut.csv")
    dcir_df = CSV.read(lut_path, DataFrame)
    sort!(dcir_df, :SoC) 
    dcir_interp = LinearInterpolation(dcir_df.DCIR_Ohms, dcir_df.SoC)
    
    # Extract physical dimensions and instantiate state tracking vectors
    A_c = therm_params.tms_geometry.channel_width * therm_params.tms_geometry.channel_height * therm_params.tms_geometry.number_of_channels
    rho = therm_params.fluid.density
    num_cells = cosim.num_cells
    V_cores = zeros(num_cells)
    Q_gens = zeros(num_cells)
    last_Q_gens = zeros(num_cells)
    T_cores = zeros(num_cells)
    
    for i in 1:num_cells
        core_T_var = getproperty(getproperty(cosim.sys_therm.therm_pack_base, Symbol("cell_$i")), :core_cap).T
        T_cores[i] = cosim.therm_integrator[core_T_var]
    end
    prev_T_cores = copy(T_cores)
    
    I_branches = zeros(num_cells)
    G_branches = zeros(num_cells)
    
    # Distribute capacity multipliers using geometric standard deviation
    cap_multipliers = ones(num_cells)
    if geom_sigma > 0.0
        rng = Random.Xoshiro(CAPACITY_SEED)
        for i in 1:num_cells cap_multipliers[i] = 1.0 / max(0.5, 1.0 + randn(rng) * geom_sigma) end
    end
    
    # Preallocate history buffers for high speed data accumulation
    chunk_limit = 5000
    chunk_idx = 1
    history_t = Float64[]; sizehint!(history_t, chunk_limit)
    history_pack_v = Float64[]; sizehint!(history_pack_v, chunk_limit)
    history_pack_soc = Float64[]; sizehint!(history_pack_soc, chunk_limit)
    history_I = Float64[]; sizehint!(history_I, chunk_limit)
    history_T_max = Float64[]; sizehint!(history_T_max, chunk_limit)
    history_vel = Float64[]; sizehint!(history_vel, chunk_limit)
    history_dt = Float64[]; sizehint!(history_dt, chunk_limit) 
    history_limit_chg = Float64[]; sizehint!(history_limit_chg, chunk_limit)
    history_limit_dsg = Float64[]; sizehint!(history_limit_dsg, chunk_limit)
    
    history_SoH = zeros(num_cells, chunk_limit)
    history_I_branches = zeros(num_cells, chunk_limit)
    history_T_cores = zeros(num_cells, chunk_limit)
    history_V_cells = zeros(num_cells, chunk_limit)

    buffer_idx = 0
    in_memory_master_df = DataFrame() 
    last_save_t = -999.0; last_save_I = 0.0; last_save_V = 0.0; last_save_T = 0.0; last_save_flow = 0.0
    
    initial_T_C = maximum(T_cores) - 273.15
    prev_step_T_max = initial_T_C
    T_max_C = initial_T_C
    
    current_t_therm = start_t
    current_t_elec = start_t
    prev_dt_therm = 0.1; current_dt_therm = 0.1
    is_done = false
    last_pack_v_elec = 4.0 * cosim.rows_series; last_I_total_elec = 0.0; dt_elec = 0.1
    
    # Find critical state symbols for degradation monitoring
    R_sei_sym = find_state(cosim.sys_cell_template, "R_sei")
    L_sei_sym = find_state(cosim.sys_cell_template, "L_sei")
    SoH_sym = get_observed(cosim.sys_cell_template, "SoH")
    if isnothing(SoH_sym) SoH_sym = try getproperty(cosim.sys_cell_template, :SoH) catch; nothing end end
    
    R_sei_update_interval = 3600.0; last_R_sei_update_t = -3600.0
    current_R_sei = zeros(num_cells)
    last_print_wall_time = time()

    # Define inline closure to offload chunks to disk and clear buffers
    function flush_chunk!(b_idx)
        if b_idx == 0 return end
        df = DataFrame(
            Time_s = history_t[1:b_idx], 
            Pack_Voltage_V = history_pack_v[1:b_idx], 
            Pack_Current_A = history_I[1:b_idx],
            BMS_Limit_Chg_A = history_limit_chg[1:b_idx],
            BMS_Limit_Dsg_A = history_limit_dsg[1:b_idx],
            Pack_SoC = history_pack_soc[1:b_idx],
            Max_Temp_C = history_T_max[1:b_idx], 
            Velocity_ms = history_vel[1:b_idx], 
            dt_sync_s = history_dt[1:b_idx]
        )
        for i in 1:num_cells 
            df[!, Symbol("V_Cell_$i")] = history_V_cells[i, 1:b_idx]
            df[!, Symbol("SoH_Cell_$i")] = history_SoH[i, 1:b_idx] 
            df[!, Symbol("I_Branch_Cell_$i")] = history_I_branches[i, 1:b_idx]
            df[!, Symbol("T_Core_Cell_$i")] = history_T_cores[i, 1:b_idx]
        end
        if save_csv CSV.write(joinpath(run_dir, "chunk_$(lpad(chunk_idx, 4, '0')).csv"), df)
        else append!(in_memory_master_df, df) end
        
        empty!(history_t); empty!(history_pack_v); empty!(history_pack_soc); empty!(history_I); 
        empty!(history_T_max); empty!(history_vel); empty!(history_dt); 
        empty!(history_limit_chg); empty!(history_limit_dsg)
        
        chunk_idx += 1
        return 0 
    end

    # Execute primary simulation loop alternating between thermal and electrical advances
    while current_t_therm < abs_total_time && !is_done
        # Retrieve regime and calculate thermal time step using gradient heuristic
        _, _, current_regime = get_load_func(current_t_therm - start_t, 4.0 * cosim.rows_series, 0.5)
        
        curr_alpha_therm = haskey(opt_alpha_therm, current_regime) ? opt_alpha_therm[current_regime] : 1.0
        curr_dt_max_therm = haskey(opt_dt_max_therm, current_regime) ? opt_dt_max_therm[current_regime] : 200.0

        dT_dt = abs(T_max_C - prev_step_T_max) / prev_dt_therm
        epsilon_therm = 1e-4 
        
        if curr_alpha_therm == 0.0
            current_dt_therm = curr_dt_max_therm
        else
            proposed_dt_therm = curr_alpha_therm / (dT_dt + epsilon_therm)
            current_dt_therm = clamp(proposed_dt_therm, 0.1, curr_dt_max_therm)
        end
        
        t_target_therm = current_t_therm + current_dt_therm
        if t_target_therm > abs_total_time t_target_therm = abs_total_time; current_dt_therm = t_target_therm - current_t_therm end

        # Predict future load requirement and evaluate active cooling demands
        window = get_lookahead(tms_strategy, current_regime)
        if window > 0.0
            future_I, _, _ = get_load_func_future((current_t_therm + window) - start_t, 4.0 * cosim.rows_series, 0.5)
            cell_future_load = abs(future_I) / cosim.cols_parallel
        else
            cell_future_load = 0.0
        end
        
        safe_m_passive = max(m_passive, 1e-6)

        if m_active == 0.0 && safe_m_passive == 1e-6
            target_flow = safe_m_passive
            target_T = therm_params.inlet_temperature
        else
            current_flow = cosim.therm_integrator.ps[cosim.m_flow_sym]
            target_flow, dynamic_T = evaluate_tms_state(tms_strategy, current_regime, T_max_C, cell_future_load, current_flow, safe_m_passive, m_active)
            target_T = (target_flow >= (m_active * 0.99) && m_active > 0.0) ? dynamic_T : therm_params.inlet_temperature
        end
        
        cosim.therm_integrator.ps[cosim.m_flow_sym] = target_flow
        cosim.therm_integrator.ps[cosim.T_inlet_sym] = target_T
        
        # Step electrical cells forward to catch up with thermal target
        current_t_elec = current_t_therm
        while current_t_elec < t_target_therm
            
            if t_target_therm - current_t_elec < 1e-6 break end
            
            # Assess cell limits to safely cap incoming current target
            min_cell_v = Inf; max_cell_v = -Inf; min_soc = Inf; max_soc = -Inf
            for i in 1:num_cells
                v_c = cosim.cell_integrators[i][cosim.sys_cell_template.cell.v]
                soc_c = cosim.cell_integrators[i][cosim.sys_cell_template.cell.soc]
                if v_c < min_cell_v min_cell_v = v_c end
                if v_c > max_cell_v max_cell_v = v_c end
                if soc_c < min_soc min_soc = soc_c end
                if soc_c > max_soc max_soc = soc_c end
            end
            
            approx_pack_v = length(history_pack_v) > 0 ? history_pack_v[end] : (4.0 * cosim.rows_series)
            intended_I, _, _ = get_load_func(current_t_elec - start_t, (4.0 * cosim.rows_series), 0.5)
            
            worst_v = intended_I > 0 ? (max_cell_v * cosim.rows_series) : (intended_I < 0 ? (min_cell_v * cosim.rows_series) : approx_pack_v)
            worst_soc = intended_I > 0 ? max_soc : min_soc
            
            I_total, is_done, regime = get_load_func(current_t_elec - start_t, worst_v, worst_soc)
            if is_done break end
            
            # Periodically update internal resistance degradation factor to derate battery management limits
            if current_t_elec - last_R_sei_update_t >= R_sei_update_interval
                worst_derate_factor = 1.0
                
                for i in 1:num_cells
                    if !isnothing(R_sei_sym) current_R_sei[i] = cosim.cell_integrators[i][R_sei_sym]
                    elseif !isnothing(L_sei_sym) current_R_sei[i] = cosim.cell_integrators[i][L_sei_sym] / 5e-6
                    else current_R_sei[i] = 0.002 end
                    
                    cell_soc = clamp(cosim.cell_integrators[i][cosim.sys_cell_template.cell.soc], 0.0, 1.0)
                    R_fresh = dcir_interp(cell_soc)
                    R_aged = R_fresh + current_R_sei[i]
                    
                    derate = (R_fresh / R_aged)^derate_exponent
                    if derate < worst_derate_factor
                        worst_derate_factor = derate
                    end
                end
                
                dyn_chg_limit = MAX_CHG_BASE * worst_derate_factor 
                dyn_dsg_limit = MAX_DSG_BASE * worst_derate_factor
                
                last_R_sei_update_t = current_t_elec
            end
            
            I_total = clamp(I_total, dyn_dsg_limit, dyn_chg_limit)
            
            v_cell_worst = worst_v / cosim.rows_series
            taper_top_start = v_nom_max - 0.05 
            taper_bot_start = v_nom_min + 0.05 
            
            if I_total > 0.0 
                if v_cell_worst > taper_top_start
                    taper_factor = (v_nom_max - v_cell_worst) / (v_nom_max - taper_top_start)
                    I_total *= clamp(taper_factor, 0.0, 1.0)
                end
            elseif I_total < 0.0 
                if v_cell_worst < taper_bot_start
                    taper_factor = (v_cell_worst - v_nom_min) / (taper_bot_start - v_nom_min)
                    I_total *= clamp(taper_factor, 0.0, 1.0)
                end
            end
            
            # Fetch branch conductances to solve parallel nodal circuit
            for i in 1:num_cells
                V_term = cosim.cell_integrators[i][cosim.sys_cell_template.cell.v]
                I_prev = cosim.cell_integrators[i].ps[cosim.I_sym]
                cell_soc = clamp(cosim.cell_integrators[i][cosim.sys_cell_template.cell.soc], 0.0, 1.0)
                R_virt = dcir_interp(cell_soc) + current_R_sei[i]
                G_branches[i] = 1.0 / R_virt
                V_cores[i] = V_term - (I_prev * R_virt)
            end

            # Distribute currents based on resistance parameters for numerical stability
            pack_voltage = solve_pack_currents!(I_branches, V_cores, G_branches, cosim.rows_series, cosim.cols_parallel, I_total)
            
            # Conditionally flush metrics to history buffers based on temporal activity
            time_since_save = current_t_elec - last_save_t
            is_transient   = abs(I_total - last_save_I) > 1e-3 || abs(target_flow - last_save_flow) > 1e-6
            is_fast_moving = abs(pack_voltage - last_save_V) > 0.01 || abs(T_max_C - last_save_T) > 0.05
            is_active      = abs(I_total) > 1e-3 
            
            should_log = false
            if dense_logging
                should_log = true
            elseif ultra_sparse_logging
                should_log = time_since_save >= 3600.0 
            elseif sparse_logging
                should_log = (is_active && time_since_save >= 60.0) || (!is_active && time_since_save >= 600.0)
            else
                should_log = is_transient || (is_fast_moving && time_since_save >= 0.2) || (is_active && time_since_save >= 1.0) || (!is_active && time_since_save >= 10.0)
            end
            
            if should_log
                buffer_idx += 1

                T_slope = (T_max_C - prev_step_T_max) / (prev_dt_therm + 1e-6)
                T_continuous = T_max_C + T_slope * (current_t_elec - current_t_therm)
                
                push!(history_t, current_t_elec)
                push!(history_pack_v, pack_voltage)
                push!(history_T_max, T_continuous)
                push!(history_vel, target_flow / (rho * A_c))
                push!(history_I, I_total)
                push!(history_dt, dt_elec)
                push!(history_limit_chg, dyn_chg_limit)
                push!(history_limit_dsg, dyn_dsg_limit)
                
                pack_soc_val = clamp(cosim.cell_integrators[1][cosim.sys_cell_template.cell.soc], 0.0, 1.0) * 100.0
                push!(history_pack_soc, pack_soc_val)
                
                T_core_slopes = [(T_cores[i] - prev_T_cores[i]) / (prev_dt_therm + 1e-6) for i in 1:num_cells]
                
                for i in 1:num_cells 
                    history_V_cells[i, buffer_idx] = cosim.cell_integrators[i][cosim.sys_cell_template.cell.v]
                    history_SoH[i, buffer_idx] = !isnothing(SoH_sym) ? cosim.cell_integrators[i][SoH_sym] : 100.0
                    history_I_branches[i, buffer_idx] = I_branches[i]
                    history_T_cores[i, buffer_idx] = (T_cores[i] + T_core_slopes[i] * (current_t_elec - current_t_therm)) - 273.15
                end
                
                last_save_t = current_t_elec
                last_save_I = I_total
                last_save_V = pack_voltage
                last_save_T = T_max_C
                last_save_flow = target_flow
                
                if time() - last_print_wall_time >= 0.5 && verbose
                    elapsed_total = round(time() - start_wall_time, digits=1)
                    if total_time > 0.0
                        pct = clamp(((current_t_elec - start_t) / total_time) * 100.0, 0.0, 100.0)
                        filled = round(Int, 30 * (pct / 100.0))
                        bar = "[" * repeat("=", filled) * repeat(" ", 30 - filled) * "]"
                        print("\r$(print_prefix)$(bar) $(round(pct, digits=1))% | dt: $(round(dt_elec, digits=2))s | Max T: $(round(T_max_C, digits=2))°C | Wall: $(elapsed_total)s        ")
                    else
                        print("\r$(print_prefix)[Running] dt: $(round(dt_elec, digits=2))s | Max T: $(round(T_max_C, digits=2))°C | Wall: $(elapsed_total)s        ")
                    end
                    last_print_wall_time = time()
                end
                if buffer_idx >= chunk_limit buffer_idx = flush_chunk!(buffer_idx) end
            end

            # Determine next electrical time step respecting dynamics limits and physical boundaries
            signed_dV_dt = (pack_voltage - last_pack_v_elec) / dt_elec
            dV_dt = abs(signed_dV_dt)
            dI_dt = abs(I_total - last_I_total_elec) / dt_elec
            
            curr_w = haskey(opt_w_ratio, regime) ? opt_w_ratio[regime] : 2.0
            curr_alpha = haskey(opt_alpha_elec, regime) ? opt_alpha_elec[regime] : 10.0
            curr_dt_max = haskey(opt_dt_max_elec, regime) ? opt_dt_max_elec[regime] : 10.0
            
            if curr_alpha == 0.0
                dt_elec = curr_dt_max
            else
                metric = max(curr_w * (dV_dt / V_MAX), (dI_dt / I_MAX))
                dt_elec = clamp(curr_alpha / (metric + EPSILON_ELEC), DT_MIN_ELEC, curr_dt_max)
            end
            
            time_to_impact = Inf
            if signed_dV_dt < -1e-3
                time_to_impact = (pack_voltage - V_MIN_BMS) / abs(signed_dV_dt)
            elseif signed_dV_dt > 1e-3 
                time_to_impact = (V_MAX_BMS - pack_voltage) / signed_dV_dt
            end
            
            if time_to_impact > 0.0 && time_to_impact < dt_elec
                dt_elec = max(time_to_impact, DT_MIN_ELEC)
            end
            
            in_transient = false
            for t_b in boundaries
                if current_t_elec >= t_b - 1e-6 && current_t_elec < t_b + 1.0 - 1e-6
                    in_transient = true
                    break
                end
            end
            if in_transient dt_elec = min(dt_elec, 0.1) end
            
            for t_b in boundaries
                if t_b > current_t_elec + 1e-6
                    dt_elec = min(dt_elec, t_b - current_t_elec)
                    break
                end
            end
            if current_t_elec + dt_elec > t_target_therm dt_elec = t_target_therm - current_t_elec end
            
            if dt_elec < 1e-6 dt_elec = 1e-6 end
            
            last_pack_v_elec = pack_voltage
            last_I_total_elec = I_total

            # Predict load at next step and project target currents onto constituent cells
            next_I_total, _, _ = get_load_func(current_t_elec + dt_elec - start_t, worst_v, worst_soc)
            next_I_total = clamp(next_I_total, dyn_dsg_limit, dyn_chg_limit)
            
            if next_I_total > 0.0
                if v_cell_worst > taper_top_start
                    taper_factor = (v_nom_max - v_cell_worst) / (v_nom_max - taper_top_start)
                    next_I_total *= clamp(taper_factor, 0.0, 1.0)
                end
            elseif next_I_total < 0.0
                if v_cell_worst < taper_bot_start
                    taper_factor = (v_cell_worst - v_nom_min) / (taper_bot_start - v_nom_min)
                    next_I_total *= clamp(taper_factor, 0.0, 1.0)
                end
            end
            
            dI_total_dt = (next_I_total - I_total) / dt_elec
            t_target_elec = current_t_elec + dt_elec
            
            Threads.@threads for i in 1:num_cells 
                branch_fraction = abs(I_total) > 1e-6 ? (I_branches[i] / I_total) : (1.0 / cosim.cols_parallel)
                dI_branch_dt = dI_total_dt * branch_fraction
                
                int_c = cosim.cell_integrators[i]
                int_c.ps[cosim.I_sym] = I_branches[i] * cap_multipliers[i]
                int_c.ps[cosim.dI_dt_sym] = dI_branch_dt * cap_multipliers[i]
                int_c.ps[cosim.t_start_sym] = current_t_elec
                int_c.ps[cosim.T_sym] = T_cores[i] 
                
                SciMLBase.add_tstop!(int_c, t_target_elec)
                while int_c.t < t_target_elec
                    SciMLBase.step!(int_c)
                    if int_c.sol.retcode != SciMLBase.ReturnCode.Success && int_c.sol.retcode != SciMLBase.ReturnCode.Default
                        break
                    end
                end
            end
            
            fatal_crash = false
            for i in 1:num_cells
                rc = cosim.cell_integrators[i].sol.retcode
                if rc != SciMLBase.ReturnCode.Success && rc != SciMLBase.ReturnCode.Default
                    fatal_crash = true
                end
            end
            if fatal_crash
                is_done = true
                break
            end
            
            current_t_elec += dt_elec
        end 
        if is_done break end
        
        # Map cell heat generation inputs onto thermal model boundaries maintaining consistency during isothermal runs
        for i in 1:num_cells
            raw_Q = cosim.cell_integrators[i][cosim.sys_cell_template.cell.Q_total] * heat_multiplier
            Q_gens[i] = is_isothermal ? 0.0 : raw_Q
            
            dQ_dt_val = (Q_gens[i] - last_Q_gens[i]) / current_dt_therm
            cosim.therm_integrator.ps[cosim.Q_syms[i]] = last_Q_gens[i]
            cosim.therm_integrator.ps[cosim.dQ_dt_syms[i]] = dQ_dt_val
        end
        cosim.therm_integrator.ps[cosim.t_start_therm_sym] = current_t_therm
        
        SciMLBase.add_tstop!(cosim.therm_integrator, t_target_therm)
        while cosim.therm_integrator.t < t_target_therm
            SciMLBase.step!(cosim.therm_integrator)
        end
        
        for i in 1:num_cells
            last_Q_gens[i] = Q_gens[i]
            prev_T_cores[i] = T_cores[i]
            core_T_var = getproperty(getproperty(cosim.sys_therm.therm_pack_base, Symbol("cell_$i")), :core_cap).T
            T_cores[i] = is_isothermal ? therm_params.ambient_temperature : cosim.therm_integrator[core_T_var]
        end
        
        prev_step_T_max = T_max_C
        T_max_C = maximum(T_cores) - 273.15
        prev_dt_therm = current_dt_therm
        current_t_therm = t_target_therm
    end 
    
    # Conclude operations writing final chunk outputs mapping everything into singular master dataframe
    elapsed_total = round(time() - start_wall_time, digits=1)
    if verbose
        pct = is_done ? 100.0 : 100.0
        filled = 30
        bar = "[" * repeat("=", filled) * repeat(" ", 30 - filled) * "]"
        print("\r$(print_prefix)$(bar) $(round(pct, digits=1))% | dt: $(round(dt_elec, digits=2))s | Max T: $(round(T_max_C, digits=2))°C | Wall: $(elapsed_total)s        \n")
    end
    
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