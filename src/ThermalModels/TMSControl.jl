# ==============================================================================
# TMSControl.jl
# Pure Julia evaluation functions for Battery Thermal Management System logic
# Contains:
# 1. velocity_to_mass_flow: Convert coolant velocity to mass flow rate
# 2. get_future_load_avg: Calculate average expected current load over future time window
# 3. get_lookahead: Retrieve lookahead window duration for specific thermal management strategy
# 4. evaluate_tms_state: Determine active or passive cooling state based on current conditions
# ==============================================================================

export velocity_to_mass_flow, get_future_load_avg
export TMSStrategy, ReactiveTMS, AnticipativeTMS, HybridTMS, evaluate_tms_state, get_lookahead

abstract type TMSStrategy end

"""
    velocity_to_mass_flow(v::Float64, rho::Float64, A_c::Float64)

Convert coolant velocity to mass flow rate.

# Arguments
- `v::Float64`: Coolant velocity
- `rho::Float64`: Fluid density
- `A_c::Float64`: Cross sectional area

# Returns
- Mass flow rate
"""
function velocity_to_mass_flow(v::Float64, rho::Float64, A_c::Float64)
    return rho * A_c * v
end

"""
    get_future_load_avg(current_t::Float64, get_load_func::Function, window_size::Float64, v_pack::Float64, soc::Float64, samples::Int=10)

Calculate average expected current load over future time window.

# Arguments
- `current_t::Float64`: Current absolute simulation time
- `get_load_func::Function`: Experimental load closure
- `window_size::Float64`: Temporal duration to scan ahead
- `v_pack::Float64`: Latest evaluated pack voltage
- `soc::Float64`: Latest evaluated state of charge
- `samples::Int`: Number of discrete checks within window

# Returns
- Average load magnitude
"""
function get_future_load_avg(current_t::Float64, get_load_func::Function, window_size::Float64, v_pack::Float64, soc::Float64, samples::Int=10)
    # Aggregate load samples over future interval
    total_load = 0.0
    dt = window_size / samples
    for i in 1:samples
        total_load += abs(get_load_func(current_t + i*dt, v_pack, soc)[1])
    end
    return total_load / samples
end

# Dynamically link coolant temperatures to ambient baseline
Base.@kwdef struct ReactiveTMS <: TMSStrategy
    ambient_temp::Float64 = 298.15
    T_high::Float64 = 35.0
    T_low::Float64  = 32.0
    T_coolant_passive::Float64 = ambient_temp
    T_coolant_active::Float64  = ambient_temp - 10.0
end

Base.@kwdef struct AnticipativeTMS <: TMSStrategy
    ambient_temp::Float64 = 298.15
    T_high::Float64 = 35.0
    T_low::Float64  = 32.0
    cell_load_threshold::Float64 = 7.5  
    lookahead_window::Float64 = 300.0   
    T_coolant_passive::Float64 = ambient_temp
    T_coolant_active::Float64  = ambient_temp - 10.0
end

Base.@kwdef struct HybridTMS <: TMSStrategy
    reactive::ReactiveTMS
    anticipative::AnticipativeTMS
    regime_map::Dict{Symbol, Symbol} 
end

"""
    get_lookahead(strategy::TMSStrategy, regime::Symbol)

Retrieve lookahead window duration for specific thermal management strategy.

# Arguments
- `strategy::TMSStrategy`: Assigned thermal management struct
- `regime::Symbol`: Current operational regime mapped from experiment

# Returns
- Window duration
"""
get_lookahead(::ReactiveTMS, regime::Symbol) = 0.0
get_lookahead(s::AnticipativeTMS, regime::Symbol) = s.lookahead_window

function get_lookahead(s::HybridTMS, regime::Symbol)
    # Route lookahead retrieval to mapped internal strategy component
    mode = get(s.regime_map, regime, :reactive)
    if mode == :anticipative
        return s.anticipative.lookahead_window
    else
        return 0.0
    end
end

"""
    evaluate_tms_state(strategy::TMSStrategy, regime::Symbol, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)

Determine active or passive cooling state based on current conditions.

# Arguments
- `strategy::TMSStrategy`: Active thermal management struct
- `regime::Symbol`: Identifier dictating internal route map
- `T_max_C::Float64`: Peak tracked cell core temperature
- `cell_future_load::Float64`: Evaluated average lookahead current
- `current_flow::Float64`: Monitored active mass flow rate
- `m_passive::Float64`: Resting flow configuration
- `m_active::Float64`: Maximum pumping configuration

# Returns
- Tuple containing required flow rate and inlet temperature
"""
function evaluate_tms_state(strategy::ReactiveTMS, regime::Symbol, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)
    # Assess hysteresis thresholds against immediate thermal conditions
    if T_max_C > strategy.T_high
        return m_active, strategy.T_coolant_active
    elseif T_max_C < strategy.T_low
        return m_passive, strategy.T_coolant_passive
    else
        if current_flow >= (m_active * 0.99)
            return m_active, strategy.T_coolant_active
        else
            return m_passive, strategy.T_coolant_passive
        end
    end
end

function evaluate_tms_state(strategy::AnticipativeTMS, regime::Symbol, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)
    # Monitor predicted future loads alongside current thermal conditions
    is_hot = T_max_C > strategy.T_high
    is_cool = T_max_C < strategy.T_low
    is_heavy_load = cell_future_load > strategy.cell_load_threshold

    if is_hot || is_heavy_load
        return m_active, strategy.T_coolant_active
    elseif is_cool && !is_heavy_load
        return m_passive, strategy.T_coolant_passive
    else
        if current_flow >= (m_active * 0.99) 
            return m_active, strategy.T_coolant_active
        else
            return m_passive, strategy.T_coolant_passive
        end
    end
end

function evaluate_tms_state(strategy::HybridTMS, regime::Symbol, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)
    # Extract mapped regime key and invoke nested evaluation function
    mode = get(strategy.regime_map, regime, :reactive)
    if mode == :anticipative
        return evaluate_tms_state(strategy.anticipative, regime, T_max_C, cell_future_load, current_flow, m_passive, m_active)
    else
        return evaluate_tms_state(strategy.reactive, regime, T_max_C, cell_future_load, current_flow, m_passive, m_active)
    end
end