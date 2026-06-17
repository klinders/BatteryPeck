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
export TMSStrategy, ReactiveTMS, AnticipativeTMS, evaluate_tms_state, get_lookahead

"""
    TMSStrategy

Abstract base type for thermal management system control logic.
"""
abstract type TMSStrategy end

"""
    velocity_to_mass_flow(v::Float64, rho::Float64, A_c::Float64)

Convert coolant velocity to mass flow rate.

Multiplies velocity by fluid density and cross sectional area.

# Arguments
- `v::Float64`: Coolant flow velocity
- `rho::Float64`: Coolant fluid density
- `A_c::Float64`: Cross sectional area of cooling channel

# Returns
- Calculated mass flow rate
"""
function velocity_to_mass_flow(v::Float64, rho::Float64, A_c::Float64)
    return rho * A_c * v
end

"""
    get_future_load_avg(current_t::Float64, get_load_func::Function, window_size::Float64, v_pack::Float64, soc::Float64, samples::Int=10)

Calculate average expected current load over future time window.

Iterates through future time steps and queries load function. Averages absolute current magnitude across sampled points.

# Arguments
- `current_t::Float64`: Current simulation time
- `get_load_func::Function`: Function returning expected load current
- `window_size::Float64`: Duration to look ahead
- `v_pack::Float64`: Current pack voltage
- `soc::Float64`: Current state of charge
- `samples::Int`: Number of points to sample within window

# Returns
- Average absolute current load
"""
function get_future_load_avg(current_t::Float64, get_load_func::Function, window_size::Float64, v_pack::Float64, soc::Float64, samples::Int=10)
    # Initialise load accumulator and calculate time step size
    total_load = 0.0
    dt = window_size / samples
    
    # Aggregate absolute load values across future sample points
    for i in 1:samples
        total_load += abs(get_load_func(current_t + i*dt, v_pack, soc)[1])
    end
    
    return total_load / samples
end

"""
    ReactiveTMS(; T_high=35.0, T_low=32.0, T_coolant_passive=298.15, T_coolant_active=288.15)

Define reactive thermal management strategy parameters.

Triggers active cooling only when maximum temperature exceeds high threshold. Reverts to passive cooling when temperature drops below low threshold.
"""
Base.@kwdef struct ReactiveTMS <: TMSStrategy
    T_high::Float64 = 35.0
    T_low::Float64  = 32.0
    T_coolant_passive::Float64 = 298.15
    T_coolant_active::Float64  = 288.15
end

"""
    AnticipativeTMS(; T_high=35.0, T_low=32.0, cell_load_threshold=7.5, lookahead_window=300.0, T_coolant_passive=298.15, T_coolant_active=288.15)

Define anticipative thermal management strategy parameters.

Triggers active cooling preemptively based on expected future load or when temperature exceeds high threshold. 
"""
Base.@kwdef struct AnticipativeTMS <: TMSStrategy
    T_high::Float64 = 35.0
    T_low::Float64  = 32.0
    cell_load_threshold::Float64 = 7.5 
    lookahead_window::Float64 = 300.0  
    T_coolant_passive::Float64 = 298.15
    T_coolant_active::Float64  = 288.15
end

"""
    get_lookahead(strategy)

Retrieve lookahead window duration for specific thermal management strategy.

Returns zero for reactive strategies to prevent unnecessary CPU load. Returns defined window size for anticipative strategies.

# Arguments
- `strategy`: Configured thermal management strategy object

# Returns
- Lookahead time window in seconds
"""
get_lookahead(::ReactiveTMS) = 0.0
get_lookahead(s::AnticipativeTMS) = s.lookahead_window

"""
    evaluate_tms_state(strategy::ReactiveTMS, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)

Determine active or passive cooling state based on reactive logic.

Evaluates current maximum temperature against thresholds to set hysteresis state.

# Arguments
- `strategy::ReactiveTMS`: Reactive management configuration
- `T_max_C::Float64`: Current maximum cell temperature
- `cell_future_load::Float64`: Average expected future load (unused)
- `current_flow::Float64`: Current coolant mass flow rate
- `m_passive::Float64`: Passive resting mass flow rate
- `m_active::Float64`: Active cooling mass flow rate

# Returns
- Tuple containing target mass flow and target coolant temperature
"""
function evaluate_tms_state(strategy::ReactiveTMS, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)
    # Engage active cooling if maximum temperature exceeds high threshold
    if T_max_C > strategy.T_high
        return m_active, strategy.T_coolant_active
    # Revert to passive cooling if maximum temperature drops below low threshold
    elseif T_max_C < strategy.T_low
        return m_passive, strategy.T_coolant_passive
    # Maintain current hysteresis state if temperature resides between bounds
    else
        if current_flow >= (m_active * 0.99)
            return m_active, strategy.T_coolant_active
        else
            return m_passive, strategy.T_coolant_passive
        end
    end
end

"""
    evaluate_tms_state(strategy::AnticipativeTMS, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)

Determine active or passive cooling state based on anticipative logic.

Evaluates current temperature and predictive load against configured thresholds.

# Arguments
- `strategy::AnticipativeTMS`: Anticipative management configuration
- `T_max_C::Float64`: Current maximum cell temperature
- `cell_future_load::Float64`: Average expected future load
- `current_flow::Float64`: Current coolant mass flow rate
- `m_passive::Float64`: Passive resting mass flow rate
- `m_active::Float64`: Active cooling mass flow rate

# Returns
- Tuple containing target mass flow and target coolant temperature
"""
function evaluate_tms_state(strategy::AnticipativeTMS, T_max_C::Float64, cell_future_load::Float64, current_flow::Float64, m_passive::Float64, m_active::Float64)
    # Evaluate temperature and load threshold conditions
    is_hot = T_max_C > strategy.T_high
    is_cool = T_max_C < strategy.T_low
    is_heavy_load = cell_future_load > strategy.cell_load_threshold

    # Engage active cooling if pack is hot or heavy future load is detected
    if is_hot || is_heavy_load
        return m_active, strategy.T_coolant_active
    # Revert to passive cooling if pack is cool and no heavy load is expected
    elseif is_cool && !is_heavy_load
        return m_passive, strategy.T_coolant_passive
    # Maintain current hysteresis state if conditions fall between defined triggers
    else
        if current_flow >= (m_active * 0.99) 
            return m_active, strategy.T_coolant_active
        else
            return m_passive, strategy.T_coolant_passive
        end
    end
end