# ==============================================================================
# TMSControl.jl
# Anticipative and reactive control logic for thermal management system
# ==============================================================================
export velocity_to_mass_flow, get_max_cell_temp, get_future_power_avg
export build_baseline_callback, build_anticipative_callback

using ModelingToolkit
using SciMLBase
using DiffEqCallbacks

"""
    velocity_to_mass_flow(v, rho, A_c)

Convert coolant velocity to mass flow rate.
"""
function velocity_to_mass_flow(v::Float64, rho::Float64, A_c::Float64)
    return rho * A_c * v
end

"""
    get_max_cell_temp(integrator, thermal_sys, num_cells)

Extract maximum core temperature across all cells.
"""
function get_max_cell_temp(integrator, thermal_sys, num_cells::Int)
    temps = [integrator[getproperty(thermal_sys, Symbol("cell_$i")).core_cap.T] for i in 1:num_cells]
    return maximum(temps)
end

"""
    get_future_power_avg(t_current, steps, window_size)

Calculate exact time-weighted average of future load using step periods.
"""
function get_future_power_avg(t_current::Float64, steps, window_size::Float64)
    t_end_window = t_current + window_size
    total_energy = 0.0
    
    t_step_start = 0.0
    for step in steps
        t_step_end = t_step_start + step.period
        
        if t_step_end > t_current && t_step_start < t_end_window
            overlap_start = max(t_current, t_step_start)
            overlap_end = min(t_end_window, t_step_end)
            overlap_duration = overlap_end - overlap_start
            
            val = 0.0
            if hasproperty(step, :value)
                val = abs(step.value)
            end
            
            total_energy += val * overlap_duration
        end
        
        t_step_start = t_step_end
        if t_step_start >= t_end_window
            break
        end
    end
    return total_energy / window_size
end

"""
    build_baseline_callback(thermal_sys, num_cells, m_flow_active, m_flow_passive, T_active, T_passive)

Construct discrete callback for reactive thermal management.
"""
function build_baseline_callback(thermal_sys, num_cells::Int, m_flow_active::Float64, m_flow_passive::Float64, T_active::Float64=288.15, T_passive::Float64=298.15)
    
    # Reference for mechanical cooldown to prevent microsecond double-triggers
    last_trigger = Ref(-100.0) 
    
    function condition(u, t, integrator)
        if t - last_trigger[] < 5.0 return false end
        
        T_max = get_max_cell_temp(integrator, thermal_sys, num_cells)
        current_flow = integrator.ps[thermal_sys.fluid_inlet.m_flow_in]
        
        if T_max > 308.15 && abs(current_flow - m_flow_active) > 1e-6
            return true
        elseif T_max < 305.15 && abs(current_flow - m_flow_passive) > 1e-6
            return true
        end
        return false
    end

    function affect!(integrator)
        T_max = get_max_cell_temp(integrator, thermal_sys, num_cells)
        current_flow = integrator.ps[thermal_sys.fluid_inlet.m_flow_in]
        
        if T_max > 308.15 && abs(current_flow - m_flow_active) > 1e-6
            integrator.ps[thermal_sys.fluid_inlet.m_flow_in] = m_flow_active
            integrator.ps[thermal_sys.fluid_inlet.T_inlet] = T_active
            last_trigger[] = integrator.t 
            
            # Safely update physics without resetting DAE states
            u_modified!(integrator, false) 
            println("  [TMS] ACTIVE at t = $(round(integrator.t, digits=1))s | T_max = $(round(T_max - 273.15, digits=2))C")
            
        elseif T_max < 305.15 && abs(current_flow - m_flow_passive) > 1e-6
            integrator.ps[thermal_sys.fluid_inlet.m_flow_in] = m_flow_passive
            integrator.ps[thermal_sys.fluid_inlet.T_inlet] = T_passive
            last_trigger[] = integrator.t 
            
            # Safely update physics without resetting DAE states
            u_modified!(integrator, false) 
            println("  [TMS] PASSIVE at t = $(round(integrator.t, digits=1))s | T_max = $(round(T_max - 273.15, digits=2))C")
        end
    end

    return DiscreteCallback(condition, affect!; save_positions=(false, true))
end

"""
    build_anticipative_callback(thermal_sys, num_cells, exp_steps, m_flow_active, m_flow_passive, power_threshold, window_size, T_active, T_passive)

Construct discrete callback for predictive thermal management.
"""
function build_anticipative_callback(thermal_sys, num_cells::Int, exp_steps, m_flow_active::Float64, m_flow_passive::Float64, power_threshold::Float64, window_size::Float64=300.0, T_active::Float64=288.15, T_passive::Float64=298.15)
    
    # Reference for mechanical cooldown to prevent microsecond double-triggers
    last_trigger = Ref(-100.0)
    
    function condition(u, t, integrator)
        if t - last_trigger[] < 5.0 return false end
        
        T_max = get_max_cell_temp(integrator, thermal_sys, num_cells)
        P_future_avg = get_future_power_avg(integrator.t, exp_steps, window_size)
        current_flow = integrator.ps[thermal_sys.fluid_inlet.m_flow_in]
        
        if (P_future_avg > power_threshold || T_max > 308.15) && abs(current_flow - m_flow_active) > 1e-6
            return true
        elseif T_max < 305.15 && P_future_avg <= power_threshold && abs(current_flow - m_flow_passive) > 1e-6
            return true
        end
        return false
    end

    function affect!(integrator)
        T_max = get_max_cell_temp(integrator, thermal_sys, num_cells)
        P_future_avg = get_future_power_avg(integrator.t, exp_steps, window_size)
        current_flow = integrator.ps[thermal_sys.fluid_inlet.m_flow_in]
        
        if (P_future_avg > power_threshold || T_max > 308.15) && abs(current_flow - m_flow_active) > 1e-6
            integrator.ps[thermal_sys.fluid_inlet.m_flow_in] = m_flow_active
            integrator.ps[thermal_sys.fluid_inlet.T_inlet] = T_active
            last_trigger[] = integrator.t
            
            # Safely update physics without resetting DAE states
            u_modified!(integrator, false) 
            println("  [TMS] Anticipative ACTIVE at t = $(round(integrator.t, digits=1))s")
            
        elseif T_max < 305.15 && P_future_avg <= power_threshold && abs(current_flow - m_flow_passive) > 1e-6
            integrator.ps[thermal_sys.fluid_inlet.m_flow_in] = m_flow_passive
            integrator.ps[thermal_sys.fluid_inlet.T_inlet] = T_passive
            last_trigger[] = integrator.t
            
            # Safely update physics without resetting DAE states
            u_modified!(integrator, false) 
            println("  [TMS] Anticipative PASSIVE at t = $(round(integrator.t, digits=1))s")
        end
    end

    return DiscreteCallback(condition, affect!; save_positions=(false, true))
end