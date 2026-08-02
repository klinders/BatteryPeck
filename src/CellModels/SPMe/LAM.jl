module LAM

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    StressDriven(; name, p::SideReactionParameters, s::SolidParticleParameters, V, g)

Create a Loss of Active Material (LAM) model driven by mechanical stress.

Models the isolation and detachment of active material particles resulting from 
swelling-induced mechanical stress during intercalation. Computes hydrostatic stress
and isolates the tensile component to drive the LAM reaction rate.

# Output Variables
- `j_lam`: Rate of active material loss
- `lli`: Loss of lithium inventory (trapped in isolated particles)
- `Q_loss`: Capacity loss (Ah)
"""
function StressDriven(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named σᵣ = RealInput()
    @named σₜ = RealInput()
    @named c_r = RealInput()

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        j_lam(t), [guess=0.0]
        lli(t) = 0.0
        Q_loss(t), [guess=0.0]
    end

    # obtain the rate of loss of active materials (LAM) by stress
    # compute hydrostatic stress
    σₕ = (σᵣ.u + 2 * σₜ.u) / 3

    # separate compressive and tensile stresses
    σₕ_c = ifelse(σₕ < 0 , σₕ, 0.0)
    σₕ_t = ifelse(σₕ > 0 , σₕ, 0.0)
    
    # assuming that only tensile stress contributes and that the minimum
    # (tensile) hydrostatic stress is zero for full cycles
    σₕ_min = σₕ * 0.0

    eqns = [
        j_lam ~ -s.β_LAM*((σₕ_t - σₕ_min) / s.stress_critical)^s.m_LAM,
        Dt(lli) ~ -V * c_r.u * j_lam,
        Q_loss ~ lli * F / 3600
    ]

    System(eqns, t; name=name, systems=[σᵣ, σₜ, c_r])
end

end