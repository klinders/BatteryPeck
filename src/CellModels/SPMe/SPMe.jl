# =====================================================================================================================
# SPMe.jl
#
# SPMe implementation
# Source [1]: https://doi.org/10.1016/j.apm.2022.12.009       # SPMe model
# Source [2]: https://doi.org/10.1016/j.electacta.2022.140700 # Regan2022 stoichiometry vs. SoC
# =====================================================================================================================

# Import packages
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")
include("SEI.jl")
include("LithiumPlating.jl")

import NaNMath

# Polynomial fit parameters entropic term positive electrode
const P_coeff_pos = (
    a1 = 0.04006, b1 = 0.2828, c1 = 0.0009855,
    a2 = -0.06656, b2 = 0.8032, c2 = 0.02179
)

# Polynomial fit parameters entropic term negative electrode
const P_coeff_neg = (
    a0 = -0.111, b0 = 0.02901,
    a1 = 0.3562, b1 = 0.08308, c1 = 0.004621
)

# Entropic term positive electrode
function dUp_dT_f(z)
    p = P_coeff_pos
    val_mV = p.a1 * exp(-((z - p.b1)^2) / p.c1) + p.a2 * exp(-((z - p.b2)^2) / p.c2)
    return val_mV * 1e-3
end

# Entropic term negative electrode
function dUn_dT_f(z)
    p = P_coeff_neg
    val_mV = p.a0 * z + p.b0 + p.a1 * exp(-((z - p.b1)^2) / p.c1)
    return val_mV * 1e-3
end

# Initiate end of experiment
function abort!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    t = round(int.t,digits=2)
    @warn "Simulation step terminated at t=$t"
    return (;)
end

"""
The SPMe model implemented based on [MarquisEtAl2019](@citet) and [BrosaPlanellaWidanage2023](@citet)
"""
function SPMe(; name="SPMe", params::BatteryParameters, Q=0, N=Dict(:Nₓ=>[10,10,10], :Nᵣ=>[10,10]), side_reactions=true)
    @parameters begin
        t # Time variable
    end
    Dt = Differential(t)

    g = build_fvm_geometry(params, N)
    
    # Electrical ports
    @named p = Pin()
    @named n = Pin()
    @named T = RealInput(guess=298)

    @named pe = SolidParticle(p=params.p, g=g.pe)
    @named ne = SolidParticle(p=params.n, g=g.ne)
    @named el = Electrolyte(p=params.e, g=g.el)
    
    submodels = Any[p, n, T, pe, ne, el]
    
    if side_reactions
        # Using Branch 2's advanced degradation models
        @named sei = SolventDiffusionLimitedSEI(p=params.n.side_reactions[1], s=params.n, g=g) 
        @named plating = PartiallyReversiblePlating(p=params.n.side_reactions[1], s=params.n, g=g)
        push!(submodels, sei, plating)
    end

    @variables begin
        # Terminal voltage and current
        v(t), [guess=4.0]
        i(t)
        soc(t)

        U₀(t)
        ηᵣ(t)
        (ηₙ(t))[1:g.el.Nx[1]]
        (ηₚ(t))[1:g.el.Nx[3]]
        Δϕₛ(t), [guess=0]
        (ϕₙ(t))[1:g.el.Nx[1]]
        (ϕₚ(t))[1:g.el.Nx[3]]
        (jₙ0(t))[1:g.el.Nx[1]]
        (jₚ0(t))[1:g.el.Nx[3]]
        j̄ₙ0(t)
        j̄ₚ0(t)
    
        ϕ̄ₙ(t)
        ϕ̄ₚ(t)
        aₙ(t), [guess=3*(1-params.e.ϵₙ)/params.n.Rₖ]
        aₚ(t)

        Rᵢ(t)

        # Heat generation sources
        Qᵢ(t)       # Reaction heat
        Qₛ(t)       # Solid Ohmic heat
        Qf(t)       # Film Ohmic heat
        Q_rev(t)    # Reversible heat
        Q_total(t)  # Total heat
        Q_sei(t)
        Q_plating(t)
        
        # Degradation trackers
        Q_loss_Ah(t)
        SoH(t)
    end

    if Q==0
        Q = params.Q₀
    end
    
    # Scale the current density to the electrode area
    A = params.Hcc*params.Wcc*params.n_el*(Q/params.Q₀)
    i_app = i/A

    # X-average
    x = g.el.x_centers
    L = sum(g.el.Ls)
    
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    Nn = length(g.el.ixₙ)
    Np = length(g.el.ixₚ)
    
    ηᵣn = 2*R*T.u/F*asinh(ne.J.u / (2*j̄ₙ0))
    ηᵣp = 2*R*T.u/F*asinh(pe.J.u / (2*j̄ₚ0))

    # X-average of the electrolyte potential
    ϕₛ_n = [i_app*(x[i] - 2*params.e.Lₙ)*x[i]/2/params.n.σₖ/params.e.Lₙ for i in g.el.ixₙ]
    ϕₛ_p = [i_app*(x[i] + (x[i] - L)^2/(2*params.e.Lₚ))/params.p.σₖ for i in g.el.ixₚ]

    eqns = Equation[
        # Temperatures
        el.T.u ~ T.u,
        pe.T.u ~ T.u,
        ne.T.u ~ T.u,
        soc ~ ne.z,

        # Potentials
        U₀ ~ pe.U₀ - ne.U₀,
        ηᵣ ~ ηᵣp - ηᵣn,
        Δϕₛ ~ -i_app/3*(params.e.Lₚ/params.p.σₖ + params.e.Lₙ/params.n.σₖ),

        # Exchange current densities
        [jₙ0[i] ~ params.n.mₖ*sqrt(el.cₑ[g.el.ixₙ[i]]*ne.c_surf*(params.n.c₊-ne.c_surf)) for i in 1:Nn]...,
        [jₚ0[i] ~ params.p.mₖ*sqrt(el.cₑ[g.el.ixₚ[i]]*pe.c_surf*(params.p.c₊-pe.c_surf)) for i in 1:Np]...,
        j̄ₙ0 ~ sum(jₙ0)/Nn,
        j̄ₚ0 ~ sum(jₚ0)/Np,

        [ϕₙ[i] ~ n.v + ϕₛ_n[i] for i in 1:Nn]...,
        [ϕₚ[i] ~ p.v - ϕₛ_p[i] for i in 1:Np]...,
        [ηₙ[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
        [ηₚ[i] ~ ϕₚ[i] - el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        ϕ̄ₙ ~ sum(ϕₙ)/Nn,
        ϕ̄ₚ ~ sum(ϕₚ)/Np,
        
        # Apply conditional constraints for side reactions
        v ~ U₀ + ηᵣ + el.ηₑ + el.Δϕₑ + Δϕₛ + (side_reactions ? sei.ϕf_x : 0.0),
        el.Δϕₙ.u ~ ne.U₀ + ηᵣn - (side_reactions ? sei.ϕf_x : 0.0),
        
        Rᵢ ~ (U₀-v)/i, 
        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        el.i_app.u ~ i_app,
        el.ϕₛn.u ~ ϕ̄ₙ,

        # Current density in the negative electrode with conditional side reactions
        ne.J.u ~ (i_app/params.e.Lₙ)/aₙ - (side_reactions ? sei.j_sei_x : 0.0) - (side_reactions ? plating.j_stripping_x : 0.0),
        pe.J.u ~  -i_app/params.e.Lₚ/aₚ, 
        
        aₙ ~ 3*(1-el.ϵ̄ₙ)/params.n.Rₖ,
        aₚ ~ params.p.aₖ,

        # Porosity (assumed constant, conditional SEI growth)
        [el.ϵ[i] ~ params.e.ϵₙ - (side_reactions ? (aₙ*(sei.L_sei[i] - params.n.L_sei₀)) : 0.0) for i in g.el.ixₙ]...,
        [el.ϵ[i] ~ params.e.ϵₛ for i in g.el.ixₛ]...,
        [el.ϵ[i] ~ params.e.ϵₚ for i in g.el.ixₚ]...,

        # Heat sources
        Qᵢ ~ -i_app * ηᵣ / L,
        Qₛ ~ -i_app * Δϕₛ / L,
        Qf ~ -i_app * (side_reactions ? sei.ϕf_x : 0.0) / L,
        Q_rev ~ (i_app / L) * T.u * (dUn_dT_f(ne.z) - dUp_dT_f(pe.z)),
        Q_sei ~ -(side_reactions ? sei.j_sei_x * sei.ϕf_x : 0.0) / L,
        Q_plating ~ -(side_reactions ? plating.j_stripping_x * plating.ϕf_x : 0.0) / L,
        Q_total ~ el.Qₑ + Qᵢ + Qₛ + Qf + Q_rev + Q_sei + Q_plating,
        
        SoH ~ 100.0 * (1.0 - (Q_loss_Ah / params.Q₀))
    ]

    # Only bind SEI boundaries if they are active
    if side_reactions
        Volume_n = params.e.Lₙ * A
        
        # Calculate initial SEI concentration to offset Q_loss at t=0
        c_sei_0 = params.n.side_reactions[1].Lf₀ / params.n.side_reactions[1].V̄ * params.n.aₖ
        
        # Moles of Lithium = Concentration (mol/m^3) * Volume (m^3) * stoichiometric ratio (z)
        moles_Li_sei = (sei.c_sei_x - c_sei_0) * Volume_n * params.n.side_reactions[1].z
        moles_Li_plating = (plating.c_dead_x + plating.c_plating_x) * Volume_n # z=1 for pure Li plating
        
        append!(eqns, [
            # Convert lost moles of Lithium back into lost Ampere-hours
            Q_loss_Ah ~ (moles_Li_sei + moles_Li_plating) * 96485.0 / 3600.0,
            
            sei.J.u ~ ne.J.u, 
            sei.T.u ~ T.u,
            sei.aₖ.u ~ aₙ,
            [sei.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[i] for i in 1:Nn]...,

            plating.J.u ~ ne.J.u,
            plating.T.u ~ T.u,
            plating.aₖ.u ~ aₙ,
            [plating.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[i] for i in 1:Nn]...,
            plating.η_sei.u ~ sei.ϕf,
            [plating.cₑ.u[i] ~ params.e.cₜ for i in 1:Nn]...,
            
            # Wire L_sei into the plating model for coupled structural damage
            [plating.L_sei.u[i] ~ sei.L_sei[i] for i in 1:Nn]...
        ])
    else
        append!(eqns, [Q_loss_Ah ~ 0.0])
    end

    events = [
        [
            v ~ params.Vmin,
            v ~ params.Vmax,
            pe.c_surf ~ params.p.c₊*0.999,
            ne.c_surf ~ params.n.c₊*0.999,
            pe.c_surf ~ params.p.c₊*0.001,
            ne.c_surf ~ params.n.c₊*0.001
        ]=>(abort!,(;))
    ]

    return System(eqns, t; name=name, systems=submodels, continuous_events=events)
end