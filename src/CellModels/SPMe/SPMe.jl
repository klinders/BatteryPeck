# =====================================================================================================================
# SPMe.jl
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")
include("SEI.jl")
include("LithiumPlating.jl")
include("ParticleCracking.jl")
include("LAM.jl")

import NaNMath

const P_coeff_pos = (a1 = 0.04006, b1 = 0.2828, c1 = 0.0009855, a2 = -0.06656, b2 = 0.8032, c2 = 0.02179)
const P_coeff_neg = (a0 = -0.111, b0 = 0.02901, a1 = 0.3562, b1 = 0.08308, c1 = 0.004621)

function dUp_dT_f(z)
    p = P_coeff_pos
    val_mV = p.a1 * exp(-((z - p.b1)^2) / p.c1) + p.a2 * exp(-((z - p.b2)^2) / p.c2)
    return val_mV * 1e-3
end

function dUn_dT_f(z)
    p = P_coeff_neg
    val_mV = p.a0 * z + p.b0 + p.a1 * exp(-((z - p.b1)^2) / p.c1)
    return val_mV * 1e-3
end

# Diagnostic abort callbacks
function abort_vmin!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Voltage hit Vmin."
    return (;)
end
function abort_vmax!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Voltage hit Vmax."
    return (;)
end
function abort_pe_max!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Cathode surface concentration hit 99.99% (Saturated)."
    return (;)
end
function abort_ne_max!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Anode surface concentration hit 99.99% (Saturated)."
    return (;)
end
function abort_pe_min!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Cathode surface concentration hit 0.01% (Depleted)."
    return (;)
end
function abort_ne_min!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Anode surface concentration hit 0.01% (Depleted)."
    return (;)
end
function abort_ce_min!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Electrolyte concentration depleted (< 1.0 mol/m³)."
    return (;)
end
function abort_porosity_min!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    @warn "Simulation terminated at t=$(round(int.t,digits=2)): Anode porosity critically low (< 1%)."
    return (;)
end

"""
The SPMe model implemented based on [MarquisEtAl2019](@citet) and [BrosaPlanellaWidanage2023](@citet)

**Arguments**
- `name` (optional) defaults to SPMe
- `params ::BatteryParameters` Parameters from the given parameter set
- `Q ::Real` (optional) The capacity of the cell in Ah
- `N ::Dict(:Nₓ=>[::Int, ::Int, ::Int], :Nᵣ=>[::Int, ::Int])` (optional) number of mesh nodes in the particles and electrolyte
- `side_reactions ::Bool` (optional) enable side reactions
"""
function SPMe(; name="SPMe", params::BatteryParameters, Q=0, N=Dict(:Nₓ=>[10,10,10], :Nᵣ=>[10,10]), side_reactions=true)
    @parameters t 
    Dt = Differential(t)
    g = build_fvm_geometry(params, N)
    
    if Q==0
        Q = params.Q₀
    end

    Nn = length(g.el.ixₙ)
    Np = length(g.el.ixₚ)

    A = params.Hcc*params.Wcc*params.n_el*(Q/params.Q₀)
    
    # MOVED CONSTANTS UP HERE
    R = 8.314 
    F = 96485 

    @named p = Pin()
    @named n = Pin()
    @named T = RealInput(guess=298.15)

    @named pe = SolidParticle(p=params.p, g=g.pe)
    @named ne = SolidParticle(p=params.n, g=g.ne)
    @named el = Electrolyte(p=params.e, g=g.el)
    
    submodels = Any[p, n, T, pe, ne, el]
    
    if side_reactions
        @named sei = SolventDiffusionLimitedSEI(p=params.n.side_reactions[1], s=params.n, g=g) 
        @named plating = PartiallyReversiblePlating(p=params.n.side_reactions[1], s=params.n, g=g)
        @named cracking_n = ParticleCracking.SwellingAndCracking(p=params.n.side_reactions[1], s=params.n, V=params.e.Lₙ*A, g=g)
        @named cracking_p = ParticleCracking.SwellingOnly(p=params.n.side_reactions[1], s=params.p, V=params.e.Lₚ*A, g=g)
        @named lam_n = LAM.StressDriven(p=params.n.side_reactions[1], s=params.n, V=params.e.Lₙ*A, g=g)
        @named lam_p = LAM.StressDriven(p=params.n.side_reactions[1], s=params.p, V=params.e.Lₚ*A, g=g)
        push!(submodels, sei, plating, cracking_n, cracking_p, lam_n, lam_p)
    end

    @variables begin
        v(t), [guess=4.0]
        i(t)
        soc(t)

        U₀(t)
        ηᵣ(t)
        (ηₙ(t))[1:Nn]
        (ηₚ(t))[1:Np]
        Δϕₛ(t), [guess=0.0]
        
        (Δϕₙ(t))[1:Nn], [guess=fill(params.n.Uₖ(params.n.c₀/params.n.c₊), Nn)]
        (Δϕₚ(t))[1:Np], [guess=fill(params.p.Uₖ(params.p.c₀/params.p.c₊), Np)]
        Δϕₙ_x(t), [guess=params.n.Uₖ(params.n.c₀/params.n.c₊)]

        (ϕₙ(t))[1:g.el.Nx[1]]
        (ϕₚ(t))[1:g.el.Nx[3]]
        (jₙ0(t))[1:g.el.Nx[1]]
        (jₚ0(t))[1:g.el.Nx[3]]
        j̄ₙ0(t)
        j̄ₚ0(t)
    
        ϕ̄ₙ(t)
        ϕ̄ₚ(t)
        Rᵢ(t)

        Cₙ(t), [guess=params.e.Lₙ * A * params.n.c₊ * F / 3600]
        Cₚ(t), [guess=params.e.Lₚ * A * params.p.c₊ * F / 3600]
        C_cell(t), [guess=Q]
        
        # Heat generation sources
        Qᵢ(t)       
        Qₛ(t)       
        Qf(t)       
        Q_rev(t)    
        Q_total(t)  
        Q_sei(t)
        Q_cracking_sei(t), [guess=0.0]
        Q_plating(t)
        
        Q_loss_Ah(t)
        SoH(t)
        Q_Ah(t) = 0.0
        Qt_Ah(t) = 0.0
    end

    i_app = i/A
    x = g.el.x_centers
    L = sum(g.el.Ls)
    
    # Pure intercalation overpotential ignores parasitic side-reaction currents
    ηᵣn = [2*R*T.u/F*asinh((i_app/params.e.Lₙ/ne.aₖ)/(max(2*jₙ0[i], 1e-10))) for i in 1:Nn]
    ηᵣp = [2*R*T.u/F*asinh((-i_app/params.e.Lₚ/pe.aₖ)/(max(2*jₚ0[i], 1e-10))) for i in 1:Np]

    ϕₛ_n = [i_app*(x[i] - 2*params.e.Lₙ)*x[i]/2/params.n.σₖ/params.e.Lₙ for i in g.el.ixₙ]
    ϕₛ_p = [i_app*(x[i] + (x[i] - L)^2/(2*params.e.Lₚ))/params.p.σₖ for i in g.el.ixₚ]

    eqns = Equation[
        el.T.u ~ T.u,
        pe.T.u ~ T.u,
        ne.T.u ~ T.u,
        soc ~ ne.z,

        U₀ ~ pe.U₀ - ne.U₀,
        ηᵣ ~ sum(ηᵣp)/Np - sum(ηᵣn)/Nn,
        Δϕₛ ~ -i_app/3*(params.e.Lₚ/params.p.σₖ + params.e.Lₙ/params.n.σₖ),

        [jₙ0[i] ~ params.n.mₖ*NaNMath.sqrt(max(el.cₑ[g.el.ixₙ[i]]*ne.c_surf*(params.n.c₊-ne.c_surf), 1)) for i in 1:Nn]...,
        [jₚ0[i] ~ params.p.mₖ*NaNMath.sqrt(max(el.cₑ[g.el.ixₚ[i]]*pe.c_surf*(params.p.c₊-pe.c_surf), 1)) for i in 1:Np]...,
        j̄ₙ0 ~ sum(jₙ0)/Nn,
        j̄ₚ0 ~ sum(jₚ0)/Np,

        [ϕₙ[i] ~ n.v + ϕₛ_n[i] for i in 1:Nn]...,
        [ϕₚ[i] ~ p.v - ϕₛ_p[i] for i in 1:Np]...,
        [ηₙ[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
        [ηₚ[i] ~ ϕₚ[i] - el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        ϕ̄ₙ ~ sum(ϕₙ)/Nn,
        ϕ̄ₚ ~ sum(ϕₚ)/Np,
        
        v ~ U₀ + ηᵣ + el.ηₑ + el.Δϕₑ + Δϕₛ + (side_reactions ? sei.ϕf_x : 0.0) + (side_reactions ? cracking_n.ϕf_x : 0.0),
        
        Rᵢ ~ (U₀-v)/i, 
        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        [Δϕₙ[i] ~ ne.U₀ + ηᵣn[i] - (side_reactions ? sei.ϕf_x : 0.0) for i in 1:Nn]...,
        [Δϕₚ[i] ~ ϕₚ[i] - el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        Δϕₙ_x ~ sum(Δϕₙ)/Nn,

        el.i_app.u ~ i_app,
        el.ϕₛn.u ~ ϕ̄ₙ,
        el.Δϕₙ.u ~ Δϕₙ_x,

        ne.J.u ~ (i_app/params.e.Lₙ)/ne.aₖ - (side_reactions ? sei.j_sei_x : 0.0) - (side_reactions ? plating.j_stripping_x : 0.0),
        pe.J.u ~  -i_app/params.e.Lₚ/pe.aₖ, 

        [el.ϵ[i] ~ params.e.ϵₙ - (side_reactions ? (ne.aₖ*(sei.L_sei[i] - params.n.L_sei₀ + plating.L_plating[i] + plating.L_dead[i] + cracking_n.L_sei[i]*(cracking_n.r_surf - 1))) : 0.0) for i in g.el.ixₙ]...,
        [el.ϵ[i] ~ params.e.ϵₛ for i in g.el.ixₛ]...,
        [el.ϵ[i] ~ params.e.ϵₚ for i in g.el.ixₚ]...,

        Cₚ ~ pe.ϵₛ*params.e.Lₚ * A * params.p.c₊ * F / 3600,
        Cₙ ~ ne.ϵₛ*params.e.Lₙ * A * params.n.c₊ * F / 3600,
        C_cell ~ Q - Q_loss_Ah,
        Dt(Q_Ah) ~ i/3600,
        Dt(Qt_Ah) ~ abs(i)/3600,

        Qᵢ ~ -i_app * ηᵣ / L,
        Qₛ ~ -i_app * Δϕₛ / L,
        Qf ~ -i_app * (side_reactions ? sei.ϕf_x : 0.0) / L,
        Q_rev ~ (i_app / L) * T.u * (dUn_dT_f(ne.z) - dUp_dT_f(pe.z)),
        Q_sei ~ -(side_reactions ? sei.j_sei_x * sei.ϕf_x : 0.0) / L,
        Q_cracking_sei ~ -(side_reactions ? cracking_n.j_sei_x * cracking_n.ϕf_x : 0.0) / L,
        Q_plating ~ -(side_reactions ? plating.j_stripping_x * plating.ϕf_x : 0.0) / L,
        Q_total ~ el.Qₑ + Qᵢ + Qₛ + Qf + Q_rev + Q_sei + Q_cracking_sei + Q_plating,
        
        SoH ~ 100.0 * (1.0 - (Q_loss_Ah / params.Q₀))
    ]

    if side_reactions
        Volume_n = params.e.Lₙ * A
        c_sei_0 = params.n.side_reactions[1].Lf₀ / params.n.side_reactions[1].V̄ * params.n.aₖ
        
        moles_Li_sei = (sei.c_sei_x - c_sei_0) * Volume_n * params.n.side_reactions[1].z
        moles_Li_plating = (plating.c_dead_x + plating.c_plating_x) * Volume_n
        
        append!(eqns, [
            Q_loss_Ah ~ (moles_Li_sei + moles_Li_plating) * 96485.0 / 3600.0 + cracking_n.Q_sei + lam_n.Q_loss + lam_p.Q_loss,
            
            sei.J.u ~ i_app/params.e.Lₙ/ne.aₖ, 
            sei.T.u ~ T.u,
            sei.aₖ.u ~ ne.aₖ,
            sei.Δϕₛ.u ~ Δϕₙ,

            plating.J.u ~ i_app/params.e.Lₙ/ne.aₖ,
            plating.T.u ~ T.u,
            plating.aₖ.u ~ ne.aₖ,
            [plating.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
            plating.η_sei.u ~ sei.ϕf,
            [plating.cₑ.u[i] ~ params.e.cₜ for i in 1:Nn]...,
            [plating.L_sei.u[i] ~ sei.L_sei[i] for i in 1:Nn]...,

            cracking_n.J.u ~ i_app/params.e.Lₙ/ne.aₖ,
            cracking_n.T.u ~ T.u,
            cracking_n.aₖ.u ~ ne.aₖ,
            [cracking_n.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
            cracking_n.c_s_r.u ~ ne.c_r,
            cracking_n.c_s_surf.u ~ ne.c_surf,

            cracking_p.J.u ~ -i_app/params.e.Lₚ/pe.aₖ,
            cracking_p.T.u ~ T.u,
            cracking_p.aₖ.u ~ pe.aₖ,
            [cracking_p.Δϕₛ.u[i] ~ ϕₚ[i] + el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
            cracking_p.c_s_r.u ~ pe.c_r,
            cracking_p.c_s_surf.u ~ pe.c_surf,

            lam_n.σₜ.u ~ cracking_n.σₜ,
            lam_n.σᵣ.u ~ cracking_n.σᵣ,
            lam_n.c_r.u ~ ne.c_r,
            Dt(ne.ϵₛ) ~ lam_n.j_lam,
            
            lam_p.σₜ.u ~ cracking_p.σₜ,
            lam_p.σᵣ.u ~ cracking_p.σᵣ,
            lam_p.c_r.u ~ pe.c_r,
            Dt(pe.ϵₛ) ~ lam_p.j_lam
        ])
    else
        append!(eqns, [
            Q_loss_Ah ~ 0.0,
            Dt(ne.ϵₛ) ~ 0.0,
            Dt(pe.ϵₛ) ~ 0.0
        ])
    end

    events = Any[
        (v ~ params.Vmin) => (abort_vmin!, (;)),
        (v ~ params.Vmax) => (abort_vmax!, (;)),
        (pe.c_surf ~ params.p.c₊*0.9999) => (abort_pe_max!, (;)),
        (ne.c_surf ~ params.n.c₊*0.9999) => (abort_ne_max!, (;)),
        (pe.c_surf ~ params.p.c₊*0.0001) => (abort_pe_min!, (;)),
        (ne.c_surf ~ params.n.c₊*0.0001) => (abort_ne_min!, (;))
    ]

    append!(events, [(el.cₑ[i] ~ 1.0) => (abort_ce_min!, (;)) for i in 1:length(g.el.x_centers)])
    append!(events, [(el.ϵ[i] ~ 0.01) => (abort_porosity_min!, (;)) for i in g.el.ixₙ])

    return System(eqns, t; name=name, systems=submodels, continuous_events=events)
end