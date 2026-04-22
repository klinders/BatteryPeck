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
include("SEIGrowth.jl")

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

**Arguments**
- `name` (optional) defaults to SPMe
- `params ::BatteryParameters` Parameters from the given parameter set
- `Q ::Real` (optional) The capacity of the cell in Ah
- `N ::Dict(:Nₓ=>[::Int, ::Int, ::Int], :Nᵣ=>[::Int, ::Int])` (optional) number of mesh nodes in the particles and electrolyte
- `side_reactions ::Bool` (optional) enable side reactions

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

    # @named Q = RealOutput()
    @named pe = SolidParticle(p=params.p, g=g.pe)
    @named ne = SolidParticle(p=params.n, g=g.ne)
    @named el = Electrolyte(p=params.e, g=g.el)
    @named sei = SEIGrowth(p=params.n.side_reactions[1],s=params.n, g=g) # Assuming first side reaction is SEI
    @named plating = LithiumPlating(p=params.n.side_reactions[1],s=params.n, g=g)

    submodels = [p,n,T,pe,ne,el,sei,plating]

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
        # Δϕf(t)
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
    # Exchange current densities

    ## Reaction overpotentials ##

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    # jₙ0 = [params.n.mₖ*sqrt(el.cₑ[i]*cₛ[i]*(params.n.c₊-cₛ[i])) for i in g.el.ixₙ]
    # jₚ0 = [params.p.mₖ*sqrt(el.cₑ[i]*cₛ[i]*(params.p.c₊-cₛ[i])) for i in g.el.ixₚ]
    Nn = length(g.el.ixₙ)
    Np = length(g.el.ixₚ)

    # asin_n = [asinh(ne.J.u/params.n.aₖ/jₙ0[i]) for i in 1:Nn]
    # asin_p = [asinh(pe.J.u/params.p.aₖ/jₚ0[i]) for i in 1:Np]
    
    ηᵣn = 2*R*T.u/F*asinh(ne.J.u / (2*j̄ₙ0))
    ηᵣp = 2*R*T.u/F*asinh(pe.J.u / (2*j̄ₚ0))

    # ηᵣ_n = ηᵣ_x(params.n, cₙ, el.cₑ[g.el.ixₙ], ne.T.u, ne.J.u)
    # ηᵣ_p = ηᵣ_x(params.p, cₚ, el.cₑ[g.el.ixₚ], pe.T.u, pe.J.u)

    # X-average of the electrolyte potential
    ϕₛ_n = [i_app*(x[i] - 2*params.e.Lₙ)*x[i]/2/params.n.σₖ/params.e.Lₙ for i in g.el.ixₙ]
    ϕₛ_p = [i_app*(x[i] + (x[i] - L)^2/(2*params.e.Lₚ))/params.p.σₖ for i in g.el.ixₚ]

    eqns = [
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
        v ~ U₀ + ηᵣ + el.ηₑ + el.Δϕₑ + Δϕₛ + sei.ϕf_x,
        Rᵢ ~ (U₀-v)/i, 

        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        # Electrolyte current density
        el.i_app.u ~ i_app,
        # el.jₙ0.u ~ jₙ0,
        el.ϕₛn.u ~ ϕ̄ₙ,
        el.Δϕₙ.u ~ ne.U₀ + ηᵣn - sei.ϕf_x,

        # Volumetric current density
        ne.J.u ~ (i_app/params.e.Lₙ)/aₙ - sei.j_sei_x - plating.j_stripping_x, # Current density in the negative electrode
        pe.J.u ~  -i_app/params.e.Lₚ/aₚ, # Current density in the positive electrode
        
        aₙ ~ 3*(1-el.ϵ̄ₙ)/params.n.Rₖ,
        aₚ ~ params.p.aₖ,#3*(1-el.ϵ̄ₚ)/params.p.Rₖ,

        # # Ne sei reaction
        sei.J.u ~ ne.J.u, # Current density for SEI side reaction
        sei.T.u ~ T.u,
        sei.aₖ.u ~ aₙ,
        [sei.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[i] for i in 1:Nn]...,

        # Li plating
        plating.J.u ~ne.J.u,
        plating.T.u ~ T.u,
        plating.aₖ.u ~ aₙ,
        [plating.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[i] for i in 1:Nn]...,
        plating.η_sei.u ~ sei.ϕf,
        [plating.cₑ.u[i] ~ params.e.cₜ for i in 1:Nn]...,
        
        # Porosity (assumed constant)
        [el.ϵ[i] ~ params.e.ϵₙ - aₙ*(sei.L_sei[i] - params.n.L_sei₀) for i in g.el.ixₙ]...,
        [el.ϵ[i] ~ params.e.ϵₛ for i in g.el.ixₛ]...,
        [el.ϵ[i] ~ params.e.ϵₚ for i in g.el.ixₚ]...,

        # Heat sources
        Qᵢ ~ -i_app * ηᵣ / L,
        Qₛ ~ -i_app * Δϕₛ / L,
        Qf ~ -i_app * sei.ϕf_x / L,
        Q_rev ~ (i_app / L) * T.u * (dUn_dT_f(ne.z) - dUp_dT_f(pe.z)),
        Q_sei ~ -sei.j_sei_x * sei.ϕf_x / L,
        Q_plating ~ -plating.j_stripping_x * plating.ϕf_x / L,
        
        # Total generated heat 
        Q_total ~ el.Qₑ + Qᵢ + Qₛ + Qf + Q_rev + Q_sei + Q_plating
    ]

    # Terminate simulation when limits reached
    events = [
        [
            v ~ params.Vmin,
            v ~ params.Vmax,
            pe.c_surf ~ params.p.c₊*0.999,   # Upper bound
            ne.c_surf ~ params.n.c₊*0.999,
            pe.c_surf ~ params.p.c₊*0.001,   # Lower bound
            ne.c_surf ~ params.n.c₊*0.001
        ]=>(abort!,(;))
    ]

    return System(eqns, t; name=name,systems=submodels, continuous_events=events)
end