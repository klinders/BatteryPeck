
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical


include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")
include("SEI.jl")
include("LithiumPlating.jl")
include("ParticleCracking.jl")

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
    @named T = RealInput(guess=298.15)

    # @named Q = RealOutput()
    @named pe = SolidParticle(p=params.p, g=g.pe)
    @named ne = SolidParticle(p=params.n, g=g.ne)
    @named el = Electrolyte(p=params.e, g=g.el)
    @named sei = SEI.SolventDiffusionLimitedSEI(p=params.n.side_reactions[1],s=params.n, g=g) # Assuming first side reaction is SEI
    @named plating = LithiumPlating.PartiallyReversiblePlating(p=params.n.side_reactions[1],s=params.n, g=g) 
    @named cracking = ParticleCracking.SwellingAndCracking(p=params.n.side_reactions[1], s=params.n, g=g)

    submodels = [p,n,T,pe,ne,el,sei,plating,cracking]

    @variables begin
        # Terminal voltage and current
        v(t)
        i(t)
        soc(t)

        U₀(t)
        ηᵣ(t)
        ηᵣ̅n(t)
        ηᵣ̅p(t)
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
        Q_loss(t)
        Q_Ah(t) = 0
        Qt_Ah(t) = 0

        #temp
        j_tot_ne(t)
        aj_tot_ne(t)

        Rᵢ(t)
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
    
    # Overpotentials from the inverse Butler-Volmer equation
    ηᵣn = [2*R*T.u/F*asinh((i_app/params.e.Lₙ/aₙ)/(2*jₙ0[i])) for i in 1:Nn]
    ηᵣp = [2*R*T.u/F*asinh((-i_app/params.e.Lₚ/aₚ)/(2*jₚ0[i])) for i in 1:Np]

    # ηᵣ_n = ηᵣ_x(params.n, cₙ, el.cₑ[g.el.ixₙ], ne.T.u, ne.J.u)
    # ηᵣ_p = ηᵣ_x(params.p, cₚ, el.cₑ[g.el.ixₚ], pe.T.u, pe.J.u)

    # X-average of the electrolyte potential
    ϕₛ_n = [i_app*(x[i] - 2*params.e.Lₙ)*x[i]/2/params.n.σₖ/params.e.Lₙ for i in g.el.ixₙ]
    ϕₛ_p = [i_app*(x[i] + (x[i] - L)^2/(2*params.e.Lₚ))/params.p.σₖ for i in g.el.ixₚ]

    eqns = [
        # Temps
        el.T.u ~ T.u,
        pe.T.u ~ T.u,
        ne.T.u ~ T.u,
        soc ~ ne.z,

        # temp
        aj_tot_ne ~ i_app/params.e.Lₙ,
        j_tot_ne ~ aj_tot_ne/aₙ,

        ## Potentials ##
        U₀ ~ pe.U₀ - ne.U₀,
        ηᵣ̅n ~ sum(ηᵣn)/Nn,
        ηᵣ̅p ~ sum(ηᵣp)/Np,
        ηᵣ ~ ηᵣ̅p - ηᵣ̅n,
        Δϕₛ ~ -i_app/3*(params.e.Lₚ/params.p.σₖ + params.e.Lₙ/params.n.σₖ),

        # Exchange current densities
        [jₙ0[i] ~ params.n.j0(el.cₑ[g.el.ixₙ[i]], ne.c_surf, params.n.c₊, T.u) for i in 1:Nn]...,
        [jₚ0[i] ~ params.p.j0(el.cₑ[g.el.ixₚ[i]], pe.c_surf, params.p.c₊, T.u) for i in 1:Np]...,
        j̄ₙ0 ~ sum(jₙ0)/Nn,
        j̄ₚ0 ~ sum(jₚ0)/Np,

        [ϕₙ[i] ~ n.v + ϕₛ_n[i] for i in 1:Nn]...,
        [ϕₚ[i] ~ p.v - ϕₛ_p[i] for i in 1:Np]...,
        [ηₙ[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
        [ηₚ[i] ~ ϕₚ[i] - el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        ϕ̄ₙ ~ sum(ϕₙ)/Nn,
        ϕ̄ₚ ~ sum(ϕₚ)/Np,
        v ~ U₀ + ηᵣ + el.ηₑ + el.Δϕₑ + Δϕₛ + sei.ϕf_x + cracking.ϕf_x,
        Rᵢ ~ (U₀-v)/i, 

        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        # Electrolyte current density
        el.i_app.u ~ i_app,
        # el.jₙ0.u ~ jₙ0,
        el.ϕₛn.u ~ ϕ̄ₙ,
        el.Δϕₙ.u ~ ne.U₀ + ηᵣ̅n - sei.ϕf_x,

        ne.J.u ~  (i_app/params.e.Lₙ - sei.j_sei_x)/aₙ, # Current density in the negative electrode
        pe.J.u ~  -i_app/params.e.Lₚ/aₚ, # Current density in the positive electrode
        
        aₙ ~ 3*(1-params.e.ϵₙ)/params.n.Rₖ,
        aₚ ~ params.p.aₖ,

        # # Ne sei reaction
        sei.J.u ~ ne.J.u, # Current density for SEI side reaction
        sei.T.u ~ T.u,
        sei.aₖ.u ~ aₙ,
        [sei.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[i] for i in 1:Nn]...,

        # # Li plating
        plating.J.u ~ne.J.u,
        plating.T.u ~ T.u,
        plating.aₖ.u ~ aₙ,
        [plating.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[i] for i in 1:Nn]...,
        plating.η_sei.u ~ sei.ϕf,
        [plating.cₑ.u[i] ~ el.cₑ[i] for i in 1:Nn]...,
        [plating.L_sei.u[i] ~ sei.L_sei[i] for i in 1:Nn]...,

        ## Cracking
        cracking.J.u ~ ne.J.u,
        cracking.T.u ~ T.u,
        cracking.aₖ.u ~ aₙ,
        [cracking.Δϕₛ.u[i] ~ ϕₙ[i] - el.ϕₑ[i] for i in 1:Nn]...,
        cracking.c_s_r.u ~ ne.c_r,
        cracking.c_s_surf.u ~ ne.c_surf,

        
        # Porosity (assumed constant)
        [el.ϵ[i] ~ params.e.ϵₙ - aₙ*(
            sei.L_sei[i] - params.n.L_sei₀ 
            + plating.L_plating[i] 
            + plating.L_dead[i]
             + cracking.L_sei[i]*(cracking.r_surf - 1)
            ) for i in g.el.ixₙ]...,
        [el.ϵ[i] ~ params.e.ϵₛ for i in g.el.ixₛ]...,
        [el.ϵ[i] ~ params.e.ϵₚ for i in g.el.ixₚ]...,

        Q_loss ~ sei.Q_loss + plating.Q_loss,

        Dt(Q_Ah) ~ i/3600,
        Dt(Qt_Ah) ~ abs(i)/3600

    ]

    # Event working
    events = [
        [
            v ~ params.Vmin,
            v ~ params.Vmax,
            pe.c_surf ~ params.p.c₊*0.99,
            ne.c_surf ~ params.n.c₊*0.99,
        ]=>(abort!,(;))
    ]

    return System(eqns, t; name=name,systems=submodels, continuous_events=events)
end