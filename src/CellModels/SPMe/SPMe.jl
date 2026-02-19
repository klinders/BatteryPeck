
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")
include("SEIGrowth.jl")

function abort!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    return (;)
end

function SPMe(; name="SPMe", params::BatteryParameters, Q=0, N=Dict(:Nₓ=>[10,10,10], :Nᵣ=>[10,10]), side_reactions=true)
    @parameters begin
        t # Time variable
    end

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

    submodels = [p,n,T,pe,ne,el,sei]

    @variables begin
        # Terminal voltage and current
        v(t)
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

        Rᵢ(t)
    end

    if Q==0
        Q = params.Q₀
    end
    
    # Scale the current density to the electrode area
    i_app = i/Q*params.i₀
    
    # Porosity per region
    ϵ = [
        [el.ϵₙ for _ in g.el.ixₙ];
        [el.ϵₛ for _ in g.el.ixₛ];
        [el.ϵₚ for _ in g.el.ixₚ];
    ]

    # X-average
    x = g.el.x_centers
    L = sum(g.el.Ls)
    # Exchange current densities
    

    # Repeat the electrode surface potential for consistency
    cₛ = [
        [ne.c_surf for _ in g.el.ixₙ]...,
        [0 for _ in g.el.ixₛ]...,
        [pe.c_surf for _ in g.el.ixₚ]...
    ]

    ## Reaction overpotentials ##

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    # jₙ0 = [params.n.mₖ*sqrt(el.cₑ[i]*cₛ[i]*(params.n.c₊-cₛ[i])) for i in g.el.ixₙ]
    # jₚ0 = [params.p.mₖ*sqrt(el.cₑ[i]*cₛ[i]*(params.p.c₊-cₛ[i])) for i in g.el.ixₚ]
    Nn = length(g.el.ixₙ)
    Np = length(g.el.ixₚ)

    # asin_n = [asinh(ne.J.u/params.n.aₖ/jₙ0[i]) for i in 1:Nn]
    # asin_p = [asinh(pe.J.u/params.p.aₖ/jₚ0[i]) for i in 1:Np]
    
    ηᵣn = 2*R*T.u/F*asinh(ne.J.u/params.n.aₖ/2/j̄ₙ0)
    ηᵣp = 2*R*T.u/F*asinh(pe.J.u/params.p.aₖ/2/j̄ₚ0)

    # ηᵣ_n = ηᵣ_x(params.n, cₙ, el.cₑ[g.el.ixₙ], ne.T.u, ne.J.u)
    # ηᵣ_p = ηᵣ_x(params.p, cₚ, el.cₑ[g.el.ixₚ], pe.T.u, pe.J.u)

    # X-average of the electrolyte potential
    ϕₛ_n = [i_app*(x[i] - 2*params.e.Lₙ)*x[i]/2/params.n.σₖ/params.e.Lₙ for i in g.el.ixₙ]
    ϕₛ_p = [i_app*(x[i] - L)*(L - 2*params.e.Lₚ - x[i])/2/params.p.σₖ/params.e.Lₚ for i in g.el.ixₚ]

    eqns = [
        # Temps
        el.T.u ~ T.u,
        pe.T.u ~ T.u,
        ne.T.u ~ T.u,
        soc ~ ne.z,

        ## Potentials ##
        U₀ ~ pe.U₀ - ne.U₀,
        ηᵣ ~ ηᵣp - ηᵣn,
        Δϕₛ ~ -i_app/3*(params.e.Lₚ/params.p.σₖ + params.e.Lₙ/params.n.σₖ),

        # Exchange current densities
        [jₙ0[i] ~ params.n.mₖ*sqrt(el.cₑ[g.el.ixₙ[i]]*ne.c_surf*(params.n.c₊-ne.c_surf)) for i in 1:Nn]...,
        [jₚ0[i] ~ params.p.mₖ*sqrt(el.cₑ[g.el.ixₚ[i]]*pe.c_surf*(params.p.c₊-pe.c_surf)) for i in 1:Np]...,
        j̄ₙ0 ~ sum(jₙ0)/Nn,
        j̄ₚ0 ~ sum(jₚ0)/Np,

        [ϕₙ[i] ~ n.v + ϕₛ_n[i] for i in 1:Nn]...,
        [ϕₚ[i] ~ p.v + ϕₛ_p[i] for i in 1:Np]...,
        [ηₙ[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
        [ηₚ[i] ~ ϕₚ[i] - el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        ϕ̄ₙ ~ sum(ϕₙ)/Nn,
        ϕ̄ₚ ~ sum(ϕₚ)/Np,
        v ~ U₀ + ηᵣ + el.ηₑ + el.Δϕₑ + Δϕₛ + sei.ϕf_av,
        Rᵢ ~ (U₀-v)/i, 

        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        # Electrolyte current density
        el.i_app.u ~ i_app,
        # el.jₙ0.u ~ jₙ0,
        el.ϕₛn.u ~ ϕ̄ₙ,
        el.Δϕₙ.u ~ ne.U₀ + ηᵣn + sei.ϕf_av ,#(ϕ̄ₙ - el.ϕ̄ₑn - ne.U₀)*log(ne.c_surf/params.n.c₊), # i_app/sqrt(j̄ₙ0^2*params.e.Lₙ^2*params.n.aₖ^2 + i_app^2)*R*T.u/F

        pe.J.u ~  -i_app/params.e.Lₚ, # Current density in the positive electrode
        ne.J.u ~  i_app/params.e.Lₙ, # Current density in the negative electrode

        # # Ne sei reaction
        sei.J.u ~ ne.J.u, # Current density for SEI side reaction
        sei.T.u ~ T.u,
        [sei.ϕₑ.u[i] ~ el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
        [sei.ϕₛ.u[i] ~ ϕₙ[i] for i in 1:Nn]...,

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