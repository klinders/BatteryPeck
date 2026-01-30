
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")
include("SEIGrowth.jl")

@register_symbolic ηᵣ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, ne, pe,i_app, T)
@register_symbolic ηₑ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, T) 
@register_symbolic Δϕₑ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, ϵ::Symbolics.AbstractArray, i_app)
@register_symbolic Δϕₛ_f(params::BatteryParameters, g::NamedTuple, i_app)
@register_symbolic Δϕf_f(params::BatteryParameters, g::NamedTuple, i_app)

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
    @named sei = SEIGrowth(p=params.n.side_reactions[1], g=g) # Assuming first side reaction is SEI

    submodels = [p,n,T,pe,ne,el,sei]

    @variables begin
        # Terminal voltage and current
        v(t)
        i(t)
        soc(t)

        U₀(t)
        ηᵣ(t)
        ηₑ(t)
        Δϕₑ(t)
        Δϕₛ(t)
        Δϕf(t)
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

    D = Differential(t)

    # X-average
    x = g.el.x_centers
    σ_f = 5e-6
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    Lsei_0 = 5e-09

    ϕₑ_n = sum([
            el.ϕₑ[1]*g.el.ixₙ[1]/2,
            [(el.ϕₑ[i] + el.ϕₑ[i-1])*(g.el.ixₙ[i] - g.el.ixₙ[i-1])/2 for i in 2:length(g.el.ixₙ)]...,
    ])/params.e.Lₙ

    j = [params.p.mₖ*sqrt(el.cₑ[i]*pe.c_surf*(params.p.c₊-pe.c_surf)) for i in 1:g.el.Nₜ]
    
    asin_p = [asinh(i_app/params.p.aₖ/params.e.Lₚ/j[i]) for i in g.el.ixₚ]
    asin_n = [asinh(i_app/params.n.aₖ/params.e.Lₙ/j[i]) for i in g.el.ixₙ]
    
    sinh_x_p = 2*R*T.u/F*∫(asin_p, x[g.el.ixₚ])/params.e.Lₚ
    
    sinh_x_n = 2*R*T.u/F*sum([
            asin_n[1]*g.el.ixₙ[1]/2,
            [(asin_n[i] + asin_n[i-1])*(g.el.ixₙ[i] - g.el.ixₙ[i-1])/2 for i in 2:length(g.el.ixₙ)]...,
    ])/params.e.Lₙ

    ϕₛ = [ne.U₀ - i_app*(2*params.e.Lₙ − x[i])*x[i]/(2*params.e.Lₙ*params.n.σₖ) + i_app*params.e.Lₙ/(3*params.n.σₖ) + el.ϕₑ[i] + sinh_x_n + i_app*Lsei_0/params.e.Lₙ/params.n.aₖ/σ_f for i in g.el.ixₙ]
    ϕₛ_n = ∫(ϕₛ, g.el.ixₙ)/params.e.Lₙ
    
    @show el.ϕₑ[2]
    @show ϕₛ[2]

    eqns = [
        # Temps
        el.T.u ~ T.u,
        pe.T.u ~ T.u,
        ne.T.u ~ T.u,
        soc ~ ne.z,

        # # Potentials
        U₀ ~ pe.U₀ - ne.U₀,
        ηᵣ ~ ηᵣ_f(params, g, el.cₑ, ne.c_surf, pe.c_surf, i_app, T.u),
        ηₑ ~ ηₑ_f(params, g, el.cₑ, T.u),
        Δϕₑ ~ Δϕₑ_f(params, g, el.cₑ, ϵ, i_app),
        Δϕₛ ~ Δϕₛ_f(params, g, i_app),
        Δϕf ~ Δϕf_f(params, g, i_app), 
        v ~ U₀, #+ ηᵣ + ηₑ + Δϕₑ + Δϕₛ + Δϕf,
        Rᵢ ~ (U₀-v)/i, 
        
        # i ~ I.u,
        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        # Electrolyte current density
        el.i_app.u ~ i_app,

        pe.J.u ~ -i_app/params.e.Lₚ, # Current density in the positive electrode
        ne.J.u ~  i_app/params.e.Lₙ, # Current density in the negative electrode

        # # Ne sei reaction
        sei.J.u ~ -ne.J.u, # Current density for SEI side reaction
        sei.T.u ~ T.u,
        [sei.ϕₑ.u[i] ~ el.ϕₑ[i] for i in g.el.ixₙ]...,
        [sei.ϕₛ.u[i] ~ ϕₛ[i] for i in g.el.ixₙ]...,

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