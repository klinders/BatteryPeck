
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")

@register_symbolic U₀_f(params::BatteryParameters, ne, pe)
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
    # @named I = RealInput(guess=0.0)

    # @named Q = RealOutput()
    @named pe = SolidParticle(p=params.p, g=g.pe)
    @named ne = SolidParticle(p=params.n, g=g.ne)
    @named el = Electrolyte(p=params.e, g=g.el)

    @variables begin
        # Terminal voltage and current
        v(t)
        i(t)
        #Jsr(t) # Side reaction current density
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
        [el.ϵₙ for _ in g.el.ixₙ]...
        [el.ϵₛ for _ in g.el.ixₛ]...
        [el.ϵₚ for _ in g.el.ixₚ]...
    ]

    D = Differential(t)
    
    eqns = [
        # Temps
        el.T.u ~ T.u,
        pe.T.u ~ T.u,
        ne.T.u ~ T.u,

        # # Potentials
        U₀ ~ U₀_f(params, ne.c_surf, pe.c_surf),
        ηᵣ ~ ηᵣ_f(params, g, el.cₑ, ne.c_surf, pe.c_surf, i_app, T.u),
        ηₑ ~ ηₑ_f(params, g, el.cₑ, T.u),
        Δϕₑ ~ Δϕₑ_f(params, g, el.cₑ, ϵ, i_app),
        Δϕₛ ~ Δϕₛ_f(params, g, i_app),
        Δϕf ~ Δϕf_f(params, g, i_app), 
        v ~ U₀ + ηᵣ + ηₑ + Δϕₑ + Δϕₛ + Δϕf,
        Rᵢ ~ (U₀-v)/i,
        
        # i ~ I.u,
        v ~ p.v + n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        # Electrolyte current density
        el.i_app.u ~ i_app,

        pe.J.u ~ -i_app/params.e.Lₚ, # Current density in the positive electrode
        ne.J.u ~  i_app/params.e.Lₙ, # Current density in the negative electrode

        #Jₛᵣ ~ parameters.p.mₖ * pp.c_avr^1.5 * (pp.Uₖ - parameters.p.Uₖ) # Side reaction current density in the positive electrode
    ]

    # Event work but very slow
    events = [
        [
            v ~ params.Vmin,
            v ~ params.Vmax,
            pe.c_surf ~ params.p.c₊*0.99,
            ne.c_surf ~ params.n.c₊*0.99,
        ]=>(abort!,(;))
    ]

    return System(eqns, t; name=name,systems=[p,n, pe, ne, el, T], continuous_events=events)
end