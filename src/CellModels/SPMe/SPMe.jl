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

# Import functions and equations
include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")
include("SEIGrowth.jl")
include("Heat.jl")

# Initialise potentials as symbolic functions (system does not look inside functions immediately)
@register_symbolic U₀_f(params::BatteryParameters, ne, pe)
@register_symbolic ηᵣ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, ne, pe,i_app, T)
@register_symbolic ηₑ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, T) 
@register_symbolic Δϕₑ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, ϵ::Symbolics.AbstractArray, i_app)
@register_symbolic Δϕₛ_f(params::BatteryParameters, g::NamedTuple, i_app)
@register_symbolic Δϕf_f(params::BatteryParameters, g::NamedTuple, i_app)

# Initialise heat sources as symbolic functions
@register_symbolic Qₑ_f(params, g, cₑ::Symbolics.AbstractArray, ϵ::Symbolics.AbstractArray, i_app, ηₑ) # Heat generated from electrolyte
@register_symbolic Qᵢ_f(params::BatteryParameters, g::NamedTuple, i_app, ηᵣ)                           # Irreversible heat
@register_symbolic Qₛ_f(params::BatteryParameters, g::NamedTuple, i_app, Δϕₛ)                           # Solid phase (electrodes) Ohmic Heat Generation
@register_symbolic Qf_f(params::BatteryParameters, g::NamedTuple, i_app, Δϕf)                          # Film Ohmic heat generation
@register_symbolic Q_rev_f(params::BatteryParameters, g::NamedTuple, i_app, T, c_s_n_surf, c_s_p_surf) # Reversible (entropic) heat generation

# Initiate end of experiment
function abort!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    return (;)
end

# SPMe implementation
function SPMe(; name="SPMe", params::BatteryParameters, Q=0, N=Dict(:Nₓ=>[10,10,10], :Nᵣ=>[10,10]), side_reactions=true)
    # Independent variables
    @parameters begin
        t # Time
    end

    # Build FVM geometry
    g = build_fvm_geometry(params, N)
    
    # Components
    @named p = Pin()                # Positive terminal
    @named n = Pin()                # Negative terminal
    @named T = RealInput(guess=298) # Ambient temperature
    # @named I = RealInput(guess=0.0) # Current
    # @named Q = RealOutput() # Generated heat
    @named pe = SolidParticle(p=params.p, g=g.pe) # Positive electrode
    @named ne = SolidParticle(p=params.n, g=g.ne) # Negative electrode
    @named el = Electrolyte(p=params.e, g=g.el)   # Electrolyte
    #@named sei = SEIGrowth(p=params.n.side_reactions[1], s=params.n, g=g) # Assuming first side reaction is SEI

    # Bundle all components
    submodels = [p, n, T, pe, ne, el]
    #submodels = [p, n, T, pe, ne, el, sei]
    
    # Side reactions
    # if !isnothing(params.sei)
    #     @named sei = SideReaction(name="SEI side reaction", p=params.sei)
    #     push!(submodels, sei)
    # end

    # if !isnothing(params.li_plating)
    #     @named plating = SideReaction(name="Lithium Plating", p=params.sei)
    #     push!(submodels, plating)
    # end

    # Time-dependent state variables
    @variables begin
        v(t), [guess=4.19]  # Terminal voltage
        i(t), [guess=0]  # Current
        soc(t)                    # State of charge

        #(ϕₙ(t))[1:g.el.Nx[1]]
        #(ϕₚ(t))[1:g.el.Nx[3]]

        U₀(t)    # Open-circuit potential
        ηᵣ(t)    # Intercalation reaction overpotential
        ηₑ(t)    # Electrolyte concentration overpotential
        Δϕₑ(t)   # Electrolyte Ohmic loss
        Δϕₛ(t)    # Separator Ohmic loss
        Δϕf(t)   # Film Ohmic loss
        # Rᵢ(t)   # Equivalent resistance

        Qₑ(t)       # Electrolyte Ohmic heat
        Qᵢ(t)       # Reaction heat
        Qₛ(t)        # Solid Ohmic heat
        Qf(t)       # Film Ohmic heat
        Q_rev(t)    # Reversible heat
        Q_total(t)  # Total heat
    end

    # Initialise capacity
    if Q==0
        Q = params.Q₀
    end
    
    # Applied current density (scaled to electrode area)
    i_app = i/Q*params.i₀
    
    # Concatenate porosity vectors
    ϵ = [
        [el.ϵₙ for _ in g.el.ixₙ]; # Negative electrode
        [el.ϵₛ for _ in g.el.ixₛ]; # Separator
        [el.ϵₚ for _ in g.el.ixₚ]; # Positive electrode
    ]

    # Time derivative
    D = Differential(t)
    
    # Equations
    eqns = [
        # Temperatures
        el.T.u ~ T.u, # Electrolyte
        pe.T.u ~ T.u, # Positive electrode
        ne.T.u ~ T.u, # Negative electrode

        # State of charge, function of stoichiometry of negative electrode (Fig. 4, p. 7 [2])
        soc ~ (ne.z - params.n.z_0) / (params.n.z_100 - params.n.z_0),

        # Potentials
        U₀ ~ U₀_f(params, ne.c_surf, pe.c_surf),                        # Open-circuit potential
        ηᵣ ~ ηᵣ_f(params, g, el.cₑ, ne.c_surf, pe.c_surf, i_app, T.u),  # Intercalation reaction overpotential
        ηₑ ~ ηₑ_f(params, g, el.cₑ, T.u),                               # Electrolyte concentration overpotential
        Δϕₑ ~ Δϕₑ_f(params, g, el.cₑ, ϵ, i_app),                        # Electrolyte Ohmic loss
        Δϕₛ ~ Δϕₛ_f(params, g, i_app),                                   # Separator Ohmic loss
        Δϕf ~ Δϕf_f(params, g, i_app),                                  # Film Ohmic loss
        v ~ U₀ + ηᵣ + ηₑ + Δϕₑ + Δϕₛ + Δϕf,                              # Terminal voltage (Eq. 12)
        # Rᵢ ~ (U₀-v)/i,                                                # Equivalent resistance

        # Heat sources
        Qₑ ~ Qₑ_f(params, g, el.cₑ, ϵ, i_app, ηₑ),
        Qᵢ ~ Qᵢ_f(params, g, i_app, ηᵣ),
        Qₛ ~ Qₛ_f(params, g, i_app, Δϕₛ),
        Qf ~ Qf_f(params, g, i_app, Δϕf),
        Q_rev ~ Q_rev_f(params, g, i_app, T.u, ne.c_surf, pe.c_surf),
        Q_total ~ Qₑ + Qᵢ + Qₛ + Qf + Q_rev,
        
        # i ~ I.u,
        v ~ p.v - n.v,  # Terminal voltage (Eq. 11)
        0 ~ p.i + n.i,  # Current in = current out
        i ~ p.i,        # Current

        # Electrolyte current density
        el.i_app.u ~ i_app,

        # Volumetric current density
        pe.J.u ~ -i_app/params.e.Lₚ, # Positive electrode
        ne.J.u ~  i_app/params.e.Lₙ, # Negative electrode
        #Jₛᵣ ~ parameters.p.mₖ * pp.c_avr^1.5 * (pp.Uₖ - parameters.p.Uₖ) # Side reaction, positive electrode
    ]

    # Events abort simulation when parameters go out of bounds
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

    # Return aggregrate of symbolic system equations
    return System(eqns, t; name=name, systems=submodels, continuous_events=events)
end