
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SolidParticle.jl")
include("Electrolyte.jl")

function SPMe(; name="SPMe", params::BatteryParameters, Q=0, N=Dict(:Nₓ=>[10,10,10], :Nᵣ=>[10,10]), side_reactions=true)
    @parameters begin
        t # Time variable
    end

    g = build_fvm_geometry(params, N)
    
    # Electrical ports
    @named p = Pin()
    @named n = Pin()
    @named pe = SolidParticle(p=params.p, g=g.pe)
    @named ne = SolidParticle(p=params.n, g=g.ne)
    @named el = Electrolyte(p=params.e, g=g.el)

    @variables begin
        # Terminal voltage and current
        v(t)
        i(t)
        #Jsr(t) # Side reaction current density
    end


    if Q==0
        Q = params.Q₀
    end

    # Potentials
    ϕₙ_array = 
    
    # Scale the current density to the electrode area
    i_app = i/Q*params.i₀
    
    eqns = [
        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,

        # Electrolyte current density
        el.i_app.u ~ i_app,

        pe.J.u ~ -i_app/params.e.Lₚ, # Current density in the positive electrode
        ne.J.u ~  i_app/params.e.Lₙ, # Current density in the negative electrode

        #Jₛᵣ ~ parameters.p.mₖ * pp.c_avr^1.5 * (pp.Uₖ - parameters.p.Uₖ) # Side reaction current density in the positive electrode
    ]

    return System(eqns, t; name=name,systems=[p,n, pe, ne, el])
end


function ϕₙ(params::BatteryParameters, g, el, ne, i_app)
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    T = 298 # Temperature

    Nx = g.Nx
    ixₙ = g.ixₙ
    Δx = g.Δx
    Δxₗ = g.Δxₗ
    Δxᵣ = g.Δxᵣ

    # Bruggeman coefficients per region
    b = [
        [params.e.bₙ for _ in g.ixₙ] 
        [params.e.bₛ for _ in g.ixₛ]
        [params.e.bₚ for _ in g.ixₚ] 
    ]

    σ_eff = [params.n.σ_eff for _ in 1:sum(Nx)] # Effective conductivity per cell

    ϕₙ_array = [
        el.ϕₑ[i] + ne.U₀ + R*T/F*log(ne.c_surf/params.n.c₊) for i in ixₙ
    ]

    return ϕₙ_array
end