
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("SolidParticle.jl")
include("Electrolyte.jl")
include("../ParameterSets/Base.jl")


function SPMe(; name="SPMe", params::BatteryParameters, Q=0, N=Dict(:Nₓ=>[10,10,10], :Nᵣ=>[10,10]), side_reactions=true)
    @parameters begin
        t # Time variable
    end
    
    # Electrical ports
    @named p = Pin()
    @named n = Pin()
    @named pe = SolidParticle(p=params.p, Nᵣ=N[:Nᵣ][1])
    @named ne = SolidParticle(p=params.n, Nᵣ=N[:Nᵣ][2])
    @named el = Electrolyte(p=params.e, Nₓ=N[:Nₓ])

    @variables begin
        # Terminal voltage and current
        v(t)
        i(t)
        #Jsr(t) # Side reaction current density
    end


    if Q==0
        Q = params.Q₀
    end

    # Spatial variable in x direction
    # Δx = [
    #     params.e.Lₙ/sum(N[:Nₓ])*ones(N[:Nₓ][1]),
    #     params.e.Lₛ/sum(N[:Nₓ])*ones(N[:Nₓ][2]),
    #     params.e.Lₚ/sum(N[:Nₓ])*ones(N[:Nₓ][3])
    # ]
    # x_nodes = [(i-0.5)*]
    # # Negative electrode potential
    # term1 =  (i_app*(2*params.e.Lₙ-x)*x)

    # ϕₙ = ne.U₀ 
    
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