# =====================================================================================================================
# Potentials.jl
#
# Potential equations for SPMe
# Source: https://doi.org/10.1016/j.apm.2022.12.009
# =====================================================================================================================

# Intercalation reaction overpotential
function ηᵣ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, ne, pe, i_app, T)

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday constant
    T = 298 # Ambient temperature
    
    # Retrieve FVM geometry
    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    
    # Concentrations
    cₚ = pe
    cₙ = ne
    cₑ = el
    
    # Volumetric current densities
    jₚ = params.p.mₖ.*sqrt.(max(cₑ[g.el.ixₚ].*cₚ.*(params.p.c₊-cₚ), 1e-10))
    jₙ = params.n.mₖ.*sqrt.(max(cₑ[g.el.ixₙ].*cₙ.*(params.n.c₊-cₙ), 1e-10))
    
    # Hyperbolic arcsine terms for inverse Butler-Volmer equation
    asin_p = asinh.(i_app./(max(params.p.aₖ.*Lₚ.*jₚ, 1e-10)))
    asin_n = asinh.(i_app./(max(params.n.aₖ.*Lₙ.*jₙ, 1e-10)))
    
    # Integral terms
    sum_p = ∫(asin_p, x[g.el.ixₚ])
    sum_n = ∫(asin_n, x[g.el.ixₙ])
    
    ηᵣ = -2*R*T/F*(sum_p/Lₚ + sum_n/Lₙ)

    return ηᵣ
end

# Electrolyte concentration overpotential
function ηₑ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, T)

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday constant
    # T = 298 # Ambient temperature

    # Retrieve FVM geometry
    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    Nₜ = g.el.Nₜ

    # Concentration in electrolyte
    cₑ = el

    # Electrolyte potential drop
    df_fac = ones(length(cₑ))

    logc = log.(max(cₑ, 1e-10))

    # Derivative using central difference 
    dlogc_dx = [
        (logc[2]-logc[1])/(x[2]-x[1]),
        [(logc[i+1]-logc[i-1])/(x[i+1]-x[i-1]) for i in 2:Nₜ-1]...,
        (logc[end]-logc[end-1])/(x[end]-x[end-1])
    ]

    # Integrand
    f1 = [(1 - params.e.t₊(cₑ[i]))*df_fac[i]*dlogc_dx[i] for i in 1:Nₜ]

    # Inner integral domain for current density
    int1 = cumsum([
        0,
        [(f1[i] + f1[i-1])*(x[i]-x[i-1])/2 for i in 2:Nₜ]...,
    ])

    # Outer integral domain on electrodes
    int1_p = ∫(int1, x[g.el.ixₚ])
    int1_n = ∫(int1, x[g.el.ixₙ])

    ηₑ = 2*R*T/F*(int1_p/Lₚ - int1_n/Lₙ)
end

# Electrolyte Ohmic loss
function Δϕₑ_f(params::BatteryParameters, g::NamedTuple, cₑ::Symbolics.AbstractArray, ϵ::Symbolics.AbstractArray, i_app)

    # Retrieve FVM geometry
    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    L = sum(g.el.Ls)
    Nₜ = g.el.Nₜ

    # Bruggeman coefficients per region
    b = [
        [params.e.bₙ for _ in g.el.ixₙ]
        [params.e.bₛ for _ in g.el.ixₛ]
        [params.e.bₚ for _ in g.el.ixₚ]
    ]

    # Transport efficiency multiplied by electrolyte conductivity
    B = [params.e.σₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:Nₜ]

    # Current in electrolyte
    function iₑ(x)
        if x <= Lₙ
            return i_app*x/Lₙ
        elseif x <= Lₙ + Lₛ
            return i_app
        else
            return i_app*(L-x)/Lₚ
        end
    end

    # Integrand
    f2 = [iₑ(x[i])/B[i] for i in 1:Nₜ]

    # Inner integral domain for current density
    int2 = cumsum([
        (f2[1] + iₑ(0)/B[1])*(x[1])/2,
        [(f2[i] + f2[i-1])*(x[i]-x[i-1])/2 for i in 2:Nₜ]...,
    ])

    # Outer integral domain on electrodes
    int2_p = ∫(int2, x[g.el.ixₚ])
    int2_n = ∫(int2, x[g.el.ixₙ])

    Δϕₑ = -int2_p/Lₚ + int2_n/Lₙ
end

# Separator Ohmic loss
function Δϕₛ_f(params::BatteryParameters, g::NamedTuple, i_app)

    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]

    Δϕₛ = -i_app/3*(Lₚ/params.p.σₖ + Lₙ/params.n.σₖ)
end

# Film Ohmic loss
function Δϕf_f(params::BatteryParameters, g::NamedTuple, i_app)

    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]

    Lf = params.n.L_sei₀
    Δϕf = -i_app*(Lf/Lₙ/params.n.aₖ/params.n.σₖ)
end

# Average reaction overpotential
function ηᵣ_x(params::SolidParticleParameters, cₖ, cₑ, j, T)

    @assert length(cₖ) == length(cₑ)

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday constant
    N = length(cₖ)

    jₓ0 = [params.mₖ*sqrt(max(cₑ[i]*cₖ[i]*(params.c₊-cₖ[i]), 1e-10)) for i in 1:N]
    
    asin = [asinh(j/(max(params.aₖ*jₓ0[i], 1e-10))) for i in 1:N]
    
    sinh_x = 2*R*T/F*sum(asin)/N

    return sinh_x
end

# Average electrolyte concentration overpotential
function ηₑ_x(params::ElectrolyteParameters, cₑ, x, T)
    
    @assert length(cₑ) == length(x)
    N = length(cₑ)
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday constant
    
    dcₑ_dx = [
        (cₑ[2]-cₑ[1])/(x[2]-x[1]),
        [(cₑ[i+1]-cₑ[i-1])/(x[i+1]-x[i-1]) for i in 2:N]...,
        (cₑ[end]-cₑ[end-1])/(x[end]-x[end-1])
    ]

    f = [(1 - params.t₊(cₑ[i]))*dcₑ_dx[i]/cₑ[i] for i in 1:N]

    int = cumsum([
        f[1]*(x[1])/2,
        [(f[i] + f[i-1])*(x[i]-x[i-1])/2 for i in 2:N]...,
    ])

    return 2*R*T/F*sum(int)/N
end