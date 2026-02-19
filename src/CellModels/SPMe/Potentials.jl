# =====================================================================================================================
# Potentials.jl
#
# Potential equations for SPMe (Eq. 13)
# Source: https://doi.org/10.1016/j.apm.2022.12.009
# =====================================================================================================================

# Import packages
import NaNMath # Return NaN for log or sqrt of -1 (potentially caused by electrolyte depletion at high C-rates)

# Open-circuit potential
function U₀_f(params::BatteryParameters, ne, pe)
    # Concentrations (used to calculate stoichiometry)
    cₚ = pe # Positive electrode
    cₙ = ne # Negative electrode

    # [U_p(z) - U_n(z)]
    U₀ = params.p.Uₖ(cₚ/params.p.c₊) - params.n.Uₖ(cₙ/params.n.c₊)
end

# Intercalation reaction overpotential
function ηᵣ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, ne, pe, i_app, T)
    # Constants
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    T = 298   # Ambient temperature
    
    # Retrieve FVM geometry
    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    
    # Concentrations
    cₚ = pe # Positive electrode
    cₙ = ne # Negative electrode
    cₑ = el # Electrolyte
    
    # Volumetric current densities
    jₚ = params.p.mₖ.*NaNMath.sqrt.(cₑ[g.el.ixₚ].*cₚ.*(params.p.c₊-cₚ)) # Positive electrode
    jₙ = params.n.mₖ.*NaNMath.sqrt.(cₑ[g.el.ixₙ].*cₙ.*(params.n.c₊-cₙ)) # Negative electrode
    
    # Hyperbolic arcsine terms (inverse Butler-Volmer current density equation)
    asin_p = asinh.(i_app./params.p.aₖ./Lₚ./jₚ)
    asin_n = asinh.(i_app./params.n.aₖ./Lₙ./jₙ)
    
    # Integral terms
    sum_p = ∫(asin_p, x[g.el.ixₚ])
    sum_n = ∫(asin_n, x[g.el.ixₙ])
    
    # Intercalation reaction overpotential
    ηᵣ = -2*R*T/F*(sum_p/Lₚ + sum_n/Lₙ)
end

# Electrolyte concentration overpotential
function ηₑ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, T)
    # Constants
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    # T = 298 # Ambient temperature

    # Retrieve FVM geometry
    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    Nₜ = g.el.Nₜ

    # Concentration in electrolyte
    cₑ = el
    logc = NaNMath.log.(cₑ)

    # Electrolyte potential drop (assumed ideal here)
    df_fac = ones(length(cₑ))

    # Derivative using central difference 
    dlogc_dx = [
        (logc[2]-logc[1]) / (x[2]-x[1]),
        [(logc[i+1]-logc[i-1]) / (x[i+1]-x[i-1]) for i in 2:Nₜ-1]...,
        (logc[end]-logc[end-1]) / (x[end]-x[end-1])
    ]

    # Integrand
    f1 = [(1 - params.e.t₊(cₑ[i])) * df_fac[i] * dlogc_dx[i] for i in 1:Nₜ]

    # Inner integral (domain [0, x] for current density)
    int1 = cumsum([
        0,
        [(f1[i]+f1[i-1]) * (x[i]-x[i-1]) / 2 for i in 2:Nₜ]...,
    ])

    # Outer integral (domain on electrodes)
    int2_p = ∫(int1, x[g.el.ixₚ])
    int2_n = ∫(int1, x[g.el.ixₙ])

    # Electrolyte concentration overpotential
    ηₑ = 2*R*T/F*(int2_p/Lₚ - int2_n/Lₙ)
end

# Electrolyte Ohmic loss
function Δϕₑ_f(params::BatteryParameters, g::NamedTuple, cₑ::Symbolics.AbstractArray, ϵ::Symbolics.AbstractArray, i_app)
    # Retrieve FVM geometry
    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    L = sum(g.el.Ls)
    Nₜ = g.el.Nₜ

    # Bruggeman coefficients
    b = [
        [params.e.bₙ for _ in g.el.ixₙ] # Negative electrode
        [params.e.bₛ for _ in g.el.ixₛ] # Separator
        [params.e.bₚ for _ in g.el.ixₚ] # Positive electrode
    ]

    # Transport efficiency (inverse MacMullin number) * Electolyte conductivity
    B = [params.e.σₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:Nₜ]

    # Current in electrolyte
    function iₑ(x)
        # Negative electrode
        if x <= Lₙ
            return i_app*x/Lₙ
        # Separator
        elseif x <= Lₙ + Lₛ
            return i_app
        # Positive electrode
        else
            return i_app*(L-x)/Lₚ
        end
    end
    
    # Integrand
    f1 = [iₑ(x[i]) / B[i] for i in 1:Nₜ]

    # Inner integral (domain [0, x] for current density)
    int1 = cumsum([
        (f1[1] + iₑ(0) / B[1]) * (x[1]) / 2,
        [(f1[i]+f1[i-1]) * (x[i]-x[i-1]) / 2 for i in 2:Nₜ]...,
    ])

    # Outer integral (domain on electrodes)
    int2_p = ∫(int1, x[g.el.ixₚ])
    int2_n = ∫(int1, x[g.el.ixₙ])

    # Electrolyte Ohmic loss
    Δϕₑ = -int2_p/Lₚ + int2_n/Lₙ
end

# Separator Ohmic loss
function Δϕₛ_f(params::BatteryParameters, g::NamedTuple, i_app)
    # Retrieve FVM geometry
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]

    # Separator Ohmic loss
    Δϕₛ = -i_app/3*(Lₚ/params.p.σₖ + Lₙ/params.n.σₖ)
end

# Film Ohmic loss (current only SEI, not CEI)
function Δϕf_f(params::BatteryParameters, g::NamedTuple, i_app)
    # Retrieve FVM geometry
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    Lf = params.n.L_sei₀

    # Film Ohmic loss
    Δϕf = -i_app*(Lf/Lₙ/params.n.aₖ/params.n.σₖ)
end

# Average reaction overpotential
function ηᵣ_x(params::SolidParticleParameters, cₖ, cₑ, j, T)

    @assert length(cₖ) == length(cₑ)

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = length(cₖ)

    jₓ0 = [params.mₖ*sqrt(cₑ[i]*cₖ[i]*(params.c₊-cₖ[i])) for i in 1:N]
    
    asin = [asinh(j/params.aₖ/jₓ0[i]) for i in 1:N]
    
    sinh_x = 2*R*T/F*sum(asin)/N

    return sinh_x
end

# Average electrolyte concentration overpotential
function ηₑ_x(params::ElectrolyteParameters, cₑ, x, T)
    
    @assert length(cₑ) == length(x)
    N = length(cₑ)
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    
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