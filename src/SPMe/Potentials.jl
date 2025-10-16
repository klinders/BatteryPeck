function U₀_f(params::BatteryParameters, ne, pe)

    cₚ = pe
    cₙ = ne

    U₀ = params.p.Uₖ(cₚ/params.p.c₊) - params.n.Uₖ(cₙ/params.n.c₊)
end

function ηᵣ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray, ne, pe, i_app)

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    T = 298 # Temperature

    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]

    cₚ = pe
    cₙ = ne
    cₑ = el

    jₚ = params.p.mₖ.*sqrt.(cₑ[g.el.ixₚ].*cₚ.*(params.p.c₊-cₚ))
    jₙ = params.n.mₖ.*sqrt.(cₑ[g.el.ixₙ].*cₙ.*(params.n.c₊-cₙ))

    asin_p = asinh.(i_app./params.p.aₖ./Lₚ./jₚ)
    asin_n = asinh.(i_app./params.n.aₖ./Lₙ./jₙ)

    sum_p = ∫(asin_p, x[g.el.ixₚ])
    sum_n = ∫(asin_n, x[g.el.ixₙ])

    ηᵣ = -2*R*T/F*(sum_p/Lₚ + sum_n/Lₙ)

end


function ηₑ_f(params::BatteryParameters, g::NamedTuple, el::Symbolics.AbstractArray)

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    T = 298 # Temperature

    x = g.el.x_centers
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
    Nₜ = g.el.Nₜ

    cₑ = el

    # Electrolyte potential drop
    df_fac = ones(length(cₑ))

    logc = log.(cₑ)

    # central difference 
    dlogc_dx = [
        (logc[2]-logc[1])/(x[2]-x[1]),
        [(logc[i+1]-logc[i-1])/(x[i+1]-x[i-1]) for i in 2:Nₜ-1]...,
        (logc[end]-logc[end-1])/(x[end]-x[end-1])
    ]

    f1 = [(1 - params.e.t₊(cₑ[i]))*df_fac[i]*dlogc_dx[i] for i in 1:Nₜ]

    int1 = cumsum([
        0,
        [(f1[i] + f1[i-1])*(x[i]-x[i-1])/2 for i in 2:Nₜ]...,
    ])

    int1_p = ∫(int1, x[g.el.ixₚ])
    int1_n = ∫(int1, x[g.el.ixₙ])

    ηₑ = 2*R*T/F*(int1_p/Lₚ - int1_n/Lₙ)

end


function Δϕₑ_f(params::BatteryParameters, g::NamedTuple, cₑ::Symbolics.AbstractArray, ϵ::Symbolics.AbstractArray, i_app)

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

    B = [params.e.σₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:Nₜ]

    function iₑ(x)
        if x <= Lₙ
            return i_app*x/Lₙ
        elseif x <= Lₙ + Lₛ
            return i_app
        else
            return i_app*(L-x)/Lₚ
        end
    end

    f2 = [iₑ(x[i])/B[i] for i in 1:Nₜ]

    int2 = cumsum([
        (f2[1] + iₑ(0)/B[1])*(x[1])/2,
        [(f2[i] + f2[i-1])*(x[i]-x[i-1])/2 for i in 2:Nₜ]...,
    ])

    int2_p = ∫(int2, x[g.el.ixₚ])
    int2_n = ∫(int2, x[g.el.ixₙ])

    Δϕₑ = -int2_p/Lₚ + int2_n/Lₙ

end

function Δϕₛ_f(params::BatteryParameters, g::NamedTuple, i_app)

    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]

    Δϕₛ = -i_app/3*(Lₚ/params.p.σₖ + Lₙ/params.n.σₖ)

end

function Δϕf_f(params::BatteryParameters, g::NamedTuple, i_app)

    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]

    Lf = params.n.L_sei₀
    Δϕf = -i_app*(Lf/Lₙ/params.n.aₖ/params.n.σₖ)

end
