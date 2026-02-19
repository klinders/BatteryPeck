# =====================================================================================================================
# Heat.jl
#
# Generated heat sources using SPMe
# Source[1]: https://doi.org/10.1016/j.apm.2022.12.009        # SPMe
# Source[2]: https://doi.org/10.1016/j.electacta.2021.138524  # Thermal equations for SPMe
# Source[3]: https://doi.org/10.1016/j.electacta.2022.140700  # EHC curves
#            https://zenodo.org/records/5171874               # .csv files
#
# Equations based on Potentials.jl
# =====================================================================================================================

# Import package (for modelling EHC)
using Symbolics

# Heat generated from electrolyte (~Δϕₑ*i + ηₑ*i)
function Qₑ_f(params::BatteryParameters, g::NamedTuple, cₑ::Symbolics.AbstractArray, ϵ::Symbolics.AbstractArray, i_app, ηₑ)
    # Retrieve FVM geometry
    x = g.el.x_centers
    L = sum(g.el.Ls)
    Lₙ,Lₛ,Lₚ = g.el.Ls[1], g.el.Ls[2], g.el.Ls[3]
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

    # Integrand (current squared)
    f1 = [iₑ(x[i])^2 / B[i] for i in 1:Nₜ]

    # Integral (sum instead of cumsum as total heat is scalar, not vector)
    int1 = sum([
        (f1[1] + iₑ(0)^2 / B[1]) * x[1] / 2,                   # Start > Node 1
        [(f1[i]+f1[i-1]) * (x[i]-x[i-1]) / 2 for i in 2:Nₜ]..., # Nodes
        (f1[end]+0) * (L - x[end]) / 2                         # Last node > end
    ])

    # Ohmic heat
    q_ohm = int1 / L

    # Concentration heat
    q_conc = -i_app * ηₑ / L

    # Total electrolyte heat (Eq. 5d [2], Eq. 13 [1])
    Qₑ = q_ohm + q_conc
end

# Irreversible heat (~ηᵣ*i)
function Qᵢ_f(params::BatteryParameters, g::NamedTuple, i_app, ηᵣ)
    # Total electrode length
    L = sum(g.el.Ls)

    # Irreversible heat (Eq. 5e [2], Eq. 13 [1])
    Qᵢ = -i_app * ηᵣ / L
end

# Solid phase (electrodes) Ohmic Heat Generation (~ϕₛ*i)
function Qₛ_f(params::BatteryParameters, g::NamedTuple, i_app, Δϕₛ)
    # Total electrode length
    L = sum(g.el.Ls)

    # Solid phase Ohmic Heat Generation (Eq. 5c [2], Eq. 13 [1])
    Qₛ = -i_app * Δϕₛ / L
end

# Film Ohmic heat generation (~ϕf*i)
function Qf_f(params::BatteryParameters, g::NamedTuple, i_app, Δϕf)
    # Total electrode length
    L = sum(g.el.Ls)

    # Film Ohmic heat generation (Eq. 13 [1])
    Qf = -i_app * Δϕf / L
end

# Polynomial fit parameters entropic term positive electrode (NMC811) (Tab. S7 [3])
const P_coeff_pos = (
    a1 = 0.04006,
    b1 = 0.2828,
    c1 = 0.0009855,
    a2 = -0.06656,
    b2 = 0.8032,
    c2 = 0.02179
)

# Polynomial fit parameters parameters entropic term negative electrode (Graphite-SiOx) (Tab. S7 [3])
const P_coeff_neg = (
    a0 = -0.111,
    b0 = 0.02901,
    a1 = 0.3562,
    b1 = 0.08308,
    c1 = 0.004621
)

# Entropic term positive electrode (Eqs. 16 [3])
function dUp_dT_f(z)
    p = P_coeff_pos
    # Coefficients in mV/K, multiply by 1e-3 to get V/K
    val_mV = p.a1 * exp(-((z - p.b1)^2) / p.c1) + 
             p.a2 * exp(-((z - p.b2)^2) / p.c2)
    return val_mV * 1e-3
end

# Entropic term positive electrode negative electrode (Eq. 17 [3])
function dUn_dT_f(z)
    p = P_coeff_neg
    # Coefficients in mV/K, multiply by 1e-3 to get V/K
    val_mV = p.a0 * z + p.b0 + 
             p.a1 * exp(-((z - p.b1)^2) / p.c1)
    return val_mV * 1e-3
end

# Reversible (entropic) heat generation 
function Q_rev_f(params::BatteryParameters, g::NamedTuple, i_app, T, c_s_n_surf, c_s_p_surf)
    # Retrieve geometry
    L = sum(g.el.Ls)
    
    # Surface stoichiometry
    z_n = c_s_n_surf / params.n.c₊
    z_p = c_s_p_surf / params.p.c₊

    # Entropic coefficients
    dUn = dUn_dT_f(z_n)
    dUp = dUp_dT_f(z_p)

    # Peltier coefficients
    Pi_n = T * dUn
    Pi_p = T * dUp

    # Reversible (entropic) heat generation (Eq. 5f [2])
    Q_rev = (i_app / L) * (Pi_n - Pi_p)
    
    return Q_rev
end