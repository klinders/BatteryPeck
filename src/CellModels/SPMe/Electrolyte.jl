# =====================================================================================================================
# Electrolyte.jl
#
# Electrolyte equations for SPMe
# Source: https://doi.org/10.1016/j.apm.2022.12.009
# =====================================================================================================================

# Import package
using ModelingToolkit

function Electrolyte(;name, p::ElectrolyteParameters, g)
    # Independent variables
    @parameters begin
        t # Time
    end

    # Retrievve FVM geometry
    Δx = g.Δx
    Δxₗ = g.Δxₗ
    Δxᵣ = g.Δxᵣ
    x = g.x_centers
    xₗ = x - Δx/2
    xᵣ = x + Δx/2
    
    # Time derivative
    Dt = Differential(t)
    
    # Constants
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant

    # Components
    @named i_app = RealInput()          # Electrolyte current density
    @named T = RealInput(guess=298.15)              # Temperature
    @named Δϕₙ = RealInput(guess=0.0)
    @named ϕₛn = RealInput(guess=0.0)

    # Time-dependent state variables
    @variables begin
        # Eq. 4a LHS; [Porosity] * [Electrolyte concentration]
        # Subdivide into total number of subdivisions over all regions to match FVM cells
        (ϵcₑ(t))[1:g.Nₜ] = [
            # Scalar regions (single porosity scalar value per region)
            # Initialise at c₀ (Eq. 4c)
            [p.ϵₙ*p.c₀ for _ in g.ixₙ]..., # Negative electrode
            [p.ϵₛ*p.c₀ for _ in g.ixₛ]..., # Separator
            [p.ϵₚ*p.c₀ for _ in g.ixₚ]..., # Positive electrode
        ]

        # Porosity
        (ϵ(t))[1:g.Nₜ] 
        # = [
        #     [p.ϵₙ for _ in g.ixₙ]...,       # Negative electrode
        #     [p.ϵₛ for _ in g.ixₛ]...,       # Separator
        #     [p.ϵₚ for _ in g.ixₚ]...,       # Positive electrode
        # ]
        ϵ̄ₙ(t) 
        ϵ̄ₛ(t) 
        ϵ̄ₚ(t)

        # Define concentration as symbolic state variable (calculated in "eqns" below)
        (cₑ(t))[1:g.Nₜ]
        # X average concentration
        c̄ₑ(t)
        c̄ₑn(t)
        c̄ₑs(t)
        c̄ₑp(t)

        # Electrolyte potential
        (ϕₑ(t))[1:g.Nₜ]
        ϕ̄ₑ(t)
        ϕ̄ₑn(t), [guess=0]
        ϕ̄ₑs(t)
        ϕ̄ₑp(t)

        Δϕₑ(t)
        ηₑ(t)

        # Electrolyte heat generation
        Qₑ(t)
    end

    # Eq. 54; Spatial derivative of applied current density in electrolyte
    function j(x)
        # Negative electrode
        if x <= p.Lₙ
            return i_app.u/p.Lₙ
        # Separator
        elseif x <= p.Lₙ + p.Lₛ
            return 0.0
        # Positive electrode
        else
            return -i_app.u/p.Lₚ
        end
    end

    L = sum(g.Ls)
    function iₑ(x)
        # Negative electrode
        if x <= p.Lₙ
            return i_app.u*x/p.Lₙ
        # Separator
        elseif x <= p.Lₙ + p.Lₛ
            return i_app.u
        # Positive electrode
        else
            return i_app.u*(L-x)/p.Lₚ
        end
    end

    # Concatenate Bruggeman coefficients vectors
    b = [
        [p.bₙ for _ in g.ixₙ] # Negative electrode
        [p.bₛ for _ in g.ixₛ] # Separator
        [p.bₚ for _ in g.ixₚ] # Positive electrode
    ]

    # Electrolyte potential drop
    B = [p.σₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:g.Nₜ]
    f1 = [iₑ(x[i])/B[i] for i in 1:g.Nₜ]

    ϕₑ_r = cumsum([
        (f1[1] + iₑ(0)/B[1])*(x[1])/2,
        [(f1[i] + f1[i-1])*(x[i]-x[i-1])/2 for i in 2:g.Nₜ]...,
    ])
    ϕₑ_rn = sum(ϕₑ_r[g.ixₙ])/g.Nx[1]
    ϕₑ_rp = sum(ϕₑ_r[g.ixₚ])/g.Nx[3]

    # Electrolyte reaction potential
    df_fac = ones(length(cₑ))

    # Eq. 4a RHS; Effective diffusivity (transport efficiency, a.k.a. inverse MacMullin number, calculated using Bruggeman coefficient)
    Dᵢ = [p.Dₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:g.Nₜ]
    # Face diffusivities (see "helpers.jl" for explanation and source)
    # Outer most cells have only one neighbour
    Dₗ = [nothing, [D_face(Dᵢ[i-1], Dᵢ[i], Δx[i-1], Δx[i]) for i in 2:g.Nₜ]...] # Left neighbours
    Dᵣ = [[D_face(Dᵢ[i], Dᵢ[i+1], Δx[i], Δx[i+1]) for i in 1:g.Nₜ-1]..., nothing] # Right neighbours

    Nₗ = [0, [-Dₗ[i]*(cₑ[i] - cₑ[i-1])/Δxₗ[i] + p.t₊(cₑ[i])*iₑ(xₗ[i])/F for i in 2:g.Nₜ]...]
    Nᵣ = [[-Dᵣ[i]*(cₑ[i+1] - cₑ[i])/Δxᵣ[i] + p.t₊(cₑ[i])*iₑ(xᵣ[i])/F for i in 1:g.Nₜ-1]..., 0]

    # Effective electrolyte conductivity
    κₙ = p.σₑ(c̄ₑ)*(ϵ̄ₙ^p.bₙ)
    κₛ = p.σₑ(c̄ₑ)*(ϵ̄ₛ^p.bₛ)
    κₚ = p.σₑ(c̄ₑ)*(ϵ̄ₚ^p.bₚ)

    # phi_e max 1e15
    function M(x)
        tol = 1e-15
        x = max(x, tol)
        return log(x)
    end

    χ = 2*(1 - p.t₊(c̄ₑ))*df_fac[1]

    Mₙ = sum([M(cₑ[i]/c̄ₑ) for i in g.ixₙ])/g.Nx[1]
    Mₚ = sum([M(cₑ[i]/c̄ₑ) for i in g.ixₚ])/g.Nx[3]
    ϕₑ_const = -Δϕₙ.u + ϕₛn.u - χ*R*T.u/F*Mₙ - i_app.u*p.Lₙ*(1/(3*κₙ) - 1/κₛ)

    ϕi = [
        [i_app.u/κₙ*(x[i]^2 - p.Lₙ^2)/(2*p.Lₙ) + i_app.u*p.Lₙ/κₛ for i in g.ixₙ]
        [i_app.u/κₛ*x[i] for i in g.ixₛ]
        [i_app.u/κₚ*(x[i]*(2*L - x[i]) + p.Lₚ^2 - L^2)/(2*p.Lₚ) + i_app.u*(L - p.Lₚ)/κₛ for i in g.ixₚ]
    ]

    ϕₑ_f = [ϕₑ_const + χ*R*T.u/F*M(cₑ[i]/c̄ₑ) - ϕi[i] for i in 1:g.Nₜ]

    # Calculate Ohmic heat integrand using current and conductivity
    f_ohm = [iₑ(x[i])^2 / B[i] for i in 1:g.Nₜ]

    # Calculate spatial integral of Ohmic heat across cell
    int_ohm = sum([
        (f_ohm[1] + iₑ(0)^2 / B[1]) * x[1] / 2,
        [(f_ohm[i] + f_ohm[i-1]) * (x[i] - x[i-1]) / 2 for i in 2:g.Nₜ]...,
        f_ohm[end] * (L - x[end]) / 2
    ])

    # Average Ohmic heat over total length
    q_ohm = int_ohm / L

    eqns = [
        # REPLACED Ce[i+1] with Ce[i] to avoid index out of bounds error...
        
        # Boundary condition at the start
        # Dt(ϵcₑ[1]) ~ (Dᵣ[1]*(cₑ[2] - cₑ[1])/Δxᵣ[1] + (1-p.t₊(cₑ[1]))*j(x[1])*Δx[1]/F)/Δx[1], 
        
        # Full region
        # [Dt(ϵcₑ[i]) ~ 
        # (Dᵣ[i]*(cₑ[i+1] - cₑ[i])/Δxᵣ[i] - Dₗ[i]*(cₑ[i] - cₑ[i-1])/Δxₗ[i] + 
        # (1-p.t₊(cₑ[i]))*j(x[i])*Δx[i]/F)/Δx[i]
        # for i in 2:g.Nₜ-1]...,
        [Dt(ϵcₑ[i]) ~ -(Nᵣ[i] - Nₗ[i])/Δx[i] + j(x[i])/F for i in 1:g.Nₜ]...,
        

        # Boundary condition at the end
        # Dt(ϵcₑ[end]) ~ (-Dₗ[end]*(cₑ[end] - cₑ[end-1])/Δxₗ[end] + (1-p.t₊(cₑ[end]))*j(x[end])*Δx[end]/F)/Δx[end], 
        
        # Actual concentration
        [cₑ[i] ~ ϵcₑ[i]/ϵ[i] for i in 1:g.Nₜ]...,
        c̄ₑn ~ sum([cₑ[i] for i in g.ixₙ])/g.Nx[1],
        c̄ₑs ~ sum([cₑ[i] for i in g.ixₛ])/g.Nx[2],
        c̄ₑp ~ sum([cₑ[i] for i in g.ixₚ])/g.Nx[3],
        # Weighted average of the electrodes
        c̄ₑ ~ (c̄ₑn*p.Lₙ + c̄ₑs*p.Lₛ + c̄ₑp*p.Lₚ)/L,

        # Electrolyte potential
        [ϕₑ[i] ~ ϕₑ_f[i] for i in 1:g.Nₜ]...,

        ϕ̄ₑn ~ ϕₑ_rn + χ*R*T.u/F*Mₙ,# sum([ϕₑ[i] for i in g.ixₙ])/g.Nx[1],
        ϕ̄ₑs ~ sum([ϕₑ[i] for i in g.ixₛ])/g.Nx[2],
        ϕ̄ₑp ~ sum([ϕₑ[i] for i in g.ixₚ])/g.Nx[3],
        ϕ̄ₑ ~ (ϕ̄ₑn*p.Lₙ + ϕ̄ₑs*p.Lₛ + ϕ̄ₑp*p.Lₚ)/L,

        # Porosity
        ϵ̄ₙ ~ sum([ϵ[i] for i in g.ixₙ])/g.Nx[1],
        ϵ̄ₛ ~ sum([ϵ[i] for i in g.ixₛ])/g.Nx[2],
        ϵ̄ₚ ~ sum([ϵ[i] for i in g.ixₚ])/g.Nx[3],

        # Constant
        # [Dt(ϵ[i]) ~ 0 for i in 1:g.Nₜ]...,
        # Marquis 2019
        Δϕₑ ~ -i_app.u*(p.Lₙ/(3*κₙ) + p.Lₛ/(κₛ) + p.Lₚ/(3*κₚ)),
        ηₑ ~ (Mₚ - Mₙ)χ*R*T.u/F,

        # Total electrolyte heat combines Ohmic and concentration heat
        Qₑ ~ q_ohm - i_app.u * ηₑ / L
    ]

    # Construct ODESystem with equations and child components
    System(eqns, t; name=name,systems=[i_app, T, Δϕₙ, ϕₛn])
end