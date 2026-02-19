using ModelingToolkit

function Electrolyte(;name, p::ElectrolyteParameters, g)
    @parameters begin
        t # Time variable
    end

    Δx = g.Δx
    Δxₗ = g.Δxₗ
    Δxᵣ = g.Δxᵣ
    x = g.x_centers
    
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    
    Dt = Differential(t)

    @named i_app = RealInput() # Electrolyte current density
    @named T = RealInput()
    @named Δϕₙ = RealInput(guess=0.0)
    @named ϕₛn = RealInput(guess=0.0)

    @variables begin
        # Electrolyte concentration in mol*m^-3
        (ϵcₑ(t))[1:g.Nₜ] = [
            [p.ϵₙ*p.c₀ for _ in g.ixₙ]...,  # Negative electrode
            [p.ϵₛ*p.c₀ for _ in g.ixₛ]...,  # Separator
            [p.ϵₚ*p.c₀ for _ in g.ixₚ]...,  # Positive electrode
        ]

        # Porosity
        ϵₙ(t) = p.ϵₙ # Negative electrode
        ϵₛ(t) = p.ϵₛ # Separator
        ϵₚ(t) = p.ϵₚ # Positive electrode

        # Actual concentration
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
    end

    function j(x)
        if x <= p.Lₙ
            return i_app.u/p.Lₙ
        elseif x <= p.Lₙ + p.Lₛ
            return 0.0
        else
            return -i_app.u/p.Lₚ
        end
    end

    L = sum(g.Ls)
    function iₑ(x)
        if x <= p.Lₙ
            return i_app.u*x/p.Lₙ
        elseif x <= p.Lₙ + p.Lₛ
            return i_app.u
        else
            return i_app.u*(L-x)/p.Lₚ
        end
    end

    # Porosity per region
    ϵ = [
        [ϵₙ for _ in g.ixₙ] 
        [ϵₛ for _ in g.ixₛ]
        [ϵₚ for _ in g.ixₚ] 
    ]

    # Bruggeman coefficients per region
    b = [
        [p.bₙ for _ in g.ixₙ] 
        [p.bₛ for _ in g.ixₛ]
        [p.bₚ for _ in g.ixₚ] 
    ]

    # Electrolyte potential drop
    B = [p.σₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:g.Nₜ]
    f1 = [iₑ(x[i])/B[i] for i in 1:g.Nₜ]

    ϕₑ_r = cumsum([
        (f1[1] + iₑ(0)/B[1])*(x[1])/2,
        [(f1[i] + f1[i-1])*(x[i]-x[i-1])/2 for i in 2:g.Nₜ]...,
    ])

    # Electrolyte reaction potential
    df_fac = ones(length(cₑ))

    Dᵢ = [p.Dₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:g.Nₜ] # Face diffusivities
    Dₗ = [nothing, [D_face(Dᵢ[i-1], Dᵢ[i], Δx[i-1], Δx[i]) for i in 2:g.Nₜ]...]
    Dᵣ = [[D_face(Dᵢ[i], Dᵢ[i+1], Δx[i], Δx[i+1]) for i in 1:g.Nₜ-1]..., nothing]

    # Other approach
    κₙ = p.σₑ(c̄ₑ)*(ϵₙ^p.bₙ)
    κₛ = p.σₑ(c̄ₑ)*(ϵₛ^p.bₛ)
    κₚ = p.σₑ(c̄ₑ)*(ϵₚ^p.bₚ)

    # phi_e max 1e15
    M = [log(max(cₑ[i]/p.c₀,1e-15)) for i in 1:g.Nₜ]
    M_n = sum(M[g.ixₙ])/g.Nx[1]
    M_p = sum(M[g.ixₚ])/g.Nx[3]
    ie_n = sum([ϕₑ_r[i] for i in g.ixₙ])/g.Nx[1]

    ϕₑ_const = -Δϕₙ.u + ϕₛn.u - (1 - p.t₊(c̄ₑn))*df_fac[1]*M_n*2*R*T.u/F - i_app.u*p.Lₙ*(1/(3*κₙ) - 1/κₛ)

    ϕi = [
        [i_app.u/κₙ*(x[i]^2 - p.Lₙ^2)/(2*p.Lₙ) + i_app.u*p.Lₙ/κₛ for i in g.ixₙ]
        [i_app.u/κₛ*x[i] for i in g.ixₛ]
        [i_app.u/κₚ*(x[i]*(2*L - x[i]) + p.Lₚ^2 - L^2)/(2*p.Lₚ) + i_app.u*(L - p.Lₚ)/κₛ for i in g.ixₚ]
    ]

    ϕₑ_f = [ϕₑ_const + (1 - p.t₊(c̄ₑ))*df_fac[i]*2*R*T.u/F*M[i] - ϕi[i] for i in 1:g.Nₜ]

    eqns = [
        # REPLACED Ce[i+1] with Ce[i] to avoid index out of bounds error...
        
        # Boundary condition at the start
        Dt(ϵcₑ[1]) ~ (Dᵣ[1]*(cₑ[2] - cₑ[1])/Δxᵣ[1] + (1-p.t₊(cₑ[1]))*j(x[1])*Δx[1]/F)/Δx[1], 
        
        # Full region
        [Dt(ϵcₑ[i]) ~ 
        (Dᵣ[i]*(cₑ[i+1] - cₑ[i])/Δxᵣ[i] - Dₗ[i]*(cₑ[i] - cₑ[i-1])/Δxₗ[i] + 
        (1-p.t₊(cₑ[i]))*j(x[i])*Δx[i]/F)/Δx[i]
        for i in 2:g.Nₜ-1]...,
        
        # Boundary condition at the end
        Dt(ϵcₑ[end]) ~ (-Dₗ[end]*(cₑ[end] - cₑ[end-1])/Δxₗ[end] + (1-p.t₊(cₑ[end]))*j(x[end])*Δx[end]/F)/Δx[end], 


        # Porosity (assumed constant)
        Dt(ϵₙ) ~ 0,
        Dt(ϵₛ) ~ 0,
        Dt(ϵₚ) ~ 0,
        
        # Actual concentration
        [cₑ[i] ~ ϵcₑ[i]/ϵ[i] for i in 1:g.Nₜ]...,
        c̄ₑ ~ sum(cₑ)/g.Nₜ,
        c̄ₑn ~ sum([cₑ[i] for i in g.ixₙ])/g.Nx[1],
        c̄ₑs ~ sum([cₑ[i] for i in g.ixₛ])/g.Nx[2],
        c̄ₑp ~ sum([cₑ[i] for i in g.ixₚ])/g.Nx[3],

        # Electrolyte potential
        [ϕₑ[i] ~ ϕₑ_f[i] for i in 1:g.Nₜ]...,

        ϕ̄ₑ ~ sum(ϕₑ)/g.Nₜ,
        ϕ̄ₑn ~ sum([ϕₑ[i] for i in g.ixₙ])/g.Nx[1],
        ϕ̄ₑs ~ sum([ϕₑ[i] for i in g.ixₛ])/g.Nx[2],
        ϕ̄ₑp ~ sum([ϕₑ[i] for i in g.ixₚ])/g.Nx[3],

        # Marquis 2019
        Δϕₑ ~ -i_app.u*(p.Lₙ/(3*κₙ) + p.Lₛ/κₛ + p.Lₚ/(3*κₚ)),
        ηₑ ~ (1 - p.t₊(c̄ₑ))*(M_p - M_n)*2*R*T.u/F
    ]

    System(eqns, t; name=name,systems=[i_app, T, Δϕₙ, ϕₛn])
end