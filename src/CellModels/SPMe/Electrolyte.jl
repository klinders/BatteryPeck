using ModelingToolkit

"""
    Electrolyte(; name, p::ElectrolyteParameters, g)

Create a ModelingToolkit system for electrolyte salt concentration and ionic transport.

Models lithium-ion transport through the electrolyte including diffusion in negative electrode,
separator, and positive electrode. Computes concentration profiles and electrochemical potentials.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::ElectrolyteParameters`: Electrolyte material and transport parameters
- `g`: FVM geometry object defining domain and node locations

# Input Ports
- `i_app`: Applied current density (A/m²)
- `T`: Temperature (K)
- `Δϕₙ`: Potential drop in negative electrode (V)
- `ϕₛn`: Solid potential in negative electrode (V)

# Output Variables
- `cₑ`: Electrolyte concentration profile (mol/m³)
- `c̄ₑ`: Average electrolyte concentration (mol/m³)
- `ϕₑ`: Electrolyte potential profile (V)

# Notes
Uses finite volume method with Bruggeman correlation for tortuosity in porous media.
Automatically computes diffusion and migration based on concentration gradients.
"""
function Electrolyte(;name, p::ElectrolyteParameters, g)
    @parameters begin
        t # Time variable
    end

    Δx = g.Δx
    Δxₗ = g.Δxₗ
    Δxᵣ = g.Δxᵣ
    x = g.x_centers
    xₗ = x - Δx/2
    xᵣ = x + Δx/2
    
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    
    Dt = Differential(t)

    @named i_app = RealInput() # Electrolyte current density
    @named j_n = RealInput() # Electrolyte current density in negative electrode
    @named j_p = RealInput() # Electrolyte current density in positive electrode

    @named T = RealInput(guess=298.15)
    @named Δϕₙ = RealInput()
    @named ϕₛn = RealInput()

    @variables begin
        # Electrolyte concentration in mol*m^-3
        (ϵcₑ(t))[1:g.Nₜ] = [
            [p.ϵₙ*p.c₀ for _ in g.ixₙ]...,  # Negative electrode
            [p.ϵₛ*p.c₀ for _ in g.ixₛ]...,  # Separator
            [p.ϵₚ*p.c₀ for _ in g.ixₚ]...,  # Positive electrode
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

        # Actual concentration
        (cₑ(t))[1:g.Nₜ], [guess=fill(p.c₀, g.Nₜ)]
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
            return j_n.u
        elseif x <= p.Lₙ + p.Lₛ
            return 0.0
        else
            return j_p.u
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

    # Bruggeman coefficients per region
    b = [
        [p.bₙ for _ in g.ixₙ] 
        [p.bₛ for _ in g.ixₛ]
        [p.bₚ for _ in g.ixₚ] 
    ]

    # Electrolyte potential drop
    B = [p.σₑ(cₑ[i], T.u)*(ϵ[i]^b[i]) for i in 1:g.Nₜ]
    f1 = [iₑ(x[i])/B[i] for i in 1:g.Nₜ]

    ϕₑ_r = cumsum([
        (f1[1] + iₑ(0)/B[1])*(x[1])/2,
        [(f1[i] + f1[i-1])*(x[i]-x[i-1])/2 for i in 2:g.Nₜ]...,
    ])
    ϕₑ_rn = sum(ϕₑ_r[g.ixₙ])/g.Nx[1]
    ϕₑ_rp = sum(ϕₑ_r[g.ixₚ])/g.Nx[3]

    # Electrolyte reaction potential
    df_fac = ones(length(cₑ))

    Dᵢ = [p.Dₑ(cₑ[i], T.u)*(ϵ[i]^b[i]) for i in 1:g.Nₜ] # Face diffusivities
    Dₗ = [nothing, [D_face(Dᵢ[i-1], Dᵢ[i], Δx[i-1], Δx[i]) for i in 2:g.Nₜ]...]
    Dᵣ = [[D_face(Dᵢ[i], Dᵢ[i+1], Δx[i], Δx[i+1]) for i in 1:g.Nₜ-1]..., nothing]

    Nₗ = [0, [-Dₗ[i]*(cₑ[i] - cₑ[i-1])/Δxₗ[i] + p.t₊(cₑ[i])*iₑ(xₗ[i])/F for i in 2:g.Nₜ]...]
    Nᵣ = [[-Dᵣ[i]*(cₑ[i+1] - cₑ[i])/Δxᵣ[i] + p.t₊(cₑ[i])*iₑ(xᵣ[i])/F for i in 1:g.Nₜ-1]..., 0]

    # Effective electrolyte conductivity
    κₙ = p.σₑ(c̄ₑ, T.u)*(ϵ̄ₙ^p.bₙ)
    κₛ = p.σₑ(c̄ₑ, T.u)*(ϵ̄ₛ^p.bₛ)
    κₚ = p.σₑ(c̄ₑ, T.u)*(ϵ̄ₚ^p.bₚ)

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

        ϕ̄ₑn ~ sum([ϕₑ[i] for i in g.ixₙ])/g.Nx[1],#ϕₑ_rn + χ*R*T.u/F*Mₙ,# 
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
        ηₑ ~ (Mₚ - Mₙ)χ*R*T.u/F
    ]

    # Event working
    events = [
        [
            minimum(cₑ) ~ 0,
            minimum(ϵ) ~ 0,
            maximum(ϵ) ~ 1,
        ]=>(abort!,(;))
    ]

    System(eqns, t; name=name,systems=[i_app, j_n, j_p, T, Δϕₙ, ϕₛn], continuous_events=events)
end