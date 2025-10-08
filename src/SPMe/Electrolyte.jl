using ModelingToolkit

function Electrolyte(;name, p::ElectrolyteParameters, g)
    @parameters begin
        t # Time variable
    end
    
    @constants begin
        R = 8.314 # Universal gas constant
        F = 96485 # Faraday's constant
        T = 298 # Temperature
    end
    
    Dt = Differential(t)

    @named i_app = RealInput() # Electrolyte current density

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
    end

    function iₑ(x)
        if x <= p.Lₙ
            return i_app.u/p.Lₙ
        elseif x <= p.Lₙ + p.Lₛ
            return 0.0
        else
            return -i_app.u/p.Lₚ
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
    
    Δx = g.Δx
    Δxₗ = g.Δxₗ
    Δxᵣ = g.Δxᵣ
    x = g.x_centers
        
    Dᵢ = [p.Dₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:g.Nₜ] # Face diffusivities
    Dₗ = [nothing, [D_face(Dᵢ[i-1], Dᵢ[i], Δx[i-1], Δx[i]) for i in 2:g.Nₜ]...]
    Dᵣ = [[D_face(Dᵢ[i], Dᵢ[i+1], Δx[i], Δx[i+1]) for i in 1:g.Nₜ-1]..., nothing]

    eqns = [
        # REPLACED Ce[i+1] with Ce[i] to avoid index out of bounds error...
        
        # Boundary condition at the start
        Dt(ϵcₑ[1]) ~ (Dᵣ[1]*(cₑ[2] - cₑ[1])/Δxᵣ[1] + (1-p.t₊(cₑ[1]))*iₑ(x[1])*Δx[1]/F)/Δx[1], 
        
        # Full region
        [Dt(ϵcₑ[i]) ~ 
        (Dᵣ[i]*(cₑ[i+1] - cₑ[i])/Δxᵣ[i] - Dₗ[i]*(cₑ[i] - cₑ[i-1])/Δxₗ[i] + 
        (1-p.t₊(cₑ[i]))*iₑ(x[i])*Δx[i]/F)/Δx[i]
        for i in 2:g.Nₜ-1]...,
        
        # Boundary condition at the end
        Dt(ϵcₑ[end]) ~ (-Dₗ[end]*(cₑ[end] - cₑ[end-1])/Δxₗ[end] + (1-p.t₊(cₑ[end]))*iₑ(x[end])*Δx[end]/F)/Δx[end], 


        # Porosity (assumed constant)
        Dt(ϵₙ) ~ 0,
        Dt(ϵₛ) ~ 0,
        Dt(ϵₚ) ~ 0,
        
        # Actual concentration
        [cₑ[i] ~ ϵcₑ[i]/ϵ[i] for i in 1:g.Nₜ]...,
    ]

    System(eqns, t; name=name,systems=[i_app])
end

function ϕₑ(params::BatteryParameters, g, el, i_app)
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    T = 298 # Temperature

    Nx = g.Nx
    ixₙ = g.ixₙ
    Δx = g.Δx
    Δxₗ = g.Δxₗ
    Δxᵣ = g.Δxᵣ
    x = g.x_centers

    # Bruggeman coefficients per region
    b = [
        [params.e.bₙ for _ in g.ixₙ] 
        [params.e.bₛ for _ in g.ixₛ]
        [params.e.bₚ for _ in g.ixₚ] 
    ]
    # Porosity per region
    ϵ = [
        [params.e.ϵₙ for _ in g.ixₙ] 
        [ϵₛ for _ in g.ixₛ]
        [ϵₚ for _ in g.ixₚ] 
    ]

    iₑ = [
        [i_app.u*x[i]/params.e.Lₙ for i in g.ixₙ]
        [i_app.u for i in g.ixₛ]
        [i_app.u*(g.L - x[i])/params.e.Lₚ for i in g.ixₚ]
    ]

    # Electrolyte potential drop
    thermodynamic_factor = 1
    # central difference 
    dlogc_dx = [
        (log(el.cₑ[2]) - log(el.cₑ[1]))/Δxᵣ[1], # Forward difference at the start
        [(log(el.cₑ[i+1]) - log(el.cₑ[i-1]))/(Δxᵣ[i] + Δxₗ[i]) for i in 2:g.Nₜ-1]...,
        (log(el.cₑ[end]) - log(el.cₑ[end-1]))/Δxₗ[end] # Backward difference at the end
    ]
    int1 = cumsum(iₑ./(params.e.σₑ(el.cₑ).*(el.ϵ.^b)).*Δx)
    int2 = cumsum((1 .- params.e.t₊(el.cₑ)).*thermodynamic_factor.*dlogc_dx.*Δx)
end