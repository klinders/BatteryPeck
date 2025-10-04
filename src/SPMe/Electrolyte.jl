using ModelingToolkit
include("../ParameterSets/Base.jl")


# Harmonic mean for face diffusivity
# From:
# http://dx.doi.org/10.1149/2.0291607jes
function D_face(Dleft, Dright, Δxleft, Δxright)
    # Harmonic mean of left and right diffusivities
    return (Δxleft + Δxright)/(Δxleft/Dleft + Δxright/Dright)
end

# Geometry builder for 3-region 1D FVM (cell-centered)
function build_fvm_geometry(Nx::Vector{Int}, Ls::Vector{Float64}, x0::Float64 = 0.0)
    @assert length(Nx) == length(Ls) == 3

    # number of cells total
    N = sum(Nx)

    # per-region cell widths (uniform within region)
    Δx_reg = [Ls[k] / Nx[k] for k in 1:3]

    # cell centers and Δx per cell
    x_centers = Float64[]
    Δx_cell   = Float64[]
    region_idx = Int[]  # store region id per cell (1,2,3)

    x_cursor = x0
    for k in 1:3
        dxk = Δx_reg[k]
        # centers at x_cursor + (j-0.5)*dxk for j=1..Nx[k]
        for j in 1:Nx[k]
            push!(x_centers, x_cursor + (j - 0.5)*dxk)
            push!(Δx_cell, dxk)
            push!(region_idx, k)
        end
        x_cursor += Ls[k]
    end

    # face distances (distance between adjacent centers)
    # faces exist between cell i and i+1 for i=1..N-1
    dx_face = [Δx_cell[i]/2 + Δx_cell[i+1]/2 for i in 1:(N-1)]

    # index ranges per region in global indexing
    ranges = []
    start = 1
    for k in 1:3
        rng = start:(start + Nx[k] - 1)
        push!(ranges, rng)
        start += Nx[k]
    end

    return (Nx=Nx, Ls=Ls, x_centers=x_centers, Δx_cell=Δx_cell, dx_face=dx_face, ranges=ranges, Δx_reg=Δx_reg, region_idx=region_idx)
end

function Electrolyte(;name, p::ElectrolyteParameters, Nₓ::Vector{Int})
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
    
    geometry = build_fvm_geometry(Nₓ, [p.Lₙ, p.Lₛ, p.Lₚ])
    Nₜ = sum(Nₓ) # Total number of cells
    Ls = geometry.Ls
    L = sum(Ls) # Total length

    ixₙ = geometry.ranges[1]
    ixₛ = geometry.ranges[2]
    ixₚ = geometry.ranges[3]

    @variables begin
        # Electrolyte concentration in mol*m^-3
        (ϵcₑ(t))[1:Nₜ] = [
            [p.ϵₙ*p.c₀ for _ in ixₙ]...,  # Negative electrode
            [p.ϵₛ*p.c₀ for _ in ixₛ]...,  # Separator
            [p.ϵₚ*p.c₀ for _ in ixₚ]...,  # Positive electrode
        ]

        # Porosity
        ϵₙ(t) = p.ϵₙ # Negative electrode
        ϵₛ(t) = p.ϵₛ # Separator
        ϵₚ(t) = p.ϵₚ # Positive electrode

        # Actual concentration
        (cₑ(t))[1:Nₜ]
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
        [ϵₙ for _ in ixₙ] 
        [ϵₛ for _ in ixₛ]
        [ϵₚ for _ in ixₚ] 
    ]

    # Bruggeman coefficients per region
    b = [
        [p.bₙ for _ in ixₙ] 
        [p.bₛ for _ in ixₛ]
        [p.bₚ for _ in ixₚ] 
    ]
    
    Dᵢ = [p.Dₑ(cₑ[i])*(ϵ[i]^b[i]) for i in 1:Nₜ] # Face diffusivities
    Dₗ = [nothing, [D_face(Dᵢ[i-1], Dᵢ[i], geometry.Δx_cell[i-1], geometry.Δx_cell[i]) for i in 2:Nₜ]...]
    Dᵣ = [[D_face(Dᵢ[i], Dᵢ[i+1], geometry.Δx_cell[i], geometry.Δx_cell[i+1]) for i in 1:Nₜ-1]..., nothing]
    Δxₗ = [nothing, [geometry.Δx_cell[i-1]/2 + geometry.Δx_cell[i]/2 for i in 2:Nₜ]...]
    Δxᵣ = [[geometry.Δx_cell[i]/2 + geometry.Δx_cell[i+1]/2 for i in 1:Nₜ-1]..., nothing]
    Δx = geometry.Δx_cell
    x = geometry.x_centers

    eqns = [
        # REPLACED Ce[i+1] with Ce[i] to avoid index out of bounds error...
        
        # Boundary condition at the start
        Dt(ϵcₑ[1]) ~ (Dᵣ[1]*(cₑ[2] - cₑ[1])/Δxᵣ[1] + (1-p.t₊(cₑ[1]))*iₑ(x[1])*Δx[1]/F)/Δx[1], 
        
        # Full region
        [Dt(ϵcₑ[i]) ~ 
        (Dᵣ[i]*(cₑ[i+1] - cₑ[i])/Δxᵣ[i] - Dₗ[i]*(cₑ[i] - cₑ[i-1])/Δxₗ[i] + 
        (1-p.t₊(cₑ[i]))*iₑ(x[i])*Δx[i]/F)/Δx[i]
        for i in 2:Nₜ-1]...,
        
        # Boundary condition at the end
        Dt(ϵcₑ[end]) ~ (-Dₗ[end]*(cₑ[end] - cₑ[end-1])/Δxₗ[end] + (1-p.t₊(cₑ[end]))*iₑ(x[end])*Δx[end]/F)/Δx[end], 


        # Porosity (assumed constant)
        Dt(ϵₙ) ~ 0,
        Dt(ϵₛ) ~ 0,
        Dt(ϵₚ) ~ 0,
        
        # Actual concentration
        [cₑ[i] ~ ϵcₑ[i]/ϵ[i] for i in 1:Nₜ]...,
    ]

    System(eqns, t; name=name,systems=[i_app])
end