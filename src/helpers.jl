# Harmonic mean for face diffusivity
# From:
# http://dx.doi.org/10.1149/2.0291607jes
function D_face(Dleft, Dright, Δxleft, Δxright)
    # Harmonic mean of left and right diffusivities
    return (Δxleft + Δxright)/(Δxleft/Dleft + Δxright/Dright)
end

# Geometry builder for 3-region 1D FVM (cell-centered)
function build_fvm_geometry(params::BatteryParameters, N::Dict{Symbol,Vector{Int}})

    # Solid parameters
    function solid_geometry(params::SolidParticleParameters, Nᵣ)
        Δr = params.Rₖ/Nᵣ
        r = ([Δr*(i-0.5) for i in 1:Nᵣ]) # Radial positions
        Vᵢ = 4/3*π*[(r[i] + Δr/2)^3 - (r[i] - Δr/2)^3 for i in 1:Nᵣ]
        Aₗ = 4*π*[(r[i] - Δr/2)^2 for i in 1:Nᵣ]
        Aᵣ = 4*π*[(r[i] + Δr/2)^2 for i in 1:Nᵣ]
        return (Nᵣ=Nᵣ, Δr=Δr, r=r, Vᵢ=Vᵢ, Aₗ=Aₗ, Aᵣ=Aᵣ)
    end
    ne = solid_geometry(params.n, N[:Nᵣ][1])
    pe = solid_geometry(params.p, N[:Nᵣ][2])

    Nx=N[:Nₓ]
    Nₜ=sum(Nx) # Total number of cells
    Ls=[params.e.Lₙ, params.e.Lₛ, params.e.Lₚ]
    L=sum(Ls)
    ixₙ=1:Nx[1]
    ixₛ=(Nx[1]+1):(Nx[1]+Nx[2])
    ixₚ=(Nx[1]+Nx[2]+1):(Nx[1]+Nx[2]+Nx[3])
    Δxₙ=[Ls[1]/Nx[1] for _ in ixₙ]  # Negative electrode
    Δxₛ=[Ls[2]/Nx[2] for _ in ixₛ]  # Separator
    Δxₚ=[Ls[3]/Nx[3] for _ in ixₚ]
    Δx =[Δxₙ..., Δxₛ..., Δxₚ...]  # All cells
    x_centers = [
        [Δxₙ[i]*(i-0.5) for i in 1:Nx[1]]...,  # Negative electrode
        [Ls[1] + Δxₛ[i]*(i-0.5) for i in 1:Nx[2]]...,  # Separator
        [Ls[1] + Ls[2] + Δxₚ[i]*(i-0.5) for i in 1:Nx[3]]...,  # Positive electrode
    ]
    Δxₗ = [nothing, [Δx[i-1]/2 + Δx[i]/2 for i in 2:Nₜ]...]
    Δxᵣ = [[Δx[i]/2 + Δx[i+1]/2 for i in 1:Nₜ-1]..., nothing]

    el = (Nx=Nx, Nₜ=Nₜ, Ls=Ls, L=L, ixₙ=ixₙ, ixₛ=ixₛ, ixₚ=ixₚ,
          Δx=Δx, Δxₗ=Δxₗ, Δxᵣ=Δxᵣ, x_centers=x_centers)

    return (ne=ne, pe=pe, el=el)
end