# =====================================================================================================================
# helpers.jl
#
# Auxilliary functions:
# - Harmonic mean (calculate face diffusivity between two regions with different diffusivities)
# - Geometry builder for FVM (calculate cell centres, volumes, areas, etc.)
# - Trapezoidal integration
# =====================================================================================================================

# Harmonic mean for face diffusivity
# Derived from addition of "diffusion resistances", where [R = Δx/D] and [R_eff = R_left + R_right]
# Arithmetic mean would give non-zero diffusivity at a face neighbouring a [D = 0] region
# Source: http://dx.doi.org/10.1149/2.0291607jes
function D_face(Dleft, Dright, Δxleft, Δxright)

    # Harmonic mean of left and right diffusivities
    return (Δxleft + Δxright)/(Δxleft/Dleft + Δxright/Dright)
end

# Geometry builder for 3-region 1D FVM (cell-centered)
# Inputs: ["Battery parameters"; "Number of cells per region"]
# params used to extract region lengths, N defined as dictionary consisting of symbol-vector pairs, where vector defines the number of cells in each region
function build_fvm_geometry(params::BatteryParameters, N::Dict{Symbol,Vector{Int}})

    # Determine FVM parameters for solid particles
    # Inputs: ["Particle parameters"; "Number of cells per region"]
    function solid_geometry(params::SolidParticleParameters, Nᵣ)
        # Length of subdivision
        Δr = params.Rₖ/Nᵣ
        # Centre of radially divided concentric cell
        r = ([Δr*(i-0.5) for i in 1:Nᵣ])
        # Volume of subdivision (difference of outer and inner radius spheres)
        Vᵢ = 4/3*π*[(r[i] + Δr/2)^3 - (r[i] - Δr/2)^3 for i in 1:Nᵣ]
        # Inner (left) face area
        Aₗ = 4*π*[(r[i] - Δr/2)^2 for i in 1:Nᵣ]
        # Outer (right) face area
        Aᵣ = 4*π*[(r[i] + Δr/2)^2 for i in 1:Nᵣ]

        return (Nᵣ=Nᵣ, Δr=Δr, r=r, Vᵢ=Vᵢ, Aₗ=Aₗ, Aᵣ=Aᵣ)
    end
    
    # FVM parameters determined separately for each particle, optionally with different resolutions
    ne = solid_geometry(params.n, N[:Nᵣ][1])
    pe = solid_geometry(params.p, N[:Nᵣ][2])

    # Number of subdivisions per region
    Nx=N[:Nₓ]
    # Total number of subdivisions
    Nₜ=sum(Nx)
    # Lengths of regions
    Ls=[params.e.Lₙ, params.e.Lₛ, params.e.Lₚ]
    # Total length
    L=sum(Ls)
    # Define region indices sequentially (i.e., second region continues with last index +1 of first region)
    ixₙ=1:Nx[1]
    ixₛ=(Nx[1]+1):(Nx[1]+Nx[2])
    ixₚ=(Nx[1]+Nx[2]+1):(Nx[1]+Nx[2]+Nx[3])
    # Widths of subdivisions
    Δxₙ=[Ls[1]/Nx[1] for _ in ixₙ] # Negative electrode
    Δxₛ=[Ls[2]/Nx[2] for _ in ixₛ] # Separator
    Δxₚ=[Ls[3]/Nx[3] for _ in ixₚ] # Positive electrode
    # Concatenate width lists
    Δx =[Δxₙ..., Δxₛ..., Δxₚ...]
    # Cell centres
    x_centers = [
        [Δxₙ[i]*(i-0.5) for i in 1:Nx[1]]...,                 # Negative electrode
        [Δxₛ[i]*(i-0.5) + Ls[1] for i in 1:Nx[2]]...,         # Separator
        [Δxₚ[i]*(i-0.5) + Ls[1] + Ls[2] for i in 1:Nx[3]]..., # Positive electrode
    ]
    # Distance between cell centres (outer most cells have only a single neighbour)
    # Left neighbours
    Δxₗ = [nothing, [Δx[i-1]/2 + Δx[i]/2 for i in 2:Nₜ]...]
    # Right neighbours
    Δxᵣ = [[Δx[i]/2 + Δx[i+1]/2 for i in 1:Nₜ-1]..., nothing]
    # Group FVM parameters for electrolyte
    el = (Nx=Nx, Nₜ=Nₜ, Ls=Ls, L=L, ixₙ=ixₙ, ixₛ=ixₛ, ixₚ=ixₚ,
          Δx=Δx, Δxₗ=Δxₗ, Δxᵣ=Δxᵣ, x_centers=x_centers)

    # FVM parameters for ["negative particle"; "positive particle"; "electrolyte"]
    return (ne=ne, pe=pe, el=el)
end

# Trapezoidal integral (FVM definition holds)
function ∫(f, x, f₀=0)
    # Check if number of positions is equal to number of values
    # @assert length(x) == length(f)

    return sum([
            # Area of triangle between f0 (start of domain) and f1 (first cell centre)
            f[1]*x[1]/2 + f₀,
            # Area of trapezium between cell centres
            [(f[i] + f[i-1])*(x[i] - x[i-1])/2 for i in 2:length(x)]...,
    ])
end

macro integrate(f, x, f0=0)
    return :(sum([
            $f[1]*$x[1]/2 + $f0,
            [($f[i] + $f[i-1])*($x[i] - $x[i-1])/2 for i in 2:length($x)]...,
    ]))
end
