# =====================================================================================================================
# SEI.jl
# SEI growth models for SPMe
# Contains:
# 1. NoSEI: Disabled SEI growth
# 2. ReactionLimitedSEI: Reaction-limited SEI growth with Butler-Volmer kinetics
# 3. SolventDiffusionLimitedSEI: Film-diffusion-limited SEI growth with Arrhenius correction
# =====================================================================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    NoSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a zero SEI growth model (reaction disabled).

Returns a ModelingToolkit system where SEI film thickness remains static and provides ohmic resistance
based on the initial thickness. Use this when SEI growth effects are negligible or you want to exclude 
them from the simulation while maintaining the baseline resistance.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: SEI reaction parameters 
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Output Variables
- `L_sei`: SEI film thickness (static at Lf₀)
- `j_sei`: SEI current density (always 0)
- `ϕf`: Film potential drop (acts purely as a static ohmic resistor)
"""
function NoSEI(; name, p, s, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput(guess=298.15)
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)

    # Calculate initial concentration based on initial film thickness
    c_sei₀ = p.Lf₀/p.V̄*s.aₖ
    
    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = c_sei₀
        (j_sei(t))[1:N], [guess=zeros(N)]
        (ϕf(t))[1:N], [guess=zeros(N)]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t), [guess=0.0]
        ϕf_x(t), [guess=0.0]
    end
    
    eqns = [
        [j_sei[i] ~ 0 for i in 1:N]...,
        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N
    ]

    System(eqns, t; name=name, systems=[J, T, Δϕₛ, aₖ])
end


"""
    ReactionLimitedSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a reaction-limited SEI growth model.

Models SEI film formation with reaction kinetics controlled by surface overpotential.
The SEI current density follows Butler-Volmer kinetics. Use when SEI growth is fast
(high overpotential) and film diffusion resistance is negligible.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: SEI reaction kinetic parameters
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Key Parameters Used
- `p.j_sei`: Exchange current density (A/m²)
- `p.α`: Transfer coefficient (charge transfer kinetics)
- `p.U`: SEI formation potential (V vs Li/Li⁺)

# Output Variables
- `L_sei`: SEI film thickness (grows over time)
- `j_sei`: SEI current density (determined by kinetics)
- `ϕf`: Film potential (typically small)

# Physical Assumption
Reaction rate dominates over diffusion; film acts as a perfect ionic conductor.
"""
function ReactionLimitedSEI(; name, p, s, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput(guess=298.15)
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)
    
    c_sei₀ = p.Lf₀/p.V̄*s.aₖ

    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = c_sei₀
        (j_sei(t))[1:N], [guess=zeros(N)]
        (ϕf(t))[1:N], [guess=zeros(N)]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t), [guess=0.0]
        ϕf_x(t), [guess=0.0]
    end

    η_sei = [Δϕₛ.u[i] - p.U - ϕf[i] for i in 1:N]

    eqns = [
        # Exchange current density using parameter naming from src Base.jl (j_sei)
        [j_sei[i] ~ -p.j_sei*exp(min(-p.α*F/R/T.u*η_sei[i], 300)) for i in 1:N]...,

        [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N
    ]

    System(eqns, t; name=name, systems=[J, T, Δϕₛ, aₖ])
end


"""
    SolventDiffusionLimitedSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a solvent-diffusion-limited SEI growth model.

Models SEI film formation limited by solvent diffusion through the growing film.
The SEI current density decreases as the film thickens due to increasing ionic resistance
and decreasing solvent diffusion. Includes a parameter-based Arrhenius temperature correction.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: SEI reaction and film transport parameters
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Key Parameters Used
- `p.D_sol`: Solvent diffusivity in film (m²/s)
- `p.c_sol`: Solvent concentration (mol/m³)
- `p.E_sei`: Activation energy for SEI growth (J/mol)
- `p.T_ref`: Reference temperature for Arrhenius kinetics (K)
- `p.U`: SEI formation potential (V vs Li/Li⁺)
- `p.R`: Film resistivity (Ω·m)

# Output Variables
- `L_sei`: SEI film thickness (grows over time, asymptotically)
- `j_sei`: SEI current density (decreases as L_sei increases)
- `ϕf`: Film potential drop (increases with thickness)

# Physical Assumption
Film diffusion resistance and potential drop dominate; SEI growth self-limits via thickness.
"""
function SolventDiffusionLimitedSEI(; name, p, s, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput(guess=298.15)
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)

    c_sei₀ = p.Lf₀/p.V̄*s.aₖ

    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = c_sei₀
        (j_sei(t))[1:N], [guess=zeros(N)]
        (ϕf(t))[1:N], [guess=zeros(N)]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t), [guess=0.0]
        ϕf_x(t), [guess=0.0]
    end

    # All SEI growth mechanisms assumed to have Arrhenius dependence
    arrhenius = exp(min(p.E_sei / R * (1 / p.T_ref - 1 / T.u), 300))
    
    eqns = [
        # Exchange current density
        [j_sei[i] ~ -p.D_sol*p.c_sol*F/L_sei[i]*arrhenius for i in 1:N]...,

        [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N
    ]

    System(eqns, t; name=name, systems=[J, T, Δϕₛ, aₖ])
end