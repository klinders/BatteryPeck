module ParticleCracking

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    NoCracking(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a model without cracking (reaction disabled).

Returns a ModelingToolkit system where crack thickness and SEI layer remains zero and provides no ohmic resistance.
Use this when cracking effects are negligible or you want to exclude them from the simulation.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: Cracking reaction parameters (unused in this model)
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Output Variables
- `L_sei`: SEI film thickness (always 0)
- `j_sei`: SEI current density (always 0)
- `ϕf`: Film potential (always 0)
"""
function NoMechanics(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        # Mechanics
        (a_cr(t))[1:N] = 0

        # SEI concentration
        (c_sei(t))[1:N] = 0
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N] = 0

        a_cr_x(t)
        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)

    end

  
    eqns = [
        [Dt(a_cr[i]) ~ 0 for i in 1:N]...,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ 0 for i in 1:N]...,

        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        a_cr_x ~ sum([a_cr[i] for i in 1:N])/N,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        Q_sei ~ 0,

    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
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
- `p.j_sei₀`: Exchange current density (A/m²)
- `p.α`: Transfer coefficient (charge transfer kinetics)
- `p.U`: SEI formation potential (V vs Li/Li⁺)

# Output Variables
- `L_sei`: SEI film thickness (grows over time)
- `j_sei`: SEI current density (determined by kinetics)
- `ϕf`: Film potential (typically small)

# Physical Assumption
Reaction rate dominates over diffusion; film acts as perfect ionic conductor.
"""
function SwellingOnly(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named c_s_r = RealInput()
    @named c_s_surf = RealInput()


    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        l_cr(t) = 0 
        r_surf(t)
        a_cr(t)
        σₜ(t)
        u_d(t)
        
        # SEI concentration
        (c_sei(t))[1:N] = 0#scale
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        l_cr_x(t)
        r_surf_x(t)
        a_cr_x(t)
        σₜ_x(t)
        u_d_x(t)

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_sei(t)
    end

    # Physical constants for stress/displacement calculations (Ai2019)
    stress_geometric_factor = 3.0
    displacement_geometric_factor = 3.0
    crack_roughness_factor = 2.0

    ρ_cr = 3.18e-15
    w_cr = 1.5e-8
    Ω = 3.1e-06 # pos: 1.25e-05 (parial molar volume)
    E = 1.5e10 # pos: 3.75e11 (Youngs modulus)
    ν = 0.3 # pos: 0.2 (poissons ratio)
    R0 = s.Rₖ
    c₀ = 0

    eqns = [

        # Cracking
        l_cr ~ 0,
        r_surf ~ 1 + crack_roughness_factor*l_cr*ρ_cr*w_cr,
        a_cr ~ (r_surf - 1)*aₖ.u,
        σₜ ~ Ω*E*(c_s_r.u - c_s_surf.u)/stress_geometric_factor/(1.0 - ν),
        u_d ~ Ω*R0*(c_s_r.u - c₀)/displacement_geometric_factor,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        # [j_sei[i] ~ -p.j_sei₀*exp(-p.α*F/R/T.u*η_sei[i]) for i in 1:N]...,

        # [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        # [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        # [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,

        l_cr_x ~ l_cr,
        r_surf_x ~ r_surf,
        a_cr_x ~ a_cr,
        σₜ_x ~ σₜ,
        u_d_x ~ u_d,

        # L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        # c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        # j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        # ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        # Q_sei ~ (c_sei_x-c_sei₀)*p.V̄*p.z*F/3600,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ, c_s_r, c_s_surf])
end

"""
    SolventDiffusionLimitedSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a solvent-diffusion-limited SEI growth model.

Models SEI film formation limited by solvent diffusion through the growing film.
The SEI current density decreases as the film thickens due to increasing ionic resistance
and decreasing solvent diffusion. Use when film resistance dominates over reaction kinetics.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: SEI reaction and film transport parameters
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Key Parameters Used
- `p.D_sol`: Solvent diffusivity in film (m²/s)
- `p.c_sol`: Solvent concentration (mol/m³)
- `p.U`: SEI formation potential (V vs Li/Li⁺)
- `p.R`: Film resistivity (Ω·m)

# Output Variables
- `L_sei`: SEI film thickness (grows over time, asymptotically)
- `j_sei`: SEI current density (decreases as L_sei increases)
- `ϕf`: Film potential drop (increases with thickness)

# Physical Assumption
Film diffusion resistance and potential drop dominate; SEI growth self-limits via thickness.
"""
function SwellingAndCracking(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named c_s_r = RealInput()
    @named c_s_surf = RealInput()


    # Time derivative operator
    Dt = Differential(t)

    l_cr_0 = 2e-8
    
    @variables begin
        l_cr(t) = l_cr_0 
        r_surf(t)
        a_cr(t)
        σₜ(t)
        u_d(t)
        
        # SEI concentration
        (c_sei(t))[1:N] = 0#scale
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        l_cr_x(t)
        r_surf_x(t)
        a_cr_x(t)
        σₜ_x(t)
        u_d_x(t)

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_sei(t)
    end

    # Physical constants for stress/displacement calculations (Ai2019)
    stress_geometric_factor = 3.0
    displacement_geometric_factor = 3.0
    crack_roughness_factor = 2.0

    ρ_cr = 3.18e-15
    w_cr = 1.5e-8
    Ω = 3.1e-06 # pos: 1.25e-05 (parial molar volume)
    E = 1.5e10 # pos: 3.75e11 (Youngs modulus)
    ν = 0.3 # pos: 0.2 (poissons ratio)
    R0 = s.Rₖ
    c₀ = 0

    k_cr = 3.9e-20
    b_cr = 1.12
    m_cr = 2.2


    dK_SIF = ifelse(σₜ >= 0, σₜ*b_cr * sqrt(pi*l_cr), 0)

    eqns = [

        # Cracking
        Dt(l_cr) ~ k_cr*(dK_SIF^m_cr)/3600,
        r_surf ~ 1 + crack_roughness_factor*l_cr*ρ_cr*w_cr,
        a_cr ~ (r_surf - 1)*aₖ.u,
        σₜ ~ Ω*E*(c_s_r.u - c_s_surf.u)/stress_geometric_factor/(1.0 - ν),
        u_d ~ Ω*R0*(c_s_r.u - c₀)/displacement_geometric_factor,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        # [j_sei[i] ~ -p.j_sei₀*exp(-p.α*F/R/T.u*η_sei[i]) for i in 1:N]...,

        # [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        # [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        # [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,

        l_cr_x ~ l_cr,
        r_surf_x ~ r_surf,
        a_cr_x ~ a_cr,
        σₜ_x ~ σₜ,
        u_d_x ~ u_d,

        # L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        # c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        # j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        # ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        # Q_sei ~ (c_sei_x-c_sei₀)*p.V̄*p.z*F/3600,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ, c_s_r, c_s_surf])
end

end