module ParticleCracking

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    NoMechanics(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a model without cracking (reaction disabled).

Returns a ModelingToolkit system where crack thickness and SEI layer remains zero and provides no ohmic resistance.
Use this when cracking effects are negligible or you want to exclude them from the simulation.
"""
function NoMechanics(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    @parameters t
    N = g.el.Nx[1]
    
    @named J = RealInput()
    @named T = RealInput(guess=298.15)
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    Dt = Differential(t)
    
    @variables begin
        (a_cr(t))[1:N] = 0.0
        (c_sei(t))[1:N] = 0.0
        (j_sei(t))[1:N], [guess=zeros(N)]
        (ϕf(t))[1:N], [guess=zeros(N)]
        (L_sei(t))[1:N] = 0.0

        a_cr_x(t), [guess=0.0]
        c_sei_x(t), [guess=0.0]
        L_sei_x(t), [guess=0.0]
        j_sei_x(t), [guess=0.0]
        ϕf_x(t), [guess=0.0]
        Q_sei(t), [guess=0.0]
    end

    eqns = [
        [Dt(a_cr[i]) ~ 0 for i in 1:N]...,
        [j_sei[i] ~ 0 for i in 1:N]...,
        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,
        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        
        a_cr_x ~ sum([a_cr[i] for i in 1:N])/N,
        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        Q_sei ~ 0.0
    ]
    System(eqns, t; name=name, systems=[J, T, Δϕₛ, aₖ])
end

"""
    SwellingOnly(; name, p::SideReactionParameters, s::SolidParticleParameters, V, g)

Models radial and tangential stresses within the active particle caused by lithium 
intercalation (swelling), but disables active crack propagation and SEI deposition 
on newly exposed surfaces.

Use this for positive electrodes (e.g., NMC) where stress drives LAM but 
SEI formation inside cracks is considered negligible.
"""
function SwellingOnly(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    @parameters t
    R = 8.314 
    F = 96485 
    N = g.el.Nx[1]

    stress_geometric_factor = 3.0
    displacement_geometric_factor = 3.0
    R0 = s.Rₖ
    c₀ = 0.0

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named c_s_r = RealInput()
    @named c_s_surf = RealInput()

    Dt = Differential(t)
    
    @variables begin
        l_cr(t), [guess=0.0]
        r_surf(t), [guess=1.0]
        a_cr(t), [guess=0.0]
        σₜ(t), [guess=0.0]
        σᵣ(t), [guess=0.0]
        u_d(t), [guess=0.0]
        
        (c_sei(t))[1:N] = 0.0
        (j_sei(t))[1:N], [guess=zeros(N)]
        (ϕf(t))[1:N], [guess=zeros(N)]
        (L_sei(t))[1:N], [guess=zeros(N)]

        l_cr_x(t), [guess=0.0]
        r_surf_x(t), [guess=1.0]
        a_cr_x(t), [guess=0.0]
        σₜ_x(t), [guess=0.0]
        σᵣ_x(t), [guess=0.0]
        u_d_x(t), [guess=0.0]

        c_sei_x(t), [guess=0.0]
        L_sei_x(t), [guess=0.0]
        j_sei_x(t), [guess=0.0]
        ϕf_x(t), [guess=0.0]
        Q_sei(t), [guess=0.0]
    end

    eqns = [
        l_cr ~ 0,
        r_surf ~ 1,
        a_cr ~ 0,
        σₜ ~ s.Ω*s.E*(c_s_r.u - c_s_surf.u)/stress_geometric_factor/(1.0 - s.ν),
        σᵣ ~ 0,
        u_d ~ s.Ω*R0*(c_s_r.u - c₀)/displacement_geometric_factor,

        [j_sei[i] ~ 0 for i in 1:N]...,
        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ 0 for i in 1:N]...,
        [ϕf[i] ~ 0 for i in 1:N]...,

        l_cr_x ~ l_cr,
        r_surf_x ~ r_surf,
        a_cr_x ~ a_cr,
        σₜ_x ~ σₜ,
        σᵣ_x ~ σᵣ,
        u_d_x ~ u_d,

        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        Q_sei ~ 0.0
    ]
    System(eqns, t; name=name, systems=[J, T, Δϕₛ, aₖ, c_s_r, c_s_surf])
end

"""
    SwellingAndCracking(; name, p::SideReactionParameters, s::SolidParticleParameters, V, g)

Models comprehensive mechanical degradation through particle swelling and crack propagation.
Calculates crack length and the associated increase in electrochemical surface area.

Solves secondary SEI deposition inside the newly formed cracks using an EC-diffusion
limited approach, coupling structural damage directly to cyclable lithium loss.
"""
function SwellingAndCracking(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    @parameters t
    R = 8.314 
    F = 96485 
    N = g.el.Nx[1]

    stress_geometric_factor = 3.0
    displacement_geometric_factor = 3.0
    crack_roughness_factor = 2.0
    R0 = s.Rₖ
    c₀ = 0.0
    k_cr = 3.9e-20
    b_cr = 1.12
    m_cr = 2.2

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named c_s_r = RealInput()
    @named c_s_surf = RealInput()

    Dt = Differential(t)

    l_cr_0 = 2e-8
    L₀ = 5e-13
    c_sei₀ = L₀/p.V̄*s.aₖ
    c_ec_0 = 4541.0
    D_ec = 2e-18
    k_sei = 1e-12
    
    arrhenius = exp(p.E_sei / R * (1 / p.T_ref - 1 / T.u))
    
    @variables begin
        l_cr(t) = l_cr_0 
        r_surf(t), [guess=1.0]
        a_cr(t), [guess=0.0]
        σₜ(t), [guess=0.0]
        σᵣ(t), [guess=0.0]
        u_d(t), [guess=0.0]
        
        (c_ec(t))[1:N], [guess=ones(N)*c_ec_0]
        (c_sei(t))[1:N] = c_sei₀
        (j_sei(t))[1:N], [guess=zeros(N)]
        (ϕf(t))[1:N], [guess=zeros(N)]
        (L_sei(t))[1:N], [guess=ones(N)*L₀]

        l_cr_x(t), [guess=l_cr_0]
        r_surf_x(t), [guess=1.0]
        a_cr_x(t), [guess=0.0]
        σₜ_x(t), [guess=0.0]
        σᵣ_x(t), [guess=0.0]
        u_d_x(t), [guess=0.0]

        c_ec_x(t), [guess=c_ec_0]
        c_sei_x(t), [guess=c_sei₀]
        L_sei_x(t), [guess=L₀]
        j_sei_x(t), [guess=0.0]
        ϕf_x(t), [guess=0.0]
        Q_sei(t), [guess=0.0]
    end

    dK_SIF = ifelse(σₜ >= 0, σₜ*b_cr* sqrt(pi*l_cr), 0.0)
    η_sei = [Δϕₛ.u[i] - p.U - ϕf[i] for i in 1:N]
    k_exp = [k_sei*exp(-p.α*F/R/T.u*η_sei[i]) for i in 1:N]
    L_over_D = [L_sei[i]/D_ec for i in 1:N]

    eqns = [
        Dt(l_cr) ~ k_cr*(dK_SIF^m_cr)/3600,
        r_surf ~ 1 + crack_roughness_factor*l_cr*s.ρ_cr*s.w_cr,
        a_cr ~ (r_surf - 1)*aₖ.u,
        σₜ ~ s.Ω*s.E*(c_s_r.u - c_s_surf.u)/stress_geometric_factor/(1.0 - s.ν),
        σᵣ ~ 0,
        u_d ~ s.Ω*R0*(c_s_r.u - c₀)/displacement_geometric_factor,

        [j_sei[i] ~ -F*c_ec_0*k_exp[i]/(1 + L_over_D[i]*k_exp[i])*arrhenius for i in 1:N]...,
        [c_ec[i] ~ c_ec_0/(1 + L_over_D[i]*k_exp[i]) for i in 1:N]...,
        [Dt(c_sei[i]) ~ -a_cr*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/max(a_cr, 1e-10) for i in 1:N]...,
        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,

        l_cr_x ~ l_cr,
        r_surf_x ~ r_surf,
        a_cr_x ~ a_cr,
        σₜ_x ~ σₜ,
        σᵣ_x ~ σᵣ,
        u_d_x ~ u_d,

        L_sei_x ~ sum([L_sei[i] for i in 1:N])/N,
        c_ec_x ~ sum([c_ec[i] for i in 1:N])/N,
        c_sei_x ~ sum([c_sei[i] for i in 1:N])/N,
        j_sei_x ~ sum([j_sei[i] for i in 1:N])/N,
        ϕf_x ~ sum([ϕf[i] for i in 1:N])/N,
        
        Q_sei ~ (c_sei_x-c_sei₀)*V*p.z*F/3600
    ]
    System(eqns, t; name=name, systems=[J, T, Δϕₛ, aₖ, c_s_r, c_s_surf])
end
end