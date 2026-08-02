# ==============================================================================
# build_dcir_lut.jl
# Generates a DCIR vs. SoC Lookup Table via Algebraic Pulse Testing
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using BatteryToolkit
using CSV, DataFrames
using Plots
using Logging

function get_sym(sys, param_name)
    for p in parameters(sys)
        if contains(string(p), param_name) return p end
    end
    error("Parameter '$param_name' not found")
end

function build_dcir_lut()
    println("=====================================================")
    println("--- Generating Dynamic DCIR vs. SoC Lookup Table ---")
    println("=====================================================")
    
    soc_targets = collect(1.0:-0.05:0.0)
    dcir_values = Float64[]
    pulse_current = -5.0 
    
    println("\nRunning algebraic pulse tests (Recompiling at each SoC state)...")
    
    for soc in soc_targets
        p = Chen2020()
        p.Vmin = 1.0; p.Vmax = 5.5
        p.n.c₀ = (p.n.z_0 + soc * (p.n.z_100 - p.n.z_0)) * p.n.c₊  
        p.p.c₀ = (p.p.z_0 + soc * (p.p.z_100 - p.p.z_0)) * p.p.c₊
        
        @named cell = SPMe(params=p, side_reactions=true)
        @named load = Current(); @named ground = Ground()
        @parameters t I_app=0.0
        
        eqs = [
            cell.T.u ~ 298.15,
            load.I.u ~ I_app,
            connect(load.p, cell.p),
            connect(cell.n, load.n, ground.g)
        ]
        
        @named sys = ODESystem(eqs, t, [], [I_app]; systems=[cell, load, ground])
        
        # Suppress compilation warnings to keep terminal clean
        sys_simp = with_logger(ConsoleLogger(stderr, Logging.Error)) do
            structural_simplify(sys)
        end
        
        I_sym = get_sym(sys_simp, "I_app")
        
        # Get instantaneous OCV (I = 0.0A)
        prob_ocv = ODEProblem(sys_simp, [I_sym => 0.0], (0.0, 1.0); sparse=true, jac=false)
        V_ocv = init(prob_ocv, QNDF(autodiff=true))[sys_simp.cell.v]
        
        # Get instantaneous loaded voltage (I = -5.0A)
        prob_pulse = ODEProblem(sys_simp, [I_sym => pulse_current], (0.0, 1.0); sparse=true, jac=false)
        V_loaded = init(prob_pulse, QNDF(autodiff=true))[sys_simp.cell.v]
        
        # Calculate exact resistance
        R_dcir = abs((V_ocv - V_loaded) / pulse_current)
        push!(dcir_values, R_dcir)
        
        println("  SoC: $(lpad(round(Int, soc*100), 3))% | V_ocv: $(round(V_ocv,digits=3)) V | DCIR: $(round(R_dcir*1000,digits=2)) mΩ")
    end
    
    # Save LUT
    out_path = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_dcir_lut.csv")
    df = DataFrame(SoC = soc_targets, DCIR_Ohms = dcir_values)
    CSV.write(out_path, df)
    
    println("\n[!] Success! Saved DCIR LUT to:")
    println("    $out_path")
    
    # Plot
    p1 = plot(soc_targets .* 100, dcir_values .* 1000, 
              title="Instantaneous DCIR vs. State of Charge",
              xlabel="State of Charge (%)", ylabel="Internal Resistance (mΩ)",
              lw=3, color=:purple, marker=:circle, grid=true, legend=false, xflip=true)
              
    display(p1)
end

Base.invokelatest(build_dcir_lut)