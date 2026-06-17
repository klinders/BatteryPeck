# ==============================================================================
# build_soc_lut.jl
# Generates voltage state of charge look up table
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using CSV
using DataFrames

using Revise
using BatteryToolkit

Revise.revise()

function generate_lut()
    # Define state of charge evaluation range
    soc_range = 0.0:0.005:1.0

    # Initialise empty dataframe
    df_lut = DataFrame(
        SoC = Float64[],
        Voltage = Float64[],
        z_n = Float64[],
        z_p = Float64[],
        c_n_init = Float64[],
        c_p_init = Float64[]
    )

    println("Building system and compiling ordinary differential equation problem once...")
    
    # Load parameters
    p = Chen2020()
    
    # Disable safety voltage limits
    p.Vmin = 2.0 
    p.Vmax = 5.0

    # Build mathematical system
    @mtkbuild sys = SingleCellPack(params=p, config=(1,1))
    
    # Define rest experiment
    exp = Experiment([RestStep(10.0)])
    
    # Compile base problem using zero inputs
    prob = ODEProblem(sys, [sys.Pin => 0.0, sys.Iin => 0.0], (0.0, exp.tend))

    for soc in soc_range
        println("\n--- Evaluating SoC = $(round(soc*100, digits=1))% ---")
        
        # Map state of charge to stoichiometries
        z_n_init = p.n.z_0 + soc * (p.n.z_100 - p.n.z_0)
        z_p_init = p.p.z_0 + soc * (p.p.z_100 - p.p.z_0)
        
        # Calculate solid particle concentrations
        c_n_init = z_n_init * p.n.c₊
        c_p_init = z_p_init * p.p.c₊
        
        # Create dictionary mapping new initial conditions
        u0_map = Dict()
        for i in 1:length(sys.cell.ne.c)
            u0_map[sys.cell.ne.c[i]] = c_n_init
            u0_map[sys.cell.pe.c[i]] = c_p_init
        end
        
        # Remake problem with updated initial conditions
        prob_new = remake(prob, u0=u0_map)
        
        # Run rest experiment using fast simulate method
        sol = simulate(sys, prob_new, exp)
        
        # Extract open circuit voltage
        v_ocv = sol[sys.V][end]
        println("Result: Voltage = $(round(v_ocv, digits=4)) V")
        
        # Store results
        push!(df_lut, (soc, v_ocv, z_n_init, z_p_init, c_n_init, c_p_init))
    end

    # Save look up table
    output_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
    CSV.write(output_file, df_lut)

    println("\nLook up table successfully generated and saved to: $output_file")
end

# Execute function
generate_lut()