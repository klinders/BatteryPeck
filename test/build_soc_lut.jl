# ==============================================================================
# build_soc_lut.jl
# Generates a V-SoC look-up table (SoC defined as a function of stoichiometry)
# ==============================================================================

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using CSV
using DataFrames

using Revise
using BatteryPeck

Revise.revise()

function generate_lut()
    # Define to be evaluated SoC points (e.g., 0% to 100% in 0.5% steps)
    soc_range = 0.0:0.005:1.0

    # Initialise empty DataFrame to store results
    df_lut = DataFrame(
        SoC = Float64[],
        Voltage = Float64[],
        z_n = Float64[],
        z_p = Float64[],
        c_n_init = Float64[],
        c_p_init = Float64[]
    )

    println("Starting LUT Generation, might take a few minutes due to MTK compilation...")

    for soc in soc_range
        println("\n--- Evaluating SoC = $(round(soc*100, digits=1))% ---")
        
        # Load fresh parameters for each loop
        p = Chen2020()
        
        # Map SoC to stoichiometries
        z_n_init = p.n.z_0 + soc * (p.n.z_100 - p.n.z_0)
        z_p_init = p.p.z_0 + soc * (p.p.z_100 - p.p.z_0)
        
        # Calculate concentrations
        c_n_init = z_n_init * p.n.c₊
        c_p_init = z_p_init * p.p.c₊
        
        # Assign to parameters
        p.n.c₀ = c_n_init
        p.p.c₀ = c_p_init
        
        # Disable safety voltage limits so solver doesn't abort at exactly 0% or 100%
        p.Vmin = 2.0 
        p.Vmax = 5.0

        # Build system
        @mtkbuild sys = SingleCellPack(params=p, config=(1,1))
        
        # Run 10-second rest experiment; since current is 0, voltage will be OCV.
        exp = Experiment([RestStep(10.0)])
        sol = simulate(sys, exp)
        
        # Extract resting OCV at end of experiment
        v_ocv = sol[sys.V][end]
        println("Result: Voltage = $(round(v_ocv, digits=4)) V")
        
        # Store results in DataFrame
        push!(df_lut, (soc, v_ocv, z_n_init, z_p_init, c_n_init, c_p_init))
    end

    # Save LUT to data folder
    output_file = joinpath(@__DIR__, "..", "data", "Chen2020", "soc_ocv_lut.csv")
    CSV.write(output_file, df_lut)

    println("\n LUT successfully generated and saved to: $output_file")
end

# Execute function
generate_lut()