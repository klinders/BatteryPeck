# =====================================================================================================================
# runtests.jl
#
# Testing environment for e.g.:
# - Does the model run?
# - Do variables stay within defined limits?
# - Are parameters being correctly loaded from parameter sets?
# =====================================================================================================================

# Import package manager
using Pkg
# Import module
using BatteryPeck
# Import package
using Test

# Activate Julia environment in current directory
Pkg.activate(".")

# Start tracking time
println("Starting tests")
ti = time()

# Run test
@testset "BatteryPeck test" begin
    @test 1 == 1    # placeholder test
end

# Stop tracking time
ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")