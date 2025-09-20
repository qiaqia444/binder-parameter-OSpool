#!/usr/bin/env julia

using Pkg
Pkg.activate(".")
using Printf

include("src/BinderSimForced.jl")
using .BinderSimForced

println("Edwards-Anderson Binder Parameter - Forced +1 Measurement Dynamics")
println("="^70)

# Test parameters
L_values = [12, 16]
lambda_values = [0.1, 0.3, 0.5, 0.7, 0.9]

println("Testing forced +1 measurement dynamics...")
println("Expected behavior: Since all measurements are forced to +1,")
println("the dynamics should be deterministic and different from random case.")
println()

results = []

for L in L_values
    println("\\n" * "="^50)
    println("SYSTEM SIZE L = $L")
    println("="^50)
    
    for lambda in lambda_values
        println("\\nTesting λ = $lambda:")
        println("-"^30)
        
        try
            result = ea_binder_mc_forced(L; 
                                       lambda_x=lambda, 
                                       lambda_zz=lambda, 
                                       ntrials=2,  # Run 2 trials to verify determinism
                                       maxdim=128, 
                                       cutoff=1e-10)
            
            push!(results, (L=L, lambda=lambda, B=result.B, 
                          B_std=result.B_std_of_trials, success=true))
            
            println("SUCCESS: B = $(round(result.B, digits=4))")
            println("Std between trials = $(round(result.B_std_of_trials, digits=6))")
            
            if result.B_std_of_trials < 1e-12
                println("✓ Deterministic: trials are identical")
            else
                println("⚠ Non-deterministic: trials differ")
            end
            
        catch e
            println("FAILED: $e")
            push!(results, (L=L, lambda=lambda, B=NaN, B_std=NaN, success=false))
        end
    end
end

println("\\n" * "="^70)
println("SUMMARY OF FORCED +1 MEASUREMENT RESULTS")
println("="^70)

println("L\\t\\tλ\\t\\tBinder B\\t\\tStd\\t\\tStatus")
println("-"^70)

for result in results
    status = result.success ? "✓" : "✗"
    B_str = result.success ? @sprintf("%.4f", result.B) : "FAIL"
    std_str = result.success ? @sprintf("%.2e", result.B_std) : "---"
    
    println("$(result.L)\\t\\t$(result.lambda)\\t\\t$B_str\\t\\t$std_str\\t\\t$status")
end

println("\\n" * "="^70)
println("COMPARISON NOTES:")
println("• Forced +1 measurements create deterministic trajectories")
println("• Results should be reproducible (std ≈ 0)")
println("• Physics may differ significantly from random measurement case")
println("• This represents a specific measurement protocol study")
println("="^70)