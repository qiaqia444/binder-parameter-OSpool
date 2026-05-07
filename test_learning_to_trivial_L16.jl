#!/usr/bin/env julia

"""
Test Learning-to-Trivial with L=16 (larger system)

Configuration:
  L = 16
  λ_x = 0.49 (fixed)
  λ_zz = 0.21 (fixed)
  P_x = 0.3
  P_zz = 0.3
  ntrials = 50
"""

using JSON
using Random
using Statistics
using Dates
using ITensors, ITensorMPS

# Load consolidated Rényi-2 Binder dynamics module
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_dynamics_density_matrix_1.jl")

function main()
    println("="^70)
    println("TEST: Learning-to-Trivial L=16 (larger system)")
    println("="^70)
    println()
    
    # Test parameters
    L = 16
    lambda_x = 0.49
    lambda_zz = 0.21
    P_x = 0.3
    P_zz = 0.3
    ntrials = 50
    maxdim = 256
    cutoff = 1e-12
    seed = 12345
    
    println("Test Configuration:")
    println("  L = $L (larger system)")
    println("  λ_x = $lambda_x (fixed)")
    println("  λ_zz = $lambda_zz (fixed)")
    println("  P_x = $P_x")
    println("  P_zz = $P_zz")
    println("  ntrials = $ntrials")
    println()
    
    println("Starting Rényi-2 Binder calculation...")
    t_start = time()
    
    # Run calculation using CORRECTED separate pipeline
    result = renyi2_binder_density_matrix_dynamics(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        maxdim = maxdim,
        cutoff = cutoff,
        seed = seed,
        verbose = false,
    )
    
    t_elapsed = time() - t_start
    
    println("✓ Calculation completed in $(round(t_elapsed, digits=1)) seconds")
    println()
    println("Results:")
    println("  B (ensemble) = $(round(result.B, digits=6))")
    println("  B (mean ± std) = $(round(result.B_mean_of_trials, digits=6)) ± $(round(result.B_std_of_trials, digits=6))")
    println("  M₂_bar = $(round(result.M2_bar, digits=8))")
    println("  M₄_bar = $(round(result.M4_bar, digits=8))")
    println("  Purity = $(round(result.purity_bar, digits=6))")
    println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
    println()
    
    # Validation checks
    println("Validation Checks:")
    
    # Check 1: Purity should be significantly reduced (dephasing present)
    if result.purity_bar < 0.5
        println("  ✓ Purity significantly reduced by dephasing: $(round(result.purity_bar, digits=6))")
    else
        println("  ⚠ Purity not strongly reduced: $(round(result.purity_bar, digits=6))")
    end
    
    # Check 2: Binder should be in valid range
    if 0.0 <= result.B <= 1.0
        println("  ✓ Binder in valid range [0, 1]: $(round(result.B, digits=6))")
    else
        println("  ✗ ERROR: Binder outside [0, 1]: $(round(result.B, digits=6))")
    end
    
    # Check 3: All trials should be valid
    if result.n_valid == result.ntrials
        println("  ✓ All trials valid: $(result.n_valid)/$(result.ntrials) (100.0%)")
    else
        println("  ✗ WARN: Some invalid trials: $(result.n_valid)/$(result.ntrials) ($(round(100.0*result.n_valid/result.ntrials, digits=1))%)")
    end
    
    println()
    
    # Save results
    output_dir = "test_learning_to_trivial_L16_results"
    mkpath(output_dir)
    
    output_file = joinpath(output_dir, "test_L16_Px0.3.json")
    open(output_file, "w") do f
        output_dict = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "P_x" => P_x,
            "P_zz" => P_zz,
            "ntrials" => ntrials,
            "B" => result.B,
            "B_mean" => result.B_mean_of_trials,
            "B_std" => result.B_std_of_trials,
            "M2_bar" => result.M2_bar,
            "M4_bar" => result.M4_bar,
            "purity_bar" => result.purity_bar,
            "valid_trials" => result.n_valid,
            "timestamp" => Dates.now(),
        )
        JSON.print(f, output_dict)
    end
    
    println("✓ Results saved to: $output_file")
    println()
    
    if 0.0 <= result.B <= 1.0 && result.n_valid == result.ntrials
        println("="^70)
        println("✓ TEST PASSED - L=16 case works correctly!")
        println("="^70)
        return 0
    else
        println("="^70)
        println("✗ TEST FAILED - Issues detected")
        println("="^70)
        return 1
    end
end

exit(main())
