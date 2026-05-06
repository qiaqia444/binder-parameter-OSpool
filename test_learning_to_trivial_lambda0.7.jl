#!/usr/bin/env julia

"""
Test script for Learning-to-Trivial Transition Scan (λ_x = 0.49, λ_zz = 0.21)

Quick verification that the setup works before submitting 2,200 jobs to the cluster.

Parameters:
- Small system: L = 8
- Single point: P_x = 0.3 (middle of phase transition region)
- Few trials: 50 (fast test)
- Output: test_learning_to_trivial_lambda0.7_results/
"""

using JSON
using Random
using Statistics
using Dates
using ITensors, ITensorMPS

# Load the CORRECTED separate pipeline (with channel averaging dephasing fix)
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_binder.jl")
include("src_new/renyi2_dynamics_separate.jl")

function main()
    println("="^70)
    println("TEST: Learning-to-Trivial Transition Scan (λ_x = 0.49, λ_zz = 0.21)")
    println("="^70)
    println()
    
    # Fixed test parameters
    L = 8
    lambda_x = 0.49
    lambda_zz = 0.21
    P_x = 0.3
    P_zz = 0.3
    ntrials = 50  # Quick test
    maxdim = 256
    cutoff = 1e-12
    seed = 12345
    output_dir = "test_learning_to_trivial_lambda0.7_results"
    
    println("Test Configuration:")
    println("  L = $L")
    println("  λ_x = $lambda_x (fixed)")
    println("  λ_zz = $lambda_zz (fixed)")
    println("  P_x = $P_x")
    println("  P_zz = $P_zz")
    println("  ntrials = $ntrials")
    println("  Output: $output_dir")
    println()
    
    # Create output directory
    mkpath(output_dir)
    
    println("Starting Rényi-2 Binder calculation (using separate pipeline with corrected dephasing)...")
    t_start = time()
    
    # Run calculation using CORRECTED separate pipeline with channel averaging dephasing
    result = renyi2_binder_density_matrix_separate(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        maxdim = maxdim,
        cutoff = cutoff,
        seed = seed,
        verbose = false
    )
    
    t_elapsed = time() - t_start
    
    println("✓ Calculation completed in $(round(t_elapsed, digits=1)) seconds")
    println()
    
    # Print results
    println("Results:")
    println("  B (ensemble) = $(round(result.B, digits=4))")
    println("  B (mean of trials) = $(round(result.B_mean_of_trials, digits=4)) ± $(round(result.B_std_of_trials, digits=4))")
    println("  M₂_bar = $(round(result.M2_bar, digits=6))")
    println("  M₄_bar = $(round(result.M4_bar, digits=6))")
    println("  Purity = $(round(result.purity_bar, digits=4))")
    println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
    println()
    
    # Save results
    result_dict = Dict(
        "test" => true,
        "L" => L,
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "P_x" => P_x,
        "P_zz" => P_zz,
        "B" => result.B,
        "B_mean_of_trials" => result.B_mean_of_trials,
        "B_std_of_trials" => result.B_std_of_trials,
        "M2_bar" => result.M2_bar,
        "M4_bar" => result.M4_bar,
        "purity_bar" => result.purity_bar,
        "ntrials" => result.ntrials,
        "n_valid" => result.n_valid,
        "n_invalid" => result.n_invalid,
        "maxdim" => maxdim,
        "cutoff" => cutoff,
        "time_seconds" => t_elapsed,
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    )
    
    outpath = joinpath(output_dir, "test_result.json")
    open(outpath, "w") do f
        JSON.print(f, result_dict, 2)
    end
    
    println("✓ Results saved to: $outpath")
    println()
    
    # Validation checks
    println("Validation Checks:")
    all_pass = true
    
    # Check 1: Purity should be < 1.0 due to dephasing
    if result.purity_bar < 1.0
        println("  ✓ Purity < 1.0 (dephasing is working): $(round(result.purity_bar, digits=4))")
    else
        println("  ✗ WARNING: Purity = 1.0 (dephasing may not be applied)")
        all_pass = false
    end
    
    # Check 2: Binder should be positive and < 1
    if 0 <= result.B <= 1.0
        println("  ✓ Binder in valid range [0, 1]: $(round(result.B, digits=4))")
    else
        println("  ✗ WARNING: Binder out of valid range: $(result.B)")
        all_pass = false
    end
    
    # Check 3: Most trials should be valid
    valid_fraction = result.n_valid / result.ntrials
    if valid_fraction > 0.8
        println("  ✓ Valid trials: $(result.n_valid)/$(result.ntrials) ($(round(100*valid_fraction, digits=1))%)")
    else
        println("  ✗ WARNING: Many invalid trials: $(result.n_valid)/$(result.ntrials)")
        all_pass = false
    end
    
    println()
    if all_pass
        println("="^70)
        println("✓ TEST PASSED - Setup is working correctly!")
        println("="^70)
        println("You can now submit the full batch with:")
        println("  ./submit_learning_to_trivial_lambda0.7.sh")
    else
        println("="^70)
        println("✗ TEST WARNINGS - Check issues above before submitting batch")
        println("="^70)
    end
    println()
end

main()
