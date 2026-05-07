#!/usr/bin/env julia

"""
Minimal test - just 1 trial to check if it completes at all
"""

using Random
using Statistics
using ITensors, ITensorMPS

# Load consolidated Rényi-2 Binder dynamics module
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_dynamics_density_matrix_1.jl")

function main()
    println("="^70)
    println("MINIMAL TEST: 1 trial only")
    println("="^70)
    println()
    
    L = 8
    lambda_x = 0.49
    lambda_zz = 0.21
    P_x = 0.3
    P_zz = 0.3
    ntrials = 1  # Just 1 trial
    maxdim = 256
    cutoff = 1e-12
    seed = 12345
    
    println("Configuration: L=$L, λ_x=$lambda_x, λ_zz=$lambda_zz, P_x=$P_x, ntrials=$ntrials")
    println()
    
    println("Starting single trial...")
    t_start = time()
    
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
        verbose = true
    )
    
    t_elapsed = time() - t_start
    
    println()
    println("✓ Completed in $(round(t_elapsed, digits=1)) seconds")
    println("  B = $(round(result.B, digits=4))")
    println("  Purity = $(round(result.purity_bar, digits=4))")
    println("  Valid trials = $(result.n_valid)/$(result.ntrials)")
end

main()
