#!/usr/bin/env julia
"""
Test script: learning_to_trivial_lambda0.7 with SEPARATE Rényi-2 pipeline

Parameters:
  λ_x = 0.49
  λ_zz = 0.21
  L = 8, 10
  P_x = P_zz ∈ {0.1, 0.3}
  ntrials = 50
"""

using Random, Statistics, ITensors, ITensorMPS, JSON, Dates

include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_dynamics_density_matrix_1.jl")

function test_learning_to_trivial_lambda0_7_separate()
    lambda_x = 0.49
    lambda_zz = 0.21
    ntrials = 50
    maxdim = 256
    cutoff = 1e-12
    
    test_configs = [
        (L=8, P=0.1, seed=1001),
        (L=8, P=0.3, seed=1003),
        (L=10, P=0.1, seed=2001),
        (L=10, P=0.3, seed=2003),
    ]
    
    results_list = []
    
    println("\n" * "="^80)
    println("TESTING: learning_to_trivial_lambda0.7 with SEPARATE Rényi-2 Pipeline")
    println("="^80)
    println("Parameters: λ_x=$lambda_x, λ_zz=$lambda_zz, ntrials=$ntrials")
    println("="^80 * "\n")
    
    for config in test_configs
        L = config.L
        P = config.P
        seed = config.seed
        
        print("L=$L, P=$P: ")
        flush(stdout)
        
        result = renyi2_binder_density_matrix_dynamics(
            L;
            lambda_x=lambda_x,
            lambda_zz=lambda_zz,
            P_x=P,
            P_zz=P,
            ntrials=ntrials,
            maxdim=maxdim,
            cutoff=cutoff,
            seed=seed,
            strobe=:after_full_layer,
            verbose=false,
        )
        
        # Check validity
        valid_pct = 100.0 * result.n_valid / result.ntrials
        
        println("B=$(round(result.B, digits=4)), M2=$(round(result.M2_bar, digits=6)), Purity=$(round(result.purity_bar, digits=4)), Valid=$valid_pct%")
        
        # Store result
        push!(results_list, Dict(
            :L => L,
            :P => P,
            :lambda_x => lambda_x,
            :lambda_zz => lambda_zz,
            :B => result.B,
            :B_std => result.B_std_of_trials,
            :M2 => result.M2_bar,
            :M4 => result.M4_bar,
            :purity => result.purity_bar,
            :n_valid => result.n_valid,
            :n_trials => result.ntrials,
            :timestamp => now(),
        ))
    end
    
    println("\n" * "="^80)
    println("SUMMARY")
    println("="^80)
    
    for res in results_list
        println("L=$(res[:L]), P=$(res[:P]): B=$(round(res[:B], digits=4)), Purity=$(round(res[:purity], digits=4))")
    end
    
    println("\n✓ All tests completed successfully!")
    
    return results_list
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_learning_to_trivial_lambda0_7_separate()
end
