#!/usr/bin/env julia
"""
Simple test: learning_to_trivial_lambda0.7 with separate Rényi-2 pipeline
"""

using Random, Statistics, ITensors, ITensorMPS

include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_dynamics_density_matrix_1.jl")

lambda_x = 0.49
lambda_zz = 0.21

println("Testing learning_to_trivial_lambda0.7 with separate pipeline")
println("="^60)

# Test case 1
result1 = renyi2_binder_density_matrix_separate(
    8; lambda_x=lambda_x, lambda_zz=lambda_zz, P_x=0.1, P_zz=0.1,
    ntrials=20, seed=42, verbose=false
)

println("L=8, P=0.1: B=$(round(result1.B, digits=4)), Purity=$(round(result1.purity_bar, digits=4))")

# Test case 2
result2 = renyi2_binder_density_matrix_separate(
    8; lambda_x=lambda_x, lambda_zz=lambda_zz, P_x=0.3, P_zz=0.3,
    ntrials=20, seed=43, verbose=false
)

println("L=8, P=0.3: B=$(round(result2.B, digits=4)), Purity=$(round(result2.purity_bar, digits=4))")

println("\n✓ Separate pipeline works on learning_to_trivial_lambda0.7!")
