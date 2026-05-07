#!/usr/bin/env julia
"""
Debug: Check if purity is changing correctly with dephasing
"""

using Random, Statistics, ITensors, ITensorMPS

include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_dynamics_density_matrix_1.jl")

lambda_x = 0.49
lambda_zz = 0.21

println("DEBUG: Checking purity vs P")
println("="^60)

for P in [0.0, 0.1, 0.3, 0.5]
    result = renyi2_binder_density_matrix_separate(
        8; lambda_x=lambda_x, lambda_zz=lambda_zz, P_x=P, P_zz=P,
        ntrials=20, seed=42, verbose=false
    )
    
    println("P=$P: Purity=$(round(result.purity_bar, digits=6)), Min=$(round(minimum(result.purities[isfinite.(result.purities)]), digits=6)), Max=$(round(maximum(result.purities[isfinite.(result.purities)]), digits=6))")
end

println("\nDEBUG: Single trial analysis")
println("="^60)

state, _ = evolve_renyi2_density_matrix_one_trial_separate(
    8; lambda_x=lambda_x, lambda_zz=lambda_zz, P_x=0.0, P_zz=0.0,
    maxdim=256, cutoff=1e-12, rng=MersenneTwister(42)
)
ρ0 = get_mps(state)
tr0 = doubledtrace(ρ0)
ρ0_norm = ρ0 / tr0
purity0 = real(inner(ρ0_norm, ρ0_norm))
println("P=0.0: Purity=$purity0")

state, _ = evolve_renyi2_density_matrix_one_trial_separate(
    8; lambda_x=lambda_x, lambda_zz=lambda_zz, P_x=0.3, P_zz=0.3,
    maxdim=256, cutoff=1e-12, rng=MersenneTwister(42)
)
ρ03 = get_mps(state)
tr03 = doubledtrace(ρ03)
ρ03_norm = ρ03 / tr03
purity03 = real(inner(ρ03_norm, ρ03_norm))
println("P=0.3: Purity=$purity03")
