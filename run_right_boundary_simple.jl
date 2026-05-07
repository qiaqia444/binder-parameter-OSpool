#!/usr/bin/env julia

"""
Right Boundary Scan - Rényi-2 Binder Focus

Scan P_x at fixed λ_x = 0.7, λ_zz = 0
Focus on measurement-induced entanglement transition using Rényi-2 Binder
"""

using Dates
using JSON
using Random
using Statistics
using ITensors, ITensorMPS
using LinearAlgebra

# Load modules
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_dynamics_density_matrix_1.jl")

# === EDIT THESE PARAMETERS ===
L = 8                    # System size
lambda_x = 0.7          # X measurement strength (FIXED)
lambda_zz = 0.0         # ZZ measurement strength (FIXED - no ZZ)
P_steps = 11            # Number of P values
P_min = 0.0             # Minimum P_x
P_max = 0.5             # Maximum P_x
P_zz_mode = "zero"      # Always zero for this boundary

ntrials = 500           # Monte Carlo trajectories
maxdim = 256            # Max bond dimension
seed = 42               # Random seed
output_dir = "right_boundary_results"
# =============================

println("="^70)
println("RIGHT BOUNDARY SCAN: Rényi-2 Binder Analysis")
println("="^70)
println("Physics: Measurement-induced transition (λ_x = 0.7)")
println()
println("Parameters:")
println("  L = $L")
println("  λ_x = $lambda_x (X measurement strength - FIXED)")
println("  λ_zz = $lambda_zz (no ZZ measurements)")
println("  P_x scan: [$P_min, $P_max] with $P_steps points")
println("  P_zz = 0 (no ZZ dephasing)")
println("  Trajectories per point: $ntrials")
println("  Focus: Rényi-2 Binder (M2, M4 moments)")
println("="^70)
println()

# Create output directory
mkpath(output_dir)

# Create P values
P_values = range(P_min, P_max, length=P_steps)

results = []

for (i, P_x) in enumerate(P_values)
    P_zz = 0.0  # Always zero for right boundary
    
    println("\n[$i/$P_steps] P_x = $(round(P_x, digits=3))")
    println("-"^70)
    
    t_start = time()
    
    result = renyi2_binder_density_matrix(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        maxdim = maxdim,
        cutoff = 1e-12,
        seed = seed + i,
        use_optimized = true,
        verbose = false
    )
    
    t_elapsed = time() - t_start
    
    # Store results
    result_dict = Dict(
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
        "time_seconds" => t_elapsed
    )
    
    push!(results, result_dict)
    
    # Print results
    println("  B (ensemble) = $(round(result.B, digits=4))")
    println("  B (mean of trials) = $(round(result.B_mean_of_trials, digits=4)) ± $(round(result.B_std_of_trials, digits=4))")
    println("  M2_bar = $(round(result.M2_bar, digits=4))")
    println("  M4_bar = $(round(result.M4_bar, digits=4))")
    println("  Purity = $(round(result.purity_bar, digits=4))")
    println("  Valid/Total: $(result.n_valid)/$(result.ntrials)")
    println("  Time: $(round(t_elapsed, digits=1))s")
    
    # Save intermediate results
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
    output_file = joinpath(output_dir, "right_boundary_L$(L)_lambdax$(lambda_x)_$(timestamp).json")
    open(output_file, "w") do f
        JSON.print(f, results, 4)
    end
end

println("\n" * "="^70)
println("✓ SCAN COMPLETE")
println("="^70)
println("Results saved to: $output_dir/")

# Print summary
println("\nSummary:")
for (i, r) in enumerate(results)
    println("  P_x=$(round(r["P_x"], digits=2)): B=$(round(r["B"], digits=4)), M2=$(round(r["M2_bar"], digits=4)), M4=$(round(r["M4_bar"], digits=4))")
end
