#!/usr/bin/env julia

"""
Quick Left Boundary Scan - Simple Version

Scan P_x at fixed λ_x (no ZZ measurements or dephasing)
This is the cleanest left boundary physics!
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

# === EDIT THESE PARAMETERS ===
L = 8                    # System size
lambda_x = 0.3          # X measurement strength (FIXED)
lambda_zz = 0.0         # ZZ measurement strength (0 for pure X scan)
P_steps = 11            # Number of P values
P_min = 0.0             # Minimum P_x
P_max = 0.5             # Maximum P_x
P_zz_mode = "zero"      # "zero", "match", or fixed value

ntrials = 500           # Monte Carlo trajectories
maxdim = 256            # Max bond dimension
seed = 42               # Random seed
output_dir = "left_boundary_results"
# =============================

println("="^70)
println("LEFT BOUNDARY SCAN: X Measurements vs X Dephasing")
println("="^70)
println("Physics: Pure competition - no ZZ terms")
println()
println("Parameters:")
println("  L = $L")
println("  λ_x = $lambda_x (FIXED)")
println("  λ_zz = $lambda_zz (no ZZ measurements)")
println("  P_x scan: [$P_min, $P_max] with $P_steps points")
println("  P_zz = 0 (no ZZ dephasing)")
println("  Trajectories per point: $ntrials")
println("="^70)
println()

# Create output directory
mkpath(output_dir)

# Create P values
P_values = range(P_min, P_max, length=P_steps)

results = []

for (i, P_x) in enumerate(P_values)
    P_zz = 0.0  # Always zero for pure X vs X competition
    
    println("\n[$i/$P_steps] P_x = $(round(P_x, digits=3))")
    println("-"^70)
    
    t_start = time()
    
    result = ea_binder_density_matrix(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        maxdim = maxdim,
        cutoff = 1e-12,
        seed = seed + i
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
        "B_mean" => result.B_mean_of_trials,
        "B_std" => result.B_std_of_trials,
        "S2_bar" => result.S2_bar,
        "S4_bar" => result.S4_bar,
        "ntrials" => ntrials,
        "time_seconds" => t_elapsed
    )
    
    push!(results, result_dict)
    
    # Print results
    println("  B = $(round(result.B, digits=4))")
    println("  B_mean = $(round(result.B_mean_of_trials, digits=4)) ± $(round(result.B_std_of_trials, digits=4))")
    println("  Time: $(round(t_elapsed, digits=1))s")
    
    # Save intermediate results
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
    output_file = joinpath(output_dir, "left_boundary_L$(L)_lambdax$(lambda_x)_$(timestamp).json")
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
    println("  P_x=$(round(r["P_x"], digits=2)): B=$(round(r["B"], digits=4))")
end
