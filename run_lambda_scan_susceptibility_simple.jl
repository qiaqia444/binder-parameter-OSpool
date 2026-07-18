#!/usr/bin/env julia

"""
Lambda-Path Scan (ITensorCorrelators.jl) - local convenience script

Scan lambda from 0 to 1 at fixed delta, q, sweeping the interior
phase-diagram path

    lambda_x  = delta * lambda
    lambda_zz = delta * (1 - lambda)

Observables: trajectory-averaged kappa_2 (Rényi-2 replica-overlap
susceptibility, ITensorCorrelators.jl) and kappa_EA (Edwards-Anderson
susceptibility), computed by `lambda_scan_susceptibilities_itensorcorrelators.jl`,
shared with the cluster entry point `run_lambda_scan_susceptibility_scan.jl`.
"""

using Dates
using JSON
using Statistics

include("lambda_scan_susceptibilities_itensorcorrelators.jl")

# === EDIT THESE PARAMETERS ===
L = 8                    # System size
delta = 0.7             # Path amplitude (FIXED)
q = 0.1                 # Dephasing strength (FIXED)
lambda_steps = 11       # Number of lambda values
lambda_min = 0.0        # Minimum lambda
lambda_max = 1.0        # Maximum lambda

ntrials = 20            # Born-sampled Monte Carlo trajectories
maxdim = 256            # Max bond dimension (dynamics)
cutoff = 1e-10          # SVD truncation cutoff (dynamics)
T_max_factor = 2        # T_max = T_max_factor * L
pair_batch_size = 4096  # ITensorCorrelators.jl batch size
compute_EA = true       # Also compute the Edwards-Anderson susceptibility
seed = 1234             # Random seed
output_dir = "lambda_scan_susceptibility_results"
# =============================

println("="^70)
println("LAMBDA-PATH SCAN: Rényi-2 / Edwards-Anderson Overlap Susceptibilities")
println("="^70)
println("Physics: lambda_x = delta*lambda, lambda_zz = delta*(1-lambda), q_x=q_zz=q")
println()
println("Parameters:")
println("  L = $L")
println("  delta = $delta (FIXED)")
println("  q = $q (FIXED dephasing)")
println("  lambda scan: [$lambda_min, $lambda_max] with $lambda_steps points")
println("  Trajectories per point: $ntrials")
println("  T_max = $T_max_factor * L = $(T_max_factor * L)")
println("="^70)
println()

# Create output directory
mkpath(output_dir)

# Create lambda values
lambda_values = lambda_steps <= 1 ? [lambda_min] : range(lambda_min, lambda_max, length=lambda_steps)

results = []

for (i, lambda) in enumerate(lambda_values)
    println("\n[$i/$lambda_steps] lambda = $(round(lambda, digits=4))")
    println("-"^70)

    t_start = time()

    result = ls_run_lambda_point(
        L, Float64(lambda);
        delta = delta,
        q = q,
        ntrials = ntrials,
        T_max = T_max_factor * L,
        maxdim = maxdim,
        cutoff = cutoff,
        pair_batch_size = pair_batch_size,
        seed = seed + i,
        compute_EA = compute_EA,
    )

    t_elapsed = time() - t_start

    # Store results (trial-level arrays are intentionally not persisted)
    result_dict = Dict(
        "L" => L,
        "lambda" => lambda,
        "lambda_x" => result.lambda_x,
        "lambda_zz" => result.lambda_zz,
        "delta" => delta,
        "q" => q,
        "kappa2" => result.kappa2,
        "kappa2_se" => result.kappa2_se,
        "kappaEA" => result.kappaEA,
        "kappaEA_se" => result.kappaEA_se,
        "purity" => result.purity,
        "ntrials" => result.ntrials,
        "T_max" => result.T_max,
        "max_interphysical_linkdim" => result.max_interphysical_linkdim,
        "max_trace_error" => result.max_trace_error,
        "observable_method" => "itensorcorrelators",
        "time_seconds" => t_elapsed
    )

    push!(results, result_dict)

    # Print results
    println("  kappa2  = $(round(result.kappa2, digits=4)) ± $(round(result.kappa2_se, digits=4))")
    if compute_EA
        println("  kappaEA = $(round(result.kappaEA, digits=4)) ± $(round(result.kappaEA_se, digits=4))")
    end
    println("  Purity = $(round(result.purity, digits=4))")
    println("  Time: $(round(t_elapsed, digits=1))s")

    # Save intermediate results
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
    output_file = joinpath(output_dir, "lambda_scan_susceptibility_L$(L)_delta$(delta)_q$(q)_$(timestamp).json")
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
    println("  lambda=$(round(r["lambda"], digits=2)): kappa2=$(round(r["kappa2"], digits=4)), kappaEA=$(round(r["kappaEA"], digits=4))")
end
