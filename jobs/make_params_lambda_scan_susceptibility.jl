#!/usr/bin/env julia

"""
Generate parameter file for Lambda-Path Scan (Rényi-2 / Edwards-Anderson
overlap susceptibilities, ITensorCorrelators.jl variant)

Sweeps the interior phase-diagram path

    lambda_x  = delta * lambda
    lambda_zz = delta * (1 - lambda)

Scans: lambda from 0.0 to 1.0 in steps of 0.1
Fixed: delta = 0.7, q = 0.1 (dephasing q_x = q_zz = q)
"""

using Printf

# System sizes (multiple L for finite-size scaling)
L_values = [8, 16, 24, 32]

# Path amplitude (FIXED)
delta = 0.7

# Fixed dephasing strength (q_x = q_zz = q)
q = 0.1

# Lambda values to scan: 0.0, 0.1, ..., 1.0
lambda_values = [round(l, digits=1) for l in 0.0:0.1:1.0]

# Number of samples per configuration (for error bars)
n_samples = 40  # More jobs for faster parallelization

# Number of trials per job
ntrials = 100  # Faster jobs, same total statistics (40×100=4000)

# Starting seed
seed_start = 80001

# Open output file
open("params_lambda_scan_susceptibility.txt", "w") do f
    seed = seed_start

    # Group by L first, then scan over lambda
    # This makes it easier to plot kappa2/kappaEA vs lambda for each L
    for L in L_values
        for lambda in lambda_values
            for sample in 1:n_samples
                out_prefix = @sprintf("lambda_scan_susceptibility_L%d_delta%.2f_q%.2f_lambda%.2f_s%d",
                                     L, delta, q, lambda, sample)

                # Format: L lambda delta q ntrials seed sample out_prefix
                println(f, "$L $lambda $delta $q $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

# Print summary
total_jobs = length(L_values) * length(lambda_values) * n_samples
println("Parameter file generated: params_lambda_scan_susceptibility.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("delta fixed: $delta")
println("q fixed: $q")
println("lambda scan: ", lambda_values)
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
