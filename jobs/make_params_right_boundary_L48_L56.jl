#!/usr/bin/env julia

"""
Generate parameter file for Right Boundary Scan: L = 48, 56 near-critical zoom

Adds two new, larger system sizes (L=48,56) to the right-boundary family,
restricted to a narrow q window around the expected critical region
(q=0.38:0.01:0.45) instead of the full [0,0.5] scan used by
`make_params_right_boundary.jl` / `make_params_right_boundary_L40.jl`. Same
physics (standard product-state initial condition, no GHZ) via the existing
`renyi2_right_boundary_core.jl` / `run_right_boundary_scan.jl` - no core/scan
changes needed. Kept as a SEPARATE params/job file so re-running this does
not change the seeds/output filenames of any other right-boundary family.

Scans: P_x = P_zz from 0.38 to 0.45 in steps of 0.01 (8 points)
Fixed: λ_x = 0.7 (X measurement strength), λ_zz = 0.0 (no ZZ measurements)

n_samples/ntrials split to 160x25 (was 40x100) to make more, smaller jobs -
same 4000 total trials per (L,q) point - matching request_cpus=2 (was 4) in
jobs_right_boundary_L48_L56.submit, so more jobs can run concurrently.
"""

using Printf

# Only the two new, larger system sizes.
L_values = [48, 56]

# Fixed measurement strengths
lambda_x = 0.7   # X measurement strength
lambda_zz = 0.0  # No ZZ measurements

# Narrow near-critical q window (X dephasing only, P_x = P_zz)
P_values = collect(0.38:0.01:0.45)

# Number of samples per configuration (for error bars)
n_samples = 160  # More, smaller jobs so more can run concurrently at request_cpus=2

# Number of trials per job
ntrials = 25  # Same total statistics as before (160×25=4000)

# Starting seed (kept disjoint from other families; see repo memory registry)
seed_start = 150001

# Open output file
open("params_right_boundary_L48_L56.txt", "w") do f
    seed = seed_start

    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                # P_x = P_zz (equal dephasing)
                P_x = P
                P_zz = P

                out_prefix = @sprintf("right_boundary_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                                     L, lambda_x, lambda_zz, P_x, sample)

                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

# Print summary
total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_right_boundary_L48_L56.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x (X measurement strength)")
println("λ_zz fixed: $lambda_zz (no ZZ measurements)")
println("P_x scan: ", P_values)
println("P_zz: equal to P_x")
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
println("Initial state: standard product state (NOT GHZ)")
