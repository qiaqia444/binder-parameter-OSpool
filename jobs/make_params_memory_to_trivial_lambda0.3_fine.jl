#!/usr/bin/env julia

"""
Generate parameter file for the Memory-to-Trivial EA Binder scan at lambda = 0.3
(delta = 0.7: lambda_x = 0.21, lambda_zz = 0.49) with a FINE dephasing grid around the
transition, P = 0.025:0.005:0.075.

Follows `make_params_memory_to_trivial_lx0.3.jl` (same physics slice: lambda_x = 0.21,
lambda_zz = 0.49; same 40 samples x 100 trials = 4000 trajectories per
(L, P), same driver `run_memory_to_trivial_scan.jl`, i.e. the plain EA Binder from
density-matrix evolution - NOT the dual-EA Binder and NOT the Renyi-2 Binder). The P
grid is fine (11 points with step 0.005 around pc ~ 0.05) and the sizes are L = 24, 32, 40.

Seeds start at 190001, disjoint from every other params file in this repository. Output
names carry the lambda in the prefix (`memory_to_trivial_lambda0.3_fine_...`) so other
lambda families can be told apart. The prefix differs from the coarse-grid
`memory_to_trivial_L..._P...` outputs of `params_memory_to_trivial_lx0.3.txt` (same physical
slice), so those and the Renyi-2 results are never matched; `Px` has three decimals because
of the 0.005 step.

Scans: P_x = P_zz from 0.025 to 0.075 in steps of 0.005 (coupled dephasing, 11 points)
Fixed: lambda_x = 0.21 (X measurements), lambda_zz = 0.49 (ZZ measurements)
"""

using Printf

# System sizes (L = 24, 32, 40)
L_values = [24, 32, 40]

# Fixed measurement strengths (lambda = 0.3, delta = 0.7)
lambda_x = 0.21   # X measurement strength
lambda_zz = 0.49  # ZZ measurement strength

# Dephasing probabilities to scan (P_x = P_zz coupled)
P_values = round.(collect(0.025:0.005:0.075); digits=3)

# Number of samples per configuration (for error bars)
n_samples = 40

# Number of trials per job
ntrials = 100  # 40 x 100 = 4000 trajectories per (L, P)

# Starting seed
seed_start = 190001

# Open output file
open("params_memory_to_trivial_lambda0.3_fine.txt", "w") do f
    seed = seed_start

    # Group by L first, then scan over P
    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                # P_x and P_zz are the same (coupled dephasing)
                P_x = P
                P_zz = P

                out_prefix = @sprintf("memory_to_trivial_lambda0.3_fine_L%d_lx%.2f_lzz%.2f_Px%.3f_s%d",
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
println("Parameter file generated: params_memory_to_trivial_lambda0.3_fine.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x (X measurements)")
println("λ_zz fixed: $lambda_zz (ZZ measurements)")
println("P scan (P_x = P_zz): ", P_values)
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1)")
