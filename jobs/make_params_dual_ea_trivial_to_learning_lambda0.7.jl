#!/usr/bin/env julia

"""
Generate parameter file for the Trivial-to-Learning dual Edwards-Anderson (EA) Binder scan
at lambda = 0.7 (delta = 0.7: lambda_x = 0.49, lambda_zz = 0.21; both weak-measurement
layers active).

Mirrors the layout of `make_params_dual_ea_right_boundary.jl`: one job per
(L, q, sample), 40 samples x 100 trajectories = 4000 trajectories per (L, q).
`run_dual_ea_trivial_to_learning_lambda0.7_scan.jl` reuses the Born-sampling dynamics
`lt_evolve_one_trial` of the Renyi-2 learning-to-trivial family and its per-trajectory
seed derivation, so a job with the same seed evolves the same trajectories as the
Renyi-2 scan at the same (lambda_x, lambda_zz). Output names are prefixed
`dual_ea_trivial_to_learning_lambda0.7_...` (the lambda is part of the name so other
lambda families stay separate, and nothing collides with the Renyi-2
`learning_to_trivial_lambda0.7_...` outputs).

Scans: q = P_x = P_zz from 0 to 0.5
Fixed: lambda_x = 0.49 (X measurement strength), lambda_zz = 0.21 (ZZ measurement strength)
"""

using Printf

# System sizes (L = 8 is intentionally not included)
L_values = [16, 24, 32, 40]

# Fixed measurement strengths (lambda = 0.7, delta = 0.7)
lambda_x = 0.49   # X measurement strength
lambda_zz = 0.21  # ZZ measurement strength

# Dephasing probabilities to scan (q = P_x = P_zz)
P_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

# Number of samples per configuration (for error bars)
n_samples = 40

# Number of trials per job
ntrials = 100  # 40 x 100 = 4000 trajectories per (L, q)

seed_start = 73001

open("params_dual_ea_trivial_to_learning_lambda0.7.txt", "w") do f
    seed = seed_start

    # Group by L first, then scan over P
    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                P_x = P
                P_zz = P

                out_prefix = @sprintf("dual_ea_trivial_to_learning_lambda0.7_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                                     L, lambda_x, lambda_zz, P_x, sample)

                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_dual_ea_trivial_to_learning_lambda0.7.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x (X measurement strength)")
println("λ_zz fixed: $lambda_zz (ZZ measurement strength)")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1)")
println("Observable: dual EA Binder, open-chain endpoints a = 2..L")
