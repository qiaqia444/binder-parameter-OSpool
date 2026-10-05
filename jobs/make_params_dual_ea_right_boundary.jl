#!/usr/bin/env julia

"""
Generate parameter file for the Right Boundary dual Edwards-Anderson (EA) Binder scan
(base family: L = 8, 16, 24, 32).

Mirrors `make_params_right_boundary.jl` (Renyi-2 Binder) job-for-job: same L values,
same q grid, same n_samples x ntrials split, and the SAME seed_start (9001) and seed
increments. `run_dual_ea_right_boundary_scan.jl` derives its per-point and
per-trajectory seeds exactly like `run_right_boundary_scan.jl` and reuses the same
Born-sampling dynamics, so job k here evolves the SAME trajectories as job k of
params_right_boundary.txt (given the same maxdim/cutoff/T_max, which are the driver
defaults in both). That makes the dual-EA Binder directly comparable, trajectory by
trajectory, with the Renyi-2 Binder. Only output filenames differ
(`dual_ea_right_boundary_...`), so nothing collides with the Renyi-2 results.

Scans: q = P_x = P_zz from 0 to 0.5
Fixed: lambda_x = 0.7 (X measurement strength), lambda_zz = 0.0 (no ZZ measurements)
"""

using Printf

# System sizes (multiple L for finite-size scaling)
L_values = [8, 16, 24, 32]

# Fixed measurement strengths
lambda_x = 0.7   # X measurement strength
lambda_zz = 0.0  # No ZZ measurements

# Dephasing probabilities to scan (q = P_x = P_zz)
P_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

# Number of samples per configuration (for error bars)
n_samples = 40

# Number of trials per job
ntrials = 100  # 40 x 100 = 4000 trajectories per (L, q)

# Same starting seed as make_params_right_boundary.jl -> paired trajectories
seed_start = 9001

open("params_dual_ea_right_boundary.txt", "w") do f
    seed = seed_start

    # Group by L first, then scan over P
    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                P_x = P
                P_zz = P

                out_prefix = @sprintf("dual_ea_right_boundary_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                                     L, lambda_x, lambda_zz, P_x, sample)

                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_dual_ea_right_boundary.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x (X measurement strength)")
println("λ_zz fixed: $lambda_zz (no ZZ measurements)")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1) (same as params_right_boundary.txt)")
println("Observable: dual EA Binder, open-chain endpoints a = 2..L")
