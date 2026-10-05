#!/usr/bin/env julia

"""
Generate parameter file for the Right Boundary dual EA Binder scan: L = 40 EXTENSION ONLY

Dual-EA counterpart of `make_params_right_boundary_L40.jl` (same L, q grid,
n_samples/ntrials, and the same seed_start = 130001, so each job evolves the same
Born trajectories as the corresponding Renyi-2 job). Kept as a SEPARATE params/job file so
re-running it does not change the seeds/output filenames of the L = 8..32 dual-EA jobs.

Scans: q = P_x = P_zz from 0 to 0.5
Fixed: lambda_x = 0.7, lambda_zz = 0.0
"""

using Printf

# Only the new system size.
L_values = [40]

lambda_x = 0.7
lambda_zz = 0.0

P_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

n_samples = 40
ntrials = 100   # 40 x 100 = 4000 trajectories per (L, q)

# Same as make_params_right_boundary_L40.jl -> paired trajectories
seed_start = 130001

open("params_dual_ea_right_boundary_L40.txt", "w") do f
    seed = seed_start

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
println("Parameter file generated: params_dual_ea_right_boundary_L40.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x, λ_zz fixed: $lambda_zz")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials, samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1) (same as params_right_boundary_L40.txt)")
