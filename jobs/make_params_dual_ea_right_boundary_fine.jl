#!/usr/bin/env julia

"""
Generate parameter file for the Right Boundary dual Edwards-Anderson (EA) Binder scan:
FINE q grid around the crossing region, q = 0.25:0.01:0.35.

Same observable, dynamics and driver as the base family
(`run_dual_ea_right_boundary_scan.jl`, `dual_ea_right_boundary_core.jl`, driver defaults
maxdim = 256, cutoff = 1e-12, T_max_factor = 4), with the same 40 x 100 split (4000
trajectories per (L, q)). Only the q grid differs: 11 points with step 0.01 instead of
the coarse 0.05 grid. L = 8 is not included.

Seeds start at 170001, a range disjoint from every existing dual-EA family
(9001-10760, 130001-130440, 150001-152560), so the q = 0.25, 0.30, 0.35 points are
INDEPENDENT re-samples of the coarse-grid points rather than copies of them.
Output names are prefixed `dual_ea_right_boundary_fine_...`, so nothing collides with the
coarse-grid results and the coarse collect/analysis scripts do not pick these files up.

Scans: q = P_x = P_zz from 0.25 to 0.35 in steps of 0.01 (11 points)
Fixed: lambda_x = 0.7 (X measurement strength), lambda_zz = 0.0 (no ZZ measurements)
"""

using Printf

L_values = [16, 24, 32, 40]

lambda_x = 0.7
lambda_zz = 0.0

P_values = collect(0.25:0.01:0.35)

n_samples = 40
ntrials = 100   # 40 x 100 = 4000 trajectories per (L, q)

seed_start = 170001

open("params_dual_ea_right_boundary_fine.txt", "w") do f
    seed = seed_start

    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                P_x = P
                P_zz = P

                out_prefix = @sprintf("dual_ea_right_boundary_fine_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                                     L, lambda_x, lambda_zz, P_x, sample)

                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_dual_ea_right_boundary_fine.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x, λ_zz fixed: $lambda_zz")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials, samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1)")
println("Observable: dual EA Binder, open-chain endpoints a = 2..L")
