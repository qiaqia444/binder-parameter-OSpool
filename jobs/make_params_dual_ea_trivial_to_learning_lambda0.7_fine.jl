#!/usr/bin/env julia

"""
Generate parameter file for the Trivial-to-Learning dual Edwards-Anderson (EA) Binder scan
at lambda = 0.7 (delta = 0.7: lambda_x = 0.49, lambda_zz = 0.21; both weak-measurement
layers active): FINE q grid around the crossing region, q = 0.05:0.005:0.10.

Same observable, dynamics and driver as the coarse-grid family
(`run_dual_ea_trivial_to_learning_lambda0.7_scan.jl`,
`dual_ea_trivial_to_learning_lambda0.7_core.jl`, driver defaults maxdim = 256,
cutoff = 1e-12, T_max_factor = 4), with the same 40 x 100 split (4000 trajectories per
(L, q)). Only the q grid differs: 11 points with step 0.005 instead of the coarse 0.05 grid,
containing the crossing of the coarse data at q ~ 0.078.

Seeds start at 180001, a range disjoint from every other family in this repository, so the
q = 0.05 and q = 0.10 points are INDEPENDENT re-samples of the coarse-grid points.
Output names are prefixed `dual_ea_trivial_to_learning_lambda0.7_fine_...` (note the three
decimals in `Px`, needed for the 0.005 step), so nothing collides with the coarse-grid
results, other lambda families or the Renyi-2 outputs, and the coarse collect/analysis
scripts do not pick these files up.

Scans: q = P_x = P_zz from 0.05 to 0.10 in steps of 0.005 (11 points)
Fixed: lambda_x = 0.49 (X measurement strength), lambda_zz = 0.21 (ZZ measurement strength)
"""

using Printf

L_values = [16, 24, 32, 40]

lambda_x = 0.49
lambda_zz = 0.21

P_values = round.(collect(0.05:0.005:0.10); digits=3)

n_samples = 40
ntrials = 100   # 40 x 100 = 4000 trajectories per (L, q)

seed_start = 180001

open("params_dual_ea_trivial_to_learning_lambda0.7_fine.txt", "w") do f
    seed = seed_start

    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                P_x = P
                P_zz = P

                out_prefix = @sprintf("dual_ea_trivial_to_learning_lambda0.7_fine_L%d_lx%.2f_lzz%.2f_Px%.3f_s%d",
                                     L, lambda_x, lambda_zz, P_x, sample)

                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_dual_ea_trivial_to_learning_lambda0.7_fine.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x, λ_zz fixed: $lambda_zz")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials, samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1)")
println("Observable: dual EA Binder, open-chain endpoints a = 2..L")
