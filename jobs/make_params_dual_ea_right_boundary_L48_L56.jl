#!/usr/bin/env julia

"""
Generate parameter file for the Right Boundary dual EA Binder scan:
L = 48, 56 EXTENSION ONLY, union of the coarse and fine q grids in ONE file.

Dual-EA counterpart of the L = 40 extension (`make_params_dual_ea_right_boundary_L40.jl`),
covering the union of the two existing q grids, with every q simulated only ONCE (19 values):
  - coarse prefix: the coarse-grid values outside the fine range, q = 0, 0.05, 0.10, 0.15, 0.20,
            0.40, 0.45, 0.50  (base family `make_params_dual_ea_right_boundary.jl`),
            output prefix `dual_ea_right_boundary_L...`,      seeds 150001-152560
  - fine prefix: q = 0.25:0.01:0.35 (fine family `make_params_dual_ea_right_boundary_fine.jl`;
            includes the coarse-grid points 0.25, 0.30, 0.35),
            output prefix `dual_ea_right_boundary_fine_L...`, seeds 200001-203520
Each (L, q) gets 4000 trajectories, split 160 x 25 instead of 40 x 100 so each L = 48, 56 job
stays short. The prefixes are the same as in the L = 8..40 families, so
`collect_dual_ea_right_boundary_results.sh` picks up the coarse files and
`collect_dual_ea_right_boundary_fine_results.sh` the fine ones; pooled per (L, q) they give the
full 0:0.05:0.5 + 0.25:0.01:0.35 grid. (For L <= 40 the q = 0.25, 0.30, 0.35 points were run
twice and pooled; for L = 48, 56 they are run once to save cluster time.)
Kept as a SEPARATE params/job file so re-running it does not change the seeds/output
filenames of the L = 8..40 jobs (seeds used there: 9001-10760, 130001-130440, 170001-171760).
The matching run script (`run_dual_ea_right_boundary_L48_L56.sh`) uses maxdim 512,
cutoff 1e-13, T_max_factor 4, 2 CPUs.

Fixed: lambda_x = 0.7, lambda_zz = 0.0
"""

using Printf

L_values = [48, 56]

lambda_x = 0.7
lambda_zz = 0.0

n_samples = 160
ntrials = 25   # 160 x 25 = 4000 trajectories per (L, q)

# (output prefix, q grid, first seed)
grids = [
    ("dual_ea_right_boundary",      [0.0, 0.05, 0.10, 0.15, 0.20, 0.40, 0.45, 0.50], 150001),
    ("dual_ea_right_boundary_fine", collect(0.25:0.01:0.35),                         200001),
]

total_jobs = 0

open("params_dual_ea_right_boundary_L48_L56.txt", "w") do f
    for (prefix, P_values, seed_start) in grids
        seed = seed_start

        for L in L_values
            for P in P_values
                for sample in 1:n_samples
                    P_x = P
                    P_zz = P

                    out_prefix = @sprintf("%s_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                                         prefix, L, lambda_x, lambda_zz, P_x, sample)

                    # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                    println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                    seed += 1
                    global total_jobs += 1
                end
            end
        end

        println("$prefix: q = ", P_values, ", seeds $seed_start to $(seed - 1)")
    end
end

println("Parameter file generated: params_dual_ea_right_boundary_L48_L56.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x, λ_zz fixed: $lambda_zz")
println("Trials per job: $ntrials, samples per configuration: $n_samples")
