#!/usr/bin/env julia

"""
Generate parameter file for right-boundary finite-time test.

Goal: compare Binder at T_max=4L vs T_max=8L for selected system sizes and probabilities.

Grid:
- L in {16, 24, 32}
- P_x = P_zz in {0.25, 0.30, 0.35, 0.40, 0.45}
- T_max_factor in {4, 8}
"""

using Printf

L_values = [16, 24, 32]
lambda_x = 0.7
lambda_zz = 0.0
P_values = [0.25, 0.30, 0.35, 0.40, 0.45]
T_max_factors = [4, 8]

n_samples = 40
ntrials = 100

# Keep seed range disjoint from other families.
seed_start = 140001

open("params_right_boundary_tmax.txt", "w") do f
    seed = seed_start

    for L in L_values
        for P in P_values
            for Tfac in T_max_factors
                for sample in 1:n_samples
                    P_x = P
                    P_zz = P
                    out_prefix = @sprintf(
                        "right_boundary_tmax_L%d_Tf%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                        L, Tfac, lambda_x, lambda_zz, P_x, sample,
                    )

                    # Format:
                    # L lambda_x lambda_zz P_x P_zz T_max_factor ntrials seed sample out_prefix
                    println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $Tfac $ntrials $seed $sample $out_prefix")
                    seed += 1
                end
            end
        end
    end
end

total_jobs = length(L_values) * length(P_values) * length(T_max_factors) * n_samples
println("Parameter file generated: params_right_boundary_tmax.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("P values: ", P_values)
println("T_max factors: ", T_max_factors)
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
