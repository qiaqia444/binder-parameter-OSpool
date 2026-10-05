#!/usr/bin/env julia

"""
Generate parameter file for the Right Boundary dual EA Binder scan:
L = 48, 56 near-critical zoom

Dual-EA counterpart of `make_params_right_boundary_L48_L56.jl`: same L values,
q = 0.38:0.01:0.45, 160 x 25 split (4000 trajectories per (L, q)), and the same
seed_start = 150001, so each job evolves the same Born trajectories as the
corresponding Renyi-2 job. The matching run script
(`run_dual_ea_right_boundary_L48_L56.sh`) uses the same dynamics settings as the
Renyi-2 L48/L56 run (maxdim 512, cutoff 1e-13, T_max_factor 4, 2 CPUs).

Scans: q = P_x = P_zz from 0.38 to 0.45 in steps of 0.01 (8 points)
Fixed: lambda_x = 0.7, lambda_zz = 0.0
"""

using Printf

L_values = [48, 56]

lambda_x = 0.7
lambda_zz = 0.0

P_values = collect(0.38:0.01:0.45)

n_samples = 160
ntrials = 25   # 160 x 25 = 4000 trajectories per (L, q)

# Same as make_params_right_boundary_L48_L56.jl -> paired trajectories
seed_start = 150001

open("params_dual_ea_right_boundary_L48_L56.txt", "w") do f
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
println("Parameter file generated: params_dual_ea_right_boundary_L48_L56.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x, λ_zz fixed: $lambda_zz")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials, samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1) (same as params_right_boundary_L48_L56.txt)")
