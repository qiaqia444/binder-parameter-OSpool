#!/usr/bin/env julia

"""
Generate parameter file for the Right Boundary Renyi-2 Binder scan measured AFTER THE LAST
X-DEPHASING step: L = 48, 56 near-critical zoom (q = 0.38:0.01:0.45, 160 x 25 split)

Counterpart of `make_params_right_boundary_L48_L56.jl`: same L values, same q grid, same
n_samples x ntrials split, and the SAME seed_start = 150001, so job k evolves the SAME Born
trajectories as job k of `params_right_boundary_L48_L56.txt` (the production runs, which measure after
the last ZZ dephasing). The two observation times can therefore be compared trajectory by
trajectory. Only the output filenames differ (`renyi2_xnoise_right_boundary_...`), so nothing
collides with the production Renyi-2 or the dual-EA results.

Scans: q = P_x = P_zz
Fixed: lambda_x = 0.7 (X measurement strength), lambda_zz = 0.0 (no ZZ measurements)
"""

using Printf

L_values = [48, 56]

lambda_x = 0.7
lambda_zz = 0.0

P_values = collect(0.38:0.01:0.45)

n_samples = 160
ntrials = 25

# Same starting seed as make_params_right_boundary_L48_L56.jl -> paired trajectories
seed_start = 150001

open("params_renyi2_right_boundary_xnoise_L48_L56.txt", "w") do f
    seed = seed_start

    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                P_x = P
                P_zz = P

                out_prefix = @sprintf("renyi2_xnoise_right_boundary_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                                     L, lambda_x, lambda_zz, P_x, sample)

                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_renyi2_right_boundary_xnoise_L48_L56.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x, λ_zz fixed: $lambda_zz")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials, samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1) (same as params_right_boundary_L48_L56.txt)")
