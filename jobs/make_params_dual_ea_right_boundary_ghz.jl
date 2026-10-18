#!/usr/bin/env julia

"""
Generate parameter file for the Right Boundary dual Edwards-Anderson (EA) Binder scan with
a GHZ INITIAL STATE, |GHZ> = (|0...0> + |1...1>)/sqrt(2) (instead of the all-up state).

Mirrors `make_params_right_boundary_ghz.jl` (Renyi-2 GHZ family) job-for-job: same L values,
same q grid, same n_samples x ntrials split, and the SAME seed_start (120001) and seed
increments. `run_dual_ea_right_boundary_ghz_scan.jl` derives its per-point and per-trajectory
seeds like the Renyi-2 GHZ driver and its GHZ-start evolution reproduces
`renyi2_right_boundary_ghz_core.jl` exactly (verified bit for bit), so job k here evolves the
SAME trajectories as job k of params_right_boundary_ghz.txt. Only output filenames differ
(`dual_ea_right_boundary_ghz_...`), so nothing collides with the all-up dual-EA results or
the Renyi-2 GHZ results.

Scans: q = P_x = P_zz from 0 to 0.5
Fixed: lambda_x = 0.7 (X measurement strength), lambda_zz = 0.0 (no ZZ measurements)
"""

using Printf

# System sizes (as in the Renyi-2 GHZ family; drop L = 8 in the analysis if wanted)
L_values = [8, 16, 24, 32]

lambda_x = 0.7
lambda_zz = 0.0

P_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

n_samples = 40
ntrials = 100   # 40 x 100 = 4000 trajectories per (L, q)

# Same starting seed as make_params_right_boundary_ghz.jl -> paired trajectories
seed_start = 120001

open("params_dual_ea_right_boundary_ghz.txt", "w") do f
    seed = seed_start

    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                P_x = P
                P_zz = P

                out_prefix = @sprintf("dual_ea_right_boundary_ghz_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                                     L, lambda_x, lambda_zz, P_x, sample)

                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")

                seed += 1
            end
        end
    end
end

total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_dual_ea_right_boundary_ghz.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x, λ_zz fixed: $lambda_zz")
println("q = P_x = P_zz scan: ", P_values)
println("Trials per job: $ntrials, samples per configuration: $n_samples")
println("Seeds: $seed_start to $(seed_start + total_jobs - 1) (same as params_right_boundary_ghz.txt)")
println("Initial state: GHZ. Observable: dual EA Binder, open-chain endpoints a = 2..L")
