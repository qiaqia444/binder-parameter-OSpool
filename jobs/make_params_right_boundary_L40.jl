#!/usr/bin/env julia

"""
Generate parameter file for Right Boundary Scan: L = 40 EXTENSION ONLY

Adds a fifth, larger system size (L=40) to the existing base right-boundary
dataset (L=8,16,24,32 from `make_params_right_boundary.jl`). Same physics
(no GHZ initial state - uses the standard product-state initial condition
via the existing `renyi2_right_boundary_core.jl` / `run_right_boundary_scan.jl`),
same P scan, same n_samples/ntrials. Kept as a SEPARATE params/job file (not
merged into `make_params_right_boundary.jl`) so re-running this does not
change the seeds/output filenames of the already-collected L=8,16,24,32 jobs.

Scans: P_x from 0 to 0.5 (X dephasing)
Fixed: λ_x = 0.7 (X measurement strength), λ_zz = 0.0 (no ZZ measurements)
Equal dephasing: P_x = P_zz
"""

using Printf

# Only the new system size.
L_values = [40]

# Fixed measurement strengths
lambda_x = 0.7   # X measurement strength
lambda_zz = 0.0  # No ZZ measurements

# Dephasing probabilities to scan (X dephasing only)
# Dense sampling around expected critical region (0.2-0.4)
P_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

# Number of samples per configuration (for error bars)
n_samples = 40  # More jobs for faster parallelization

# Number of trials per job
ntrials = 100  # Faster jobs, same total statistics (40×100=4000)

# Starting seed (kept disjoint from other families; see repo memory registry)
seed_start = 130001

# Open output file
open("params_right_boundary_L40.txt", "w") do f
    seed = seed_start
    
    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                # P_x = P_zz (equal dephasing)
                P_x = P
                P_zz = P
                
                out_prefix = @sprintf("right_boundary_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d", 
                                     L, lambda_x, lambda_zz, P_x, sample)
                
                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")
                
                seed += 1
            end
        end
    end
end

# Print summary
total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_right_boundary_L40.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x (X measurement strength)")
println("λ_zz fixed: $lambda_zz (no ZZ measurements)")
println("P_x scan: ", P_values)
println("P_zz: 0.0 (no ZZ dephasing)")
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
println("Initial state: standard product state (NOT GHZ)")
