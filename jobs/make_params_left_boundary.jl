#!/usr/bin/env julia

"""
Generate parameter file for Left Boundary Scan
Explores dephasing-induced phase transition at fixed measurement strength

Scans: P_x from 0 to 0.5 at fixed λ_x = 0.3
Fixed: λ_zz = 0, P_zz = 0 (no ZZ effects for cleanest physics)
"""

using Printf

# System sizes (multiple L for finite-size scaling)
L_values = [8, 10, 12, 14, 16]

# Fixed measurement strength (left boundary scan)
lambda_x = 0.3
lambda_zz = 0.0

# Dephasing probabilities to scan
# Dense sampling around expected critical region (0.2-0.4)
P_x_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

# P_zz always zero for left boundary scan
P_zz = 0.0

# Number of samples per configuration (for error bars)
n_samples = 10

# Number of trials per job
ntrials = 400  # More trials for better statistics

# Starting seed
seed_start = 7001

# Open output file
open("params_left_boundary.txt", "w") do f
    seed = seed_start
    
    # Group by L first, then scan over P_x
    # This makes it easier to plot Binder vs P_x for each L
    for L in L_values
        for P_x in P_x_values
            for sample in 1:n_samples
                out_prefix = @sprintf("left_boundary_L%d_lx%.2f_Px%.2f_s%d", 
                                     L, lambda_x, P_x, sample)
                
                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")
                
                seed += 1
            end
        end
    end
end

# Print summary
total_jobs = length(L_values) * length(P_x_values) * n_samples
println("Parameter file generated: params_left_boundary.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x")
println("λ_zz fixed: $lambda_zz")
println("P_x scan: ", P_x_values)
println("P_zz fixed: $P_zz")
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
