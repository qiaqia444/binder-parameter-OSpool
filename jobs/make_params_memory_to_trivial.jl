#!/usr/bin/env julia

"""
Generate parameter file for Memory-to-Trivial Transition Scan
Explores the phase transition from memory-keeping to trivial phase

Scans: P_x = P_zz from 0 to 0.5 (coupled dephasing)
Fixed: λ_x = 0.1 (X measurements), λ_zz = 0.7 (ZZ measurements)
"""

using Printf

# System sizes (multiple L for finite-size scaling)
L_values = [8, 10, 12, 14, 16]

# Fixed measurement strengths
lambda_x = 0.1   # X measurement strength
lambda_zz = 0.7  # ZZ measurement strength

# Dephasing probabilities to scan (P_x = P_zz coupled)
# Dense sampling around expected critical region
P_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

# Number of samples per configuration (for error bars)
n_samples = 40  # More jobs for faster parallelization

# Number of trials per job
ntrials = 100  # Faster jobs, same total statistics (40×100=4000)

# Starting seed
seed_start = 9001

# Open output file
open("params_memory_to_trivial.txt", "w") do f
    seed = seed_start
    
    # Group by L first, then scan over P
    # This makes it easier to plot Binder vs P for each L
    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                # P_x and P_zz are the same (coupled dephasing)
                P_x = P
                P_zz = P
                
                out_prefix = @sprintf("memory_to_trivial_L%d_lx%.2f_lzz%.2f_P%.2f_s%d", 
                                     L, lambda_x, lambda_zz, P, sample)
                
                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                println(f, "$L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix")
                
                seed += 1
            end
        end
    end
end

# Print summary
total_jobs = length(L_values) * length(P_values) * n_samples
println("Parameter file generated: params_memory_to_trivial.txt")
println("Total jobs: $total_jobs")
println("L values: ", L_values)
println("λ_x fixed: $lambda_x (X measurements)")
println("λ_zz fixed: $lambda_zz (ZZ measurements)")
println("P scan (P_x = P_zz): ", P_values)
println("Trials per job: $ntrials")
println("Samples per configuration: $n_samples")
