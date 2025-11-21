#!/usr/bin/env julia

"""
Generate parameter file for Dephasing P=0.2 jobs (L=8,12,16,20) - CORRECTED VERSION
Increased parallelization for faster completion
"""

using Printf

# System sizes for small systems
L_values = [8, 12, 16, 20]

# Lambda values (17 points, dense near critical point)
lambda_values = [0.1, 0.2, 0.3, 0.4, 0.45, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95]

# Dephasing probability
P_x = 0.2

# Number of samples per configuration
# CORRECTED VERSION: More parallel jobs for faster completion
n_samples = 10  # Increased from 3

# Number of trials per job
ntrials = 600  # Reduced from 2000 (10×600 = 6000 trials per point)

# Starting seed for P=0.2 dephasing (use 3000+ range)
seed_start = 3001

# Open output file
open("params_dephasing_p02.txt", "w") do f
    seed = seed_start
    
    for L in L_values
        for lambda in lambda_values
            lambda_x = lambda
            lambda_zz = round(1.0 - lambda, digits=2)
            
            for sample in 1:n_samples
                out_prefix = @sprintf("dephasing_p02_L%d_lam%.2f_P%.1f_s%d", L, lambda, P_x, sample)
                
                # Format: L P_x lambda_x lambda_zz lambda ntrials seed sample out_prefix
                println(f, "$L $P_x $lambda_x $lambda_zz $lambda $ntrials $seed $sample $out_prefix")
                
                seed += 1
            end
        end
    end
end

# Count jobs
n_jobs = length(L_values) * length(lambda_values) * n_samples
println("Generated $n_jobs jobs (CORRECTED VERSION):")
println("  L values: $(L_values)")
println("  λ values: $(length(lambda_values)) points")
println("  Samples: $n_samples")
println("  Trials per job: $ntrials")
println("  Total trials per (L,λ): $(n_samples * ntrials)")
println("  Seeds: $seed_start to $(seed_start + n_jobs - 1)")
println("")
println("CORRECTED: Proper quantum channel implementation (1-p)ρ + p·X·ρ·X")
