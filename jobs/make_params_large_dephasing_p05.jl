#!/usr/bin/env julia

"""
Generate parameter file for Large Dephasing P=0.5 jobs (L=20,24,28,32,36)
"""

using Printf

# System sizes for large systems (L=20 already complete, only need L=24,28,32,36)
L_values = [24, 28, 32, 36]

# Lambda values (same as small systems)
lambda_values = [0.1, 0.2, 0.3, 0.4, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.6, 0.7, 0.8, 0.9]

# Dephasing probability
P_x = 0.5

# Number of samples per configuration
# FAST VERSION: Even more samples with fewer trials for maximum parallelization
n_samples = 24  # Was 12, now 24 (even more parallel jobs!)

# Number of trials per job
ntrials = 250  # Was 500, now 250 (2× faster, same total: 24×250 = 6000 trials per point)

# Starting seed for large dephasing (use 7000+ range to avoid conflicts)
seed_start = 7001

# Open output file
open("params_large_dephasing_p05.txt", "w") do f
    seed = seed_start
    
    for L in L_values
        for lambda in lambda_values
            lambda_x = lambda
            lambda_zz = round(1.0 - lambda, digits=2)  # Fix floating-point precision
            
            for sample in 1:n_samples
                out_prefix = @sprintf("large_dephasing_p05_L%d_lam%.2f_P%.1f_s%d", L, lambda, P_x, sample)
                
                # Format: L P_x lambda_x lambda_zz lambda ntrials seed sample out_prefix
                println(f, "$L $P_x $lambda_x $lambda_zz $lambda $ntrials $seed $sample $out_prefix")
                
                seed += 1
            end
        end
    end
end

# Count jobs
n_jobs = length(L_values) * length(lambda_values) * n_samples
println("Generated $n_jobs jobs:")
println("  L values: $(L_values) (excluding L=20)")
println("  λ values: $(length(lambda_values)) points")
println("  Samples: $n_samples")
println("  Trials per job: $ntrials")
println("  Total trials per (L,λ): $(n_samples * ntrials)")
println("  Seeds: $seed_start to $(seed_start + n_jobs - 1)")
println("")
println("ULTRA-FAST strategy: 8× more jobs, each 8× faster, same total statistics!")
