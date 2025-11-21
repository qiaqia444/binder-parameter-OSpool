#!/usr/bin/env julia

"""
Generate parameter file for RBIM simulations
Random Bond Ising Model in strong dephasing limit (λ_x = 0)
"""

using Printf

# System sizes (start with small systems for quick results)
L_values = [8, 12, 16, 20]

# ZZ measurement strengths (fix a few representative values)
lambda_zz_values = [0.3, 0.5, 0.7]

# Dephasing probabilities to scan (looking for Nishimori point)
# Dense sampling around expected critical region
P_x_values = [0.05, 0.08, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.17, 0.20, 0.25, 0.30, 0.40, 0.50]

# Number of samples per configuration
n_samples = 10

# Number of trials per job
ntrials = 600

# Starting seed
seed_start = 9001

# Open output file
open("params_rbim.txt", "w") do f
    seed = seed_start
    
    for P_x in P_x_values
        for L in L_values
            for lambda_zz in lambda_zz_values
                for sample in 1:n_samples
                    out_prefix = @sprintf("rbim_L%d_lzz%.2f_P%.1f_s%d", L, lambda_zz, P_x, sample)
                    
                    # Format: L P_x lambda_zz ntrials seed sample out_prefix
                    println(f, "$L $P_x $lambda_zz $ntrials $seed $sample $out_prefix")
                    
                    seed += 1
                end
            end
        end
    end
end

# Count jobs
n_jobs = length(P_x_values) * length(L_values) * length(lambda_zz_values) * n_samples
println("Generated $n_jobs jobs for RBIM:")
println("  L values: $(L_values)")
println("  λ_zz values: $(length(lambda_zz_values)) points")
println("  P_x values: $(P_x_values)")
println("  Samples: $n_samples")
println("  Trials per job: $ntrials")
println("  Total trials per (L,λ_zz,P): $(n_samples * ntrials)")
println("  Seeds: $seed_start to $(seed_start + n_jobs - 1)")
println("")
println("RBIM mapping: λ_x = 0 (strong dephasing limit)")
println("Order parameter: EA correlations [⟨Z_i Z_j⟩²]_J")
