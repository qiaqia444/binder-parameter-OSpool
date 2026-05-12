#!/usr/bin/env julia

using Printf

"""
Generate parameter file for EA Binder Memory-to-Trivial (λ_x=0.21, λ_zz=0.49) simulations.

Configuration:
- L = [8, 16, 24, 32]
- λ_x = 0.21 (measuring X)
- λ_zz = 0.49 (measuring ZZ)
- P values: 20 points with dense sampling around pc ≈ 0.05
- ntrials: 100 per job
- samples: 40 per L-P combination
- Total: 4 L × 20 P × 40 samples = 3,200 jobs

Output format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
"""

output_file = "params_memory_to_trivial_lx0.3.txt"

# Configuration
L_values = [8, 16, 24, 32]
lambda_x = 0.21
lambda_zz = 0.49
ntrials = 100
n_samples = 40

# P values optimized densely around pc ≈ 0.05
# Very fine sampling near critical point (20 points)
P_values = [
    0.0,      # Far below transition
    0.01,     # Approaching
    0.02,     # Approaching
    0.03,     # Close to pc
    0.04,     # Very close
    0.045,    # Fine scale
    0.05,     # At critical point
    0.055,    # Fine scale
    0.06,     # Just past
    0.07,     # Moving away
    0.08,     # Moving away
    0.10,     # Course scale
    0.15,     # Bulk region
    0.20,     # Higher regime
    0.25,     # Higher regime
    0.30,     # Higher regime
    0.35,     # Higher regime
    0.40,     # Higher regime
    0.45,     # Higher regime
    0.50      # Full dephasing
]

open(output_file, "w") do f
    seed_counter = 20001  # Start seeds at 20001 to avoid collision with left_boundary
    
    for L in L_values
        for P in P_values
            for sample in 1:n_samples
                # Format: L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix
                out_prefix = "memory_to_trivial_L$(L)_P$(Printf.format(Printf.Format("%.2f"), P))_s$(sample)"
                
                line = "$(L) $(lambda_x) $(lambda_zz) $(P) $(P) $(ntrials) $(seed_counter) $(sample) $(out_prefix)"
                println(f, line)
                
                seed_counter += 1
            end
        end
    end
end

println("✓ Generated $output_file")
println("  L values: $L_values")
println("  λ_x = $lambda_x, λ_zz = $lambda_zz")
println("  P values: $(length(P_values)) points (dense around pc ≈ 0.05)")
println("  Samples per L-P: $n_samples")
println("  Total jobs: $(length(L_values) * length(P_values) * n_samples)")
