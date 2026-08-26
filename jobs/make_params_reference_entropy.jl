#!/usr/bin/env julia

"""
Generate parameter file for the reference-entropy diagnostic.
Follows the same conventions as make_params_right_boundary.jl /
make_params_defect_insertion.jl: L in {8,16,24,32}, q scan 0-0.5,
40 samples x 100 trials per (L,q) for parallelization, T_max = 2L
(the right_boundary baseline; left as reference_entropy_run.jl's
--T_max_factor default rather than a params column).
"""

using Printf

L_values = [8, 16, 24, 32]
q_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
n_samples = 40
ntrials = 100
seed_start = 91001

let next_seed = seed_start
    open("params_reference_entropy.txt", "w") do io
        for L in L_values
            for q in q_values
                for sample in 1:n_samples
                    out_prefix = @sprintf(
                        "reference_entropy_L%d_lx%.2f_lzz%.2f_q%.2f_s%d",
                        L, 0.7, 0.0, q, sample,
                    )
                    println(io, "$L 0.7 0.0 $q $ntrials $next_seed $sample $out_prefix")
                    next_seed += 1
                end
            end
        end
    end
end

println("Parameter file generated: params_reference_entropy.txt")
println("L values: ", L_values)
println("T_max: 2L (reference_entropy_run.jl default --T_max_factor)")
println("Total jobs: ", length(L_values) * length(q_values) * n_samples)
