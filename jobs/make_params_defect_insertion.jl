#!/usr/bin/env julia

using Printf

L_values = [16, 24, 32]
P_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
n_samples = 40
ntrials = 100
seed_start = 81001

let next_seed = seed_start
    open("params_defect_insertion.txt", "w") do io
        for L in L_values
            for P in P_values
                for sample in 1:n_samples
                    out_prefix = @sprintf(
                        "defect_insertion_L%d_lx%.2f_lzz%.2f_Px%.2f_s%d",
                        L, 0.7, 0.0, P, sample,
                    )
                    println(io, "$L 0.7 0.0 $P $P $ntrials $next_seed $sample $out_prefix")
                    next_seed += 1
                end
            end
        end
    end
end

println("Parameter file generated: params_defect_insertion.txt")
println("L values: ", L_values)
println("T_max: L")
println("Total jobs: ", length(L_values) * length(P_values) * n_samples)
