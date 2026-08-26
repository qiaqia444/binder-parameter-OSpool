#!/usr/bin/env julia

using Printf

L_values = [8, 16, 24, 32]
q_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
ntrials = 50
T_max_factor = 4
seed_start = 91001

let next_seed = seed_start
    open("params_reference_entropy.txt", "w") do io
        for L in L_values
            for q in q_values
                out_prefix = @sprintf(
                    "reference_entropy_L%d_lx0.70_lzz0.00_q%.2f",
                    L, q,
                )
                println(io, "$L 0.7 0.0 $q $ntrials $T_max_factor $next_seed $out_prefix")
                next_seed += 1
            end
        end
    end
end

println("Parameter file generated: params_reference_entropy.txt")
println("L values: ", L_values)
println("T_max_factor: ", T_max_factor, " (T_max = T_max_factor * L)")
println("Total jobs: ", length(L_values) * length(q_values))
