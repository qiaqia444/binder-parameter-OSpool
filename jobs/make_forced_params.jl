#!/usr/bin/env julia

"Format one CLI line for run_forced.jl with adaptive parameters."
function format_forced_line(L, lambda_x, lambda_zz, lambda, seed, sample, out_prefix, outdir)
    # Use adaptive parameters - the script will automatically adjust maxdim, cutoff
    ntrials = 1  # Deterministic, so only need 1 trial per job
    return "--L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir $outdir"
end

"Write parameter lines to file."
function write_params(filename, lines)
    open(filename, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    println("Wrote $(length(lines)) lines to $filename")
end

# Generate FORCED +1 production parameters
function generate_forced_production_params()
    Ls = [8, 12, 16]  # System sizes that should work with forced measurements
    
    # Your specified lambda values
    coarse_lambdas = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
    fine_lambdas = [0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54]
    lambdas = sort(unique(vcat(coarse_lambdas, fine_lambdas)))
    
    samples_per_lambda = 3  # 3 independent runs for verification
    base_seed = 1234
    outdir = "output"
    
    lines = String[]
    job_id = 0
    for L in Ls, λ in lambdas, sample in 1:samples_per_lambda
        job_id += 1
        lambda_x = round(λ, digits=6)
        lambda_zz = round(1.0 - λ, digits=6)
        seed = base_seed + job_id
        out_prefix = "forced_L$(L)_lam$(round(λ, digits=3))_s$(sample)"
        
        push!(lines, format_forced_line(L, lambda_x, lambda_zz, λ, seed, sample, out_prefix, outdir))
    end
    
    write_params("jobs/params_forced.txt", lines)
end

# Generate test with smaller set
function generate_forced_test_params()
    Ls = [8, 12, 16]  # All system sizes for testing
    lambdas = [0.1, 0.3, 0.5, 0.7, 0.9]  # Representative lambda values
    samples_per_lambda = 1  # Only 1 since it's deterministic
    base_seed = 1234
    outdir = "output"
    
    lines = String[]
    job_id = 0
    for L in Ls, λ in lambdas, sample in 1:samples_per_lambda
        job_id += 1
        lambda_x = round(λ, digits=6)
        lambda_zz = round(1.0 - λ, digits=6)
        seed = base_seed + job_id
        out_prefix = "forced_test_L$(L)_lam$(round(λ, digits=3))_s$(sample)"
        
        push!(lines, format_forced_line(L, lambda_x, lambda_zz, λ, seed, sample, out_prefix, outdir))
    end
    
    write_params("jobs/params_forced_test.txt", lines)
end

println("Generating forced +1 measurement parameters...")
generate_forced_test_params()
generate_forced_production_params()

println("Done! Files created:")
println("  jobs/params_forced_test.txt - for testing (10 jobs)")
println("  jobs/params_forced.txt - for full forced +1 simulation (153 jobs)")
println()
println("The forced +1 approach uses:")
println("  • Deterministic evolution (all measurements forced to +1)")
println("  • Only 1 trial per job (since results are reproducible)")
println("  • Adaptive maxdim and cutoff based on system size")
println("  • Should work for larger systems due to deterministic nature")