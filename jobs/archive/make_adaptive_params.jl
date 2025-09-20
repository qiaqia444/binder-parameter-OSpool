#!/usr/bin/env julia

"Format one CLI line for run_adaptive.jl with adaptive parameters."
function format_adaptive_line(L, lambda_x, lambda_zz, lambda, seed, sample, out_prefix, outdir)
    # Use adaptive parameters - the script will automatically adjust maxdim, cutoff, chunk4
    ntrials = 2000  # Your requested number of trials
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

# Generate ADAPTIVE production parameters (all system sizes with adaptive settings)
function generate_adaptive_production_params()
    Ls = [8, 12, 16]  # Your requested system sizes
    
    # Your specified lambda values
    coarse_lambdas = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
    fine_lambdas = [0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54]
    lambdas = sort(unique(vcat(coarse_lambdas, fine_lambdas)))
    
    samples_per_lambda = 3  # 3 independent runs
    base_seed = 1234
    outdir = "output"
    
    lines = String[]
    job_id = 0
    for L in Ls, λ in lambdas, sample in 1:samples_per_lambda
        job_id += 1
        lambda_x = round(λ, digits=6)
        lambda_zz = round(1.0 - λ, digits=6)
        seed = base_seed + job_id
        out_prefix = "L$(L)_lam$(round(λ, digits=3))_s$(sample)"
        
        push!(lines, format_adaptive_line(L, lambda_x, lambda_zz, λ, seed, sample, out_prefix, outdir))
    end
    
    write_params("jobs/params_adaptive.txt", lines)
end

# Generate test with your specified systems
function generate_adaptive_test_params()
    Ls = [8, 12, 16]  # Your requested system sizes
    lambdas = [0.3, 0.5, 0.7]  # Three representative lambda values
    samples_per_lambda = 1
    ntrials = 200  # Smaller for testing
    base_seed = 1234
    outdir = "output"
    
    lines = String[]
    job_id = 0
    for L in Ls, λ in lambdas, sample in 1:samples_per_lambda
        job_id += 1
        lambda_x = round(λ, digits=6)
        lambda_zz = round(1.0 - λ, digits=6)
        seed = base_seed + job_id
        out_prefix = "adaptive_test_L$(L)_lam$(round(λ, digits=3))_s$(sample)"
        
        # For testing, use the adaptive format but override ntrials
        line = "--L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $λ --ntrials 200 --seed $seed --sample $sample --out_prefix $out_prefix --outdir $outdir"
        push!(lines, line)
    end
    
    write_params("jobs/params_adaptive_test.txt", lines)
end

println("Generating adaptive parameters...")
generate_adaptive_test_params()
generate_adaptive_production_params()

println("Done! Files created:")
println("  jobs/params_adaptive_test.txt - for testing larger systems (9 jobs)")
println("  jobs/params_adaptive.txt - for full adaptive simulation (255 jobs)")
println()
println("The adaptive approach automatically adjusts:")
println("  L ≤ 12: maxdim=256, cutoff=1e-12, chunk4=50000")
println("  L ≤ 16: maxdim=128, cutoff=1e-10, chunk4=20000")
println("  L ≥ 20: maxdim=64,  cutoff=1e-8,  chunk4=10000")
