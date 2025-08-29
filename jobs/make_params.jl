#!/usr/bin/env julia

"Format one CLI line for run.jl."
format_line(L, lambda_x, lambda_zz, lambda, maxdim, cutoff, ntrials, chunk4, seed, sample, out_prefix, outdir) =
    "--L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --maxdim $maxdim --cutoff $cutoff --ntrials $ntrials --chunk4 $chunk4 --seed $seed --sample $sample --out_prefix $out_prefix --outdir $outdir"

"Write parameter lines to file."
function write_params(filename, lines)
    open(filename, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    println("Wrote $(length(lines)) lines to $filename")
end

# Generate TEST parameters (small for testing)
function generate_test_params()
    Ls = [8]  # One small system
    lambdas = [0.1, 0.5, 0.9]  # Three lambda values
    samples_per_lambda = 1
    ntrials = 100
    maxdim = 64
    chunk4 = 1000
    cutoff = 1e-12
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
        
        push!(lines, format_line(L, lambda_x, lambda_zz, λ, maxdim, cutoff, ntrials, chunk4, seed, sample, out_prefix, outdir))
    end
    
    write_params("jobs/params_test.txt", lines)
end

# Generate PRODUCTION parameters (full simulation per README)
function generate_production_params()
    Ls = [12, 16, 20, 24, 28]  # System sizes from README
    
    # Lambda values as specified in README
    coarse_lambdas = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
    fine_lambdas = [0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54]
    lambdas = sort(unique(vcat(coarse_lambdas, fine_lambdas)))
    
    samples_per_lambda = 3  # 3 independent runs
    ntrials = 1000  # Monte Carlo trials
    maxdim = 256
    chunk4 = 50000
    cutoff = 1e-12
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
        
        push!(lines, format_line(L, lambda_x, lambda_zz, λ, maxdim, cutoff, ntrials, chunk4, seed, sample, out_prefix, outdir))
    end
    
    write_params("jobs/params_production.txt", lines)
end

# Generate both files
println("Generating test parameters...")
generate_test_params()

println("Generating production parameters...")
generate_production_params()

println("Done! Use:")
println("  jobs/params_test.txt - for testing (3 jobs)")
println("  jobs/params_production.txt - for full simulation (255 jobs)")
