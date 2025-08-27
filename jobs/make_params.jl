#!/usr/bin/env julia
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using ArgParse, JSON

function build_parser()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--Ls"
            default = "12,16,20,24,28"
        "--lmin"
            arg_type = Float64
            default = 0.1
        "--lmax"
            arg_type = Float64
            default = 0.9
        "--lstep"
            arg_type = Float64
            default = 0.1
        "--lstep_fine"
            arg_type = Float64
            default = 0.02
        "--fine_range_center"
            arg_type = Float64
            default = 0.5
        "--fine_range_width"
            arg_type = Float64
            default = 0.2
        "--samples"
            arg_type = Int
            default = 3
        "--ntrials"
            arg_type = Int
            default = 1000
        "--maxdim"
            arg_type = Int
            default = 256
        "--cutoff"
            arg_type = Float64
            default = 1e-12
        "--chunk4"
            arg_type = Int
            default = 50_000
        "--seed"
            arg_type = Int
            default = 1234
        "--out"
            default = "jobs/params.txt"
    end
    return s
end

function generate_lambda_values(lmin, lmax, lstep, lstep_fine, fine_center, fine_width)
    # Generate coarse grid
    coarse_lambdas = collect(lmin:lstep:lmax)
    
    # Generate fine grid around the center
    fine_min = max(lmin, fine_center - fine_width/2)
    fine_max = min(lmax, fine_center + fine_width/2)
    fine_lambdas = collect(fine_min:lstep_fine:fine_max)
    
    # Combine and remove duplicates, then sort
    all_lambdas = unique(vcat(coarse_lambdas, fine_lambdas))
    sort!(all_lambdas)
    
    return all_lambdas
end

function main()
    args = parse_args(build_parser())
    Ls = parse.(Int, split(args["Ls"], ","))
    lambdas = generate_lambda_values(args["lmin"], args["lmax"], args["lstep"], 
                                   args["lstep_fine"], args["fine_range_center"], 
                                   args["fine_range_width"])
    samples = args["samples"]
    ntrials = args["ntrials"]
    maxdim = args["maxdim"]
    cutoff = args["cutoff"]
    chunk4 = args["chunk4"]
    seed0 = args["seed"]
    out = args["out"]
    
    println("Generated λ values: ", lambdas)
    println("System sizes L: ", Ls)
    println("Trials per job: ", ntrials)
    println("Samples per (L,λ): ", samples)
    
    open(out, "w") do io
        idx = 0
        for L in Ls, λ in lambdas, s in 1:samples
            idx += 1
            params = Dict(
                "L"=>L,
                "lambda_x"=>round(λ, digits=6),
                "lambda_zz"=>round(1.0-λ, digits=6),
                "lambda"=>round(λ, digits=6),
                "sample"=>s,
                "ntrials"=>ntrials,
                "maxdim"=>maxdim,
                "cutoff"=>cutoff,
                "chunk4"=>chunk4,
                "seed"=>seed0 + idx,
                "out_prefix"=>"L$(L)_lam$(round(λ, digits=3))_s$(s)",
            )
            println(io, JSON.json(params))
        end
    end
    @info "Wrote $(length(Ls) * length(lambdas) * samples) parameter sets to $out"
    @info "Total computational cost: $(length(Ls) * length(lambdas) * samples * ntrials) Monte Carlo trials"
end

isinteractive() || main()
