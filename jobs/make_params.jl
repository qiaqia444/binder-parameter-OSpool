#!/usr/bin/env julia
using ArgParse, JSON

function build_parser()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--Ls"
            default = "8,12,16"
        "--lmin"
            arg_type = Float64
            default = 0.0
        "--lmax"
            arg_type = Float64
            default = 1.0
        "--lstep"
            arg_type = Float64
            default = 0.1
        "--samples"
            arg_type = Int
            default = 3
        "--seed"
            arg_type = Int
            default = 1234
        "--out"
            default = "jobs/params.txt"
    end
    return s
end

function frange(a,b,s)
    n = Int(round((b-a)/s))
    [a + i*s for i in 0:n]
end

function main()
    args = parse_args(build_parser())
    Ls = parse.(Int, split(args["Ls"], ","))
    lambdas = frange(args["lmin"], args["lmax"], args["lstep"])
    samples = args["samples"]
    seed0 = args["seed"]
    out = args["out"]
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
                "seed"=>seed0 + idx,
                "out_prefix"=>"L$(L)_lam$(round(λ, digits=3))_s$(s)",
            )
            println(io, JSON.json(params))
        end
    end
    @info "Wrote params to $out"
end

isinteractive() || main()
