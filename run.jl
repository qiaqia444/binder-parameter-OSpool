
#!/usr/bin/env julia
using ArgParse, JSON, DelimitedFiles

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--params"
            help = "JSON object with parameters"
            default = ""
        "--params-file"
            help = "File with one JSON object per line"
            default = ""
        "--outdir"
            help = "Directory to write outputs"
            default = "output"
    end
    return parse_args(s)
end

function compute!(p::Dict)
    # TODO: replace this placeholder with your simulation.
    # for now we just compute a fake binder as a function of lambda
    λx = p["lambda_x"]; λzz = p["lambda_zz"]; L = p["L"]
    binder = 2/3 + 0.3*tanh( (λx-λzz) * L / 50 )
    p["binder"] = binder
    return p
end

function run_one(json_str::AbstractString, outdir::AbstractString)
    p = JSON.parse(json_str)
    res = compute!(p)
    mkpath(outdir)
    outfile = joinpath(outdir, string(p["out_prefix"], ".json"))
    open(outfile, "w") do io
        JSON.print(io, res)
    end
    println("Wrote ", outfile)
end

function main()
    args = parse_cli()
    if args["params"] != ""
        run_one(args["params"], args["outdir"])
    elseif args["params_file"] != ""
        for line in eachline(args["params_file"])
            run_one(line, args["outdir"])
        end
    else
        error("Provide --params JSON or --params-file")
    end
end

isinteractive() || main()
