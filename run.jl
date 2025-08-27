
#!/usr/bin/env julia
using ArgParse, JSON, DelimitedFiles
using Pkg; Pkg.activate(".")
using BinderSim

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
    # Extract parameters
    L = p["L"]
    lambda_x = p["lambda_x"]
    lambda_zz = p["lambda_zz"]
    seed = get(p, "seed", 1234)
    ntrials = get(p, "ntrials", 100)
    maxdim = get(p, "maxdim", 256)
    cutoff = get(p, "cutoff", 1e-12)
    chunk4 = get(p, "chunk4", 50_000)
    
    println("Computing Binder parameter for L=$L, λₓ=$lambda_x, λ_zz=$lambda_zz, seed=$seed")
    
    # Run the actual simulation
    result = ea_binder_mc(L; lambda_x=lambda_x, lambda_zz=lambda_zz, 
                         ntrials=ntrials, maxdim=maxdim, cutoff=cutoff,
                         chunk4=chunk4, seed=seed)
    
    # Store results in the parameter dictionary
    p["binder"] = result.B
    p["binder_mean_of_trials"] = result.B_mean_of_trials
    p["binder_std_of_trials"] = result.B_std_of_trials
    p["S2_bar"] = result.S2_bar
    p["S4_bar"] = result.S4_bar
    p["ntrials_completed"] = result.ntrials
    
    println("Binder parameter computed: B = $(result.B)")
    
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
