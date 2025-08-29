
#!/usr/bin/env julia
using Pkg; Pkg.activate(".")

# Install packages if not available
required_packages = ["ITensors", "ITensorMPS", "ITensorCorrelators", "ArgParse", "JSON", "DelimitedFiles", "Random", "Statistics"]

println("Checking and installing required packages...")
for pkg in required_packages
    try
        eval(Meta.parse("using $pkg"))
        println("✓ $pkg already available")
    catch
        println("Installing $pkg...")
        Pkg.add(pkg)
        eval(Meta.parse("using $pkg"))
        println("✓ $pkg installed and loaded")
    end
end

# Add the src directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
using BinderSim

function build_parser()
    s = ArgParseSettings(; description = "Run Binder parameter calculation and save results as JSON.")
    @add_arg_table s begin
        "--L"
            help = "System size"
            arg_type = Int
            required = true
        "--lambda_x"
            help = "Parameter lambda_x"
            arg_type = Float64
            required = true
        "--lambda_zz"
            help = "Parameter lambda_zz" 
            arg_type = Float64
            required = true
        "--lambda"
            help = "Parameter lambda (alternative to lambda_x)"
            arg_type = Float64
            default = 0.1
        "--maxdim"
            help = "Maximum bond dimension"
            arg_type = Int
            default = 256
        "--cutoff"
            help = "SVD cutoff"
            arg_type = Float64
            default = 1e-12
        "--ntrials"
            help = "Number of Monte Carlo trials"
            arg_type = Int
            default = 1000
        "--chunk4"
            help = "Chunk size for S4 calculation"
            arg_type = Int
            default = 50000
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 1234
        "--sample"
            help = "Sample number (for output filename)"
            arg_type = Int
            default = 1
        "--out_prefix"
            help = "Output filename prefix"
            arg_type = String
            default = "result"
        "--outdir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end
    return s
end

function main(args)
    opts = parse_args(args, build_parser())
    
    # Extract parameters
    L = opts["L"]
    lambda_x = opts["lambda_x"]
    lambda_zz = opts["lambda_zz"]
    lambda = opts["lambda"]
    maxdim = opts["maxdim"]
    cutoff = opts["cutoff"]
    ntrials = opts["ntrials"]
    chunk4 = opts["chunk4"]
    seed = opts["seed"]
    sample = opts["sample"]
    out_prefix = opts["out_prefix"]
    outdir = opts["outdir"]
    
    println("Computing Binder parameter for L=$L, λₓ=$lambda_x, λ_zz=$lambda_zz, seed=$seed")
    
    # Run the actual simulation
    result = ea_binder_mc(L; lambda_x=lambda_x, lambda_zz=lambda_zz, 
                         ntrials=ntrials, maxdim=maxdim, cutoff=cutoff,
                         chunk4=chunk4, seed=seed)
    
    # Prepare output data
    output_data = Dict(
        "L" => L,
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "lambda" => lambda,
        "maxdim" => maxdim,
        "cutoff" => cutoff,
        "ntrials" => ntrials,
        "chunk4" => chunk4,
        "seed" => seed,
        "sample" => sample,
        "out_prefix" => out_prefix,
        "binder" => result.B,
        "binder_mean_of_trials" => result.B_mean_of_trials,
        "binder_std_of_trials" => result.B_std_of_trials,
        "S2_bar" => result.S2_bar,
        "S4_bar" => result.S4_bar,
        "ntrials_completed" => result.ntrials
    )
    
    println("Binder parameter computed: B = $(result.B)")
    
    # Save results
    mkpath(outdir)
    outfile = joinpath(outdir, "$(out_prefix).json")
    open(outfile, "w") do io
        JSON.print(io, output_data)
    end
    println("Wrote ", outfile)
end

main(ARGS)
