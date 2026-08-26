#!/usr/bin/env julia

using ArgParse
using JSON
using Dates
using Random

include("defect_insertion_core.jl")

function parse_args_defect()
    settings = ArgParseSettings(description="Right-boundary defect insertion diagnostic")
    @add_arg_table! settings begin
        "--L"; arg_type=Int; default=16
        "--lambda_x"; arg_type=Float64; default=0.7
        "--lambda_zz"; arg_type=Float64; default=0.0
        "--P"; arg_type=Float64; default=0.3
        "--ntrials"; arg_type=Int; default=100
        "--maxdim"; arg_type=Int; default=256
        "--cutoff"; arg_type=Float64; default=1e-12
        "--seed"; arg_type=Int; default=42
        "--output_dir"; arg_type=String; default="defect_insertion_results"
        "--output_file"; arg_type=String; default=""
    end
    return parse_args(settings)
end

function main()
    args = parse_args_defect()
    L = args["L"]
    @assert haskey(DEFECT_INSERTION_PAIRS, L) "L must be 16, 24, or 32."
    @assert isapprox(args["lambda_x"], 0.7; atol=1e-9)
    @assert isapprox(args["lambda_zz"], 0.0; atol=1e-9)
    @assert args["ntrials"] >= 1
    @assert 0.0 <= args["P"] <= 0.5
    mkpath(args["output_dir"])

    @info "Domain-wall defect D_ij = prod_{p=i}^{j-1} X_p^bra (see " *
          "defect_insertion_core.jl docstring). Z_ij can be negative for an " *
          "individual trajectory (it is a genuine, possibly-signed " *
          "correlator, not a manifestly positive Boltzmann weight); rows " *
          "with Z_ij <= 0 report deltaF = NaN and are excluded from the " *
          "trajectory average, per the quenched-averaging convention."

    rows = Any[]
    master_rng = MersenneTwister(args["seed"])

    for trial in 1:args["ntrials"]
        trial_seed = Int(rand(master_rng, UInt32))
        evolved = defect_insertion_evolve_one_trial(
            L; lambda_x=args["lambda_x"], q=args["P"], T_max=L,
            maxdim=args["maxdim"], cutoff=args["cutoff"], seed=trial_seed,
        )
        trajectory = defect_insertion_compute_trajectory(
            evolved, L, DEFECT_INSERTION_PAIRS[L];
            maxdim=args["maxdim"], cutoff=args["cutoff"],
        )
        for pair in trajectory.rows
            push!(rows, Dict(
                "row_type"=>"pair", "L"=>L, "T_max"=>L,
                "lambda_x"=>args["lambda_x"], "q"=>args["P"],
                "control_parameter"=>args["P"], "trial"=>trial,
                "seed"=>trial_seed, "record_id"=>trial,
                "record_checksum"=>pair.record_checksum, "i"=>pair.i,
                "j"=>pair.j, "r"=>pair.r, "logZ0"=>pair.logZ0,
                "logZij"=>pair.logZij, "logR"=>pair.logR,
                "deltaF"=>pair.deltaF, "maxdim"=>args["maxdim"],
                "cutoff"=>args["cutoff"], "contraction_error"=>pair.contraction_error,
            ))
        end
        push!(rows, Dict(
            "row_type"=>"trajectory_average", "L"=>L, "T_max"=>L,
            "lambda_x"=>args["lambda_x"], "q"=>args["P"],
            "control_parameter"=>args["P"], "trial"=>trial,
            "seed"=>trial_seed, "record_id"=>trial,
            "record_checksum"=>evolved.record_checksum, "deltaF"=>trajectory.average_deltaF,
            "maxdim"=>args["maxdim"], "cutoff"=>args["cutoff"],
        ))
    end

    filename = isempty(args["output_file"]) ?
        "defect_insertion_L$(L)_P$(args["P"])_$(Dates.format(now(), "yyyymmdd_HHMMSS")).json" : args["output_file"]
    open(joinpath(args["output_dir"], filename), "w") do io
        JSON.print(io, rows, 2)
    end
    println("Saved $(length(rows)) rows to $(joinpath(args["output_dir"], filename))")
end

main()
