#!/usr/bin/env julia

"""
Right-boundary Renyi-2 disorder-string overlap diagnostic (cluster/local
entry point). See defect_insertion_core.jl for the full derivation and
the physical-limitation caveats. Writes one raw CSV row per trajectory,
endpoint pair, and separation r -- no trajectory-level averaging or
NaN-discarding happens here (that belongs to defect_insertion_analyze.jl,
which must average translated pairs within a trajectory before computing
error bars across independent trajectories).
"""

using ArgParse
using Dates
using Random

include("defect_insertion_core.jl")

function parse_args_defect()
    settings = ArgParseSettings(description="Right-boundary Renyi-2 disorder-string overlap diagnostic")
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

function csv_field(x::AbstractFloat)
    isnan(x) && return "NaN"
    isinf(x) && return x > 0 ? "Inf" : "-Inf"
    return string(x)
end
csv_field(x) = string(x)

const RAW_CSV_HEADER = [
    "L", "q", "trial", "seed", "i", "j", "r", "z0", "zij", "R", "absR", "signR",
    "string_log_magnitude", "is_numerical_zero",
    "contraction_error", "record_checksum", "max_linkdim", "max_trace_error",
]

function main()
    args = parse_args_defect()
    L = args["L"]
    @assert haskey(DISORDER_STRING_PAIRS, L) "L must be 16, 24, or 32."
    @assert isapprox(args["lambda_x"], 0.7; atol=1e-9)
    @assert isapprox(args["lambda_zz"], 0.0; atol=1e-9)
    @assert args["ntrials"] >= 1
    @assert 0.0 <= args["P"] <= 0.5
    mkpath(args["output_dir"])

    @info "R_ij = zij/z0 is a SIGNED Renyi-2 disorder-string overlap " *
          "(see defect_insertion_core.jl); it is NOT a partition-function " *
          "ratio unless separately proven. No trajectory is discarded; " *
          "R < 0 and R == 0 are recorded, not clipped."

    rows = Vector{Any}[]
    master_rng = MersenneTwister(args["seed"])
    r_values = sort(collect(keys(DISORDER_STRING_PAIRS[L])))

    for trial in 1:args["ntrials"]
        trial_seed = Int(rand(master_rng, UInt32))
        evolved = defect_insertion_evolve_one_trial(
            L; lambda_x=args["lambda_x"], q=args["P"], T_max=L,
            maxdim=args["maxdim"], cutoff=args["cutoff"], seed=trial_seed,
        )

        for r in r_values
            for (i, j) in DISORDER_STRING_PAIRS[L][r]
                value = defect_insertion_contract_pair(
                    evolved.rho, evolved.sites, i, j;
                    maxdim=args["maxdim"], cutoff=args["cutoff"],
                )
                push!(rows, Any[
                    L, args["P"], trial, trial_seed, i, j, value.r,
                    value.z0, value.zij, value.R, value.absR, value.signR,
                    value.string_log_magnitude, value.is_numerical_zero,
                    value.contraction_error, evolved.record_checksum,
                    evolved.max_interphysical_linkdim, evolved.max_trace_error,
                ])
            end
        end
    end

    filename = isempty(args["output_file"]) ?
        "defect_insertion_L$(L)_P$(args["P"])_$(Dates.format(now(), "yyyymmdd_HHMMSS")).csv" :
        args["output_file"]
    outpath = joinpath(args["output_dir"], filename)

    open(outpath, "w") do io
        println(io, join(RAW_CSV_HEADER, ","))
        for row in rows
            println(io, join(csv_field.(row), ","))
        end
    end

    println("Saved $(length(rows)) rows to $outpath")
end

main()
