#!/usr/bin/env julia

"""
Reference-entropy diagnostic (cluster/local entry point). See
reference_entropy_core.jl for the full derivation. Writes one raw CSV
row per (trial, saved time). Never averages density matrices before
computing entropy -- see the "averaging test" in
reference_entropy_validate.jl for the numerical demonstration of why
that would be wrong.

Primary scan (Section 4 of the request): lambda_x = 0.7, lambda_zz = 0
(the right-boundary slice), q_x = q_zz = q in [0, 0.5].

General-lambda scan (same slice as `lambda_scan_susceptibilities_itensorcorrelators.jl`,
NOT invented here): pass --lambda to set lambda_x = delta*lambda,
lambda_zz = delta*(1-lambda) with delta = 0.7, or pass --lambda_x /
--lambda_zz directly.
"""

using ArgParse
using Dates
using Random

include("reference_entropy_core.jl")

const REFERENCE_ENTROPY_DELTA = 0.7  # matches lambda_scan_susceptibilities_itensorcorrelators.jl's documented mapping

function parse_args_reference_entropy()
    settings = ArgParseSettings(description="Reference-entropy diagnostic")
    @add_arg_table! settings begin
        "--L"; arg_type=Int; default=16
        "--lambda_x"; arg_type=Float64; default=0.7
        "--lambda_zz"; arg_type=Float64; default=0.0
        "--lambda"; arg_type=Float64; default=NaN  # if set, overrides lambda_x/lambda_zz via delta*lambda mapping
        "--q"; arg_type=Float64; default=0.3
        "--ntrials"; arg_type=Int; default=50
        "--maxdim"; arg_type=Int; default=256
        "--cutoff"; arg_type=Float64; default=1e-12
        "--T_max_factor"; arg_type=Int; default=4
        "--save_every"; arg_type=Int; default=0  # 0 => auto: max(1, L/8)
        "--seed"; arg_type=Int; default=42
        "--output_dir"; arg_type=String; default="reference_entropy_results"
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
    "L", "q", "lambda_x", "lambda_zz", "q_x", "q_zz",
    "trial", "seed", "time", "T_max",
    "rhoR_00_real", "rhoR_00_imag",
    "rhoR_01_real", "rhoR_01_imag",
    "rhoR_10_real", "rhoR_10_imag",
    "rhoR_11_real", "rhoR_11_imag",
    "eig1", "eig2", "reference_entropy",
    "trace_error", "hermiticity_error",
    "minimum_eigenvalue", "max_linkdim",
    "record_checksum",
]

function main()
    args = parse_args_reference_entropy()
    L = args["L"]
    @assert L >= 2
    @assert args["ntrials"] >= 1
    @assert 0.0 <= args["q"] <= 0.5

    lambda_x = args["lambda_x"]
    lambda_zz = args["lambda_zz"]
    if !isnan(args["lambda"])
        lam = args["lambda"]
        @assert 0.0 <= lam <= 1.0
        lambda_x = REFERENCE_ENTROPY_DELTA * lam
        lambda_zz = REFERENCE_ENTROPY_DELTA * (1.0 - lam)
    end

    T_max = args["T_max_factor"] * L
    save_every = args["save_every"] > 0 ? args["save_every"] : max(1, L ÷ 8)
    mkpath(args["output_dir"])

    @info "Reference entropy S_R(m) is computed per Born-sampled trajectory " *
          "and then averaged over trajectories (sample mean), never by " *
          "computing the entropy of an averaged reference density matrix. " *
          "lambda_x=$lambda_x lambda_zz=$lambda_zz q_x=q_zz=$(args["q"]) T_max=$T_max"

    rows = Vector{Any}[]
    master_rng = MersenneTwister(args["seed"])

    for trial in 1:args["ntrials"]
        trial_seed = Int(rand(master_rng, UInt32))
        evolved = reference_entropy_evolve_one_trial(
            L; lambda_x=lambda_x, lambda_zz=lambda_zz,
            q_x=args["q"], q_zz=args["q"], T_max=T_max,
            maxdim=args["maxdim"], cutoff=args["cutoff"],
            seed=trial_seed, save_every=save_every,
        )

        for snap in evolved.snapshots
            M = snap.rhoR
            e = snap.entropy
            push!(rows, Any[
                L, args["q"], lambda_x, lambda_zz, args["q"], args["q"],
                trial, trial_seed, snap.time, T_max,
                real(M[1, 1]), imag(M[1, 1]),
                real(M[1, 2]), imag(M[1, 2]),
                real(M[2, 1]), imag(M[2, 1]),
                real(M[2, 2]), imag(M[2, 2]),
                e.eig1, e.eig2, e.S,
                e.trace_error, e.hermiticity_error,
                e.minimum_eigenvalue, evolved.max_interphysical_linkdim,
                evolved.record_checksum,
            ])
        end
    end

    filename = isempty(args["output_file"]) ?
        "reference_entropy_L$(L)_q$(args["q"])_$(Dates.format(now(), "yyyymmdd_HHMMSS")).csv" :
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
