#!/usr/bin/env julia

"""
Analyze Lambda-Path Susceptibility Scan Results (ITensorCorrelators.jl
variant) - Plot Generation Only

Load all JSON files, compute statistics, and generate plots.
Properly groups by (L, delta, q, lambda) to avoid mixing different
parameter sets.

Observables: Rényi-2 replica-overlap susceptibility kappa2 (ITensorCorrelators.jl)
and Edwards-Anderson susceptibility kappaEA, both vs the path parameter lambda
(lambda_x = delta*lambda, lambda_zz = delta*(1-lambda)).
"""

using JSON
using Statistics
using DataFrames
using Printf
using Plots
using LaTeXStrings

function load_results(results_dir)
    """Load JSON files and tolerate files with different schemas."""

    rows = Dict{String,Any}[]

    for L in [8, 16, 24, 32]
        L_dir = joinpath(results_dir, "L$L")

        if !isdir(L_dir)
            continue
        end

        json_files = filter(
            f -> contains(f, "delta0.70_q0.10") &&
                 endswith(f, ".json"),
            readdir(L_dir, join=true),
        )

        for file in json_files
            try
                parsed = JSON.parsefile(file)
                entries = parsed isa AbstractVector ? parsed : [parsed]

                for entry in entries
                    if !(entry isa AbstractDict)
                        @warn "Skipping non-dictionary JSON entry" file
                        continue
                    end

                    row = Dict{String,Any}(
                        String(key) => value
                        for (key, value) in entry
                    )

                    # Helpful for finding where every row came from.
                    row["_source_file"] = file

                    push!(rows, row)
                end

            catch e
                @warn(
                    "Could not read JSON file",
                    file=file,
                    exception=(e, catch_backtrace()),
                )
            end
        end
    end

    isempty(rows) && return DataFrame()

    # Construct the union of all keys.
    all_keys = Set{String}()

    for row in rows
        union!(all_keys, keys(row))
    end

    # Insert `missing` wherever an older JSON lacks a newer field.
    for row in rows
        for key in all_keys
            if !haskey(row, key)
                row[key] = missing
            end
        end
    end

    return DataFrame(rows)
end

function compute_statistics(df)
    """Compute mean and std of kappa2/kappaEA for each (L, delta, q, lambda) combination."""
    gdf = groupby(df, [:L, :delta, :q, :lambda])

    stats = combine(gdf) do group
        (
            kappa2_mean = mean(group.kappa2),
            kappa2_std = std(group.kappa2),
            kappa2_sem = std(group.kappa2) / sqrt(nrow(group)),
            kappaEA_mean = mean(group.kappaEA),
            kappaEA_std = std(group.kappaEA),
            kappaEA_sem = std(group.kappaEA) / sqrt(nrow(group)),
            purity_mean = mean(group.purity),
            n_samples = nrow(group),
            time_mean = mean(group.time_seconds)
        )
    end

    sort!(stats, [:delta, :q, :L, :lambda])
    return stats
end

function plot_kappa2_vs_lambda(stats)
    """Create plot of Rényi-2 overlap susceptibility vs lambda."""

    p1 = plot(xlabel=L"\lambda", ylabel=L"\kappa_2",
              title="Rényi-2 overlap susceptibility (ITensorCorrelators.jl)",
              legend=:outertopright, grid=true, size=(800, 600), dpi=300)

    colors = [:blue, :red, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle]

    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)

        plot!(p1, L_data.lambda, L_data.kappa2_mean,
              yerr=L_data.kappa2_sem,
              label="L = $L",
              color=colors[idx],
              marker=markers[idx],
              markersize=6,
              linewidth=2)
    end

    savefig(p1, "lambda_scan_susceptibility_delta0.7_q0.1_kappa2_vs_lambda.pdf")
    savefig(p1, "lambda_scan_susceptibility_delta0.7_q0.1_kappa2_vs_lambda.png")
end

function plot_kappaEA_vs_lambda(stats)
    """Create plot of Edwards-Anderson susceptibility vs lambda."""

    p2 = plot(xlabel=L"\lambda", ylabel=L"\kappa_{EA}",
              title="Edwards-Anderson susceptibility",
              legend=:outertopright, grid=true, size=(800, 600), dpi=300)

    colors = [:blue, :red, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle]

    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)
        plot!(p2, L_data.lambda, L_data.kappaEA_mean,
              yerr=L_data.kappaEA_sem,
              label="L = $L",
              color=colors[idx],
              marker=markers[idx],
              markersize=6,
              linewidth=2)
    end

    savefig(p2, "lambda_scan_susceptibility_delta0.7_q0.1_kappaEA_vs_lambda.pdf")
    savefig(p2, "lambda_scan_susceptibility_delta0.7_q0.1_kappaEA_vs_lambda.png")
end

function main()
    # Find the most recent results directory
    results_dirs = filter(d -> startswith(d, "lambda_scan_susceptibility_results_") && isdir(d), readdir())
    if isempty(results_dirs)
        println("ERROR: No results directory found matching 'lambda_scan_susceptibility_results_*'")
        return
    end

    results_dir = last(sort(results_dirs))  # Get most recent by timestamp
    println("Loading results from: $results_dir")

    df = load_results(results_dir)

    if nrow(df) == 0
        println("ERROR: No data loaded. Check that results directory contains delta0.70_q0.10_*.json files")
        return
    end

    println("Loaded $(nrow(df)) data points")
    stats = compute_statistics(df)

    println("\nGenerating plots...")
    plot_kappa2_vs_lambda(stats)
    plot_kappaEA_vs_lambda(stats)

    println("✓ Plots generated:")
    println("  - lambda_scan_susceptibility_delta0.7_q0.1_kappa2_vs_lambda.pdf")
    println("  - lambda_scan_susceptibility_delta0.7_q0.1_kappaEA_vs_lambda.pdf")
end

main()
