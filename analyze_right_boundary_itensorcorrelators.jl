#!/usr/bin/env julia

"""
Analyze Right Boundary Scan Results (ITensorCorrelators.jl variant) - Plot Generation Only

Load all JSON files, compute statistics, and generate plots.
Properly groups by (L, lambda_x, lambda_zz, P_x) to avoid mixing different parameter sets.

Focus: Rényi-2 Binder (M2, M4 moments), computed via ITensorCorrelators.jl
n-point correlators instead of a global replica-overlap MPO.
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
            f -> contains(f, "lx0.70_lzz0.00_P") &&
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
    """Compute mean and std of B for each (L, lambda_x, lambda_zz, P_x) combination."""
    gdf = groupby(df, [:L, :lambda_x, :lambda_zz, :P_x])
    
    stats = combine(gdf) do group
        (
            B_mean = mean(group.B),
            B_std = std(group.B),
            B_sem = std(group.B) / sqrt(nrow(group)),
            M2_mean = mean(group.M2_bar),
            M4_mean = mean(group.M4_bar),
            purity_mean = mean(group.purity_bar),
            n_samples = nrow(group),
            time_mean = mean(group.time_seconds)
        )
    end
    
    sort!(stats, [:lambda_x, :lambda_zz, :L, :P_x])
    return stats
end

function plot_renyi2_binder_vs_px(stats)
    """Create plot of Rényi-2 Binder parameter vs P for λ_x=0.7, λ_zz=0."""
    
    # Should only have one parameter set now (lambda_x=0.7, lambda_zz=0.0)
    p1 = plot(xlabel=L"P_{x} = P_{zz}", ylabel="Binder Parameter (Rényi-2)", 
              title=L"\lambda_x = 0.7, \lambda_{zz} = 0.0" * " (ITensorCorrelators.jl)",
              legend=:outertopright, grid=true, size=(800, 600), dpi=300)
    
    colors = [:blue, :red, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle]
    
    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)
        
        plot!(p1, L_data.P_x, L_data.B_mean,
              yerr=L_data.B_sem,
              label="L = $L",
              color=colors[idx],
              marker=markers[idx],
              markersize=6,
              linewidth=2)
    end
    
    hline!(p1, [2/3], linestyle=:dash, color=:black, label="B = 2/3", linewidth=2, alpha=0.7)
    
    savefig(p1, "right_boundary_itensorcorrelators_lx0.7_lzz0.0_renyi2_binder_vs_px.pdf")
    savefig(p1, "right_boundary_itensorcorrelators_lx0.7_lzz0.0_renyi2_binder_vs_px.png")
end

function plot_moments_vs_px(stats)
    """Create plots of M2 and M4 moments vs P_x."""
    
    p2 = plot(xlabel=L"P_x" * " (X dephasing probability)", ylabel=L"M_2", 
              title=L"M_2 \text{ vs } P_x \text{ (λ_x=0.7, λ_zz=0.0)}",
              legend=:outertopright, grid=true, size=(800, 600), dpi=300)
    
    colors = [:blue, :red, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle]
    
    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)
        plot!(p2, L_data.P_x, L_data.M2_mean,
              label="L = $L",
              color=colors[idx],
              marker=markers[idx],
              markersize=6,
              linewidth=2)
    end
    
    savefig(p2, "right_boundary_itensorcorrelators_lx0.7_lzz0.0_M2_vs_px.pdf")
    savefig(p2, "right_boundary_itensorcorrelators_lx0.7_lzz0.0_M2_vs_px.png")
    
    # M4 plot
    p3 = plot(xlabel=L"P_x" * " (X dephasing probability)", ylabel=L"M_4", 
              title=L"M_4 \text{ vs } P_x \text{ (λ_x=0.7, λ_zz=0.0)}",
              legend=:outertopright, grid=true, size=(800, 600), dpi=300)
    
    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)
        plot!(p3, L_data.P_x, L_data.M4_mean,
              label="L = $L",
              color=colors[idx],
              marker=markers[idx],
              markersize=6,
              linewidth=2)
    end
    
    savefig(p3, "right_boundary_itensorcorrelators_lx0.7_lzz0.0_M4_vs_px.pdf")
    savefig(p3, "right_boundary_itensorcorrelators_lx0.7_lzz0.0_M4_vs_px.png")
end

function main()
    # Find the most recent results directory
    results_dirs = filter(d -> startswith(d, "right_boundary_itensorcorrelators_results_") && isdir(d), readdir())
    if isempty(results_dirs)
        println("ERROR: No results directory found matching 'right_boundary_itensorcorrelators_results_*'")
        return
    end
    
    results_dir = last(sort(results_dirs))  # Get most recent by timestamp
    println("Loading results from: $results_dir")
    
    df = load_results(results_dir)
    
    if nrow(df) == 0
        println("ERROR: No data loaded. Check that results directory contains lx0.70_lzz0.00_P*.json files")
        return
    end
    
    println("Loaded $(nrow(df)) data points")
    stats = compute_statistics(df)
    
    println("\nGenerating plots...")
    plot_renyi2_binder_vs_px(stats)
    plot_moments_vs_px(stats)
    
    println("✓ Plots generated:")
    println("  - right_boundary_itensorcorrelators_lx0.7_lzz0.0_renyi2_binder_vs_px.pdf")
    println("  - right_boundary_itensorcorrelators_lx0.7_lzz0.0_M2_vs_px.pdf")
    println("  - right_boundary_itensorcorrelators_lx0.7_lzz0.0_M4_vs_px.pdf")
end

main()
