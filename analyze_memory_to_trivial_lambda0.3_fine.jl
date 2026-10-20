#!/usr/bin/env julia

"""
Analyze Memory-to-Trivial EA Binder FINE-grid Scan Results
(lambda = 0.3: lambda_x = 0.21, lambda_zz = 0.49; P = 0.025:0.005:0.075) - Plot Generation Only

Usage:
    julia analyze_memory_to_trivial_lambda0.3_fine.jl [results_dir]

Follows analyze_memory_to_trivial_lx0.3.jl: loads all JSON files, computes statistics and
generates plots. Properly groups by (L, lambda_x, lambda_zz, P_x) to avoid mixing different
parameter sets. Only `memory_to_trivial_lambda0.3_fine_L*.json` files are read (default
directory: the newest `memory_to_trivial_lambda0.3_fine_results_*`).

Each job file holds one (L, P) point with its own EA Binder `B`; the statistics are the mean
and standard error of `B` over the jobs (samples) at each (L, P_x). Pairwise crossings of the
mean Binder curves are located by linear interpolation on the fine grid.

Outputs: memory_to_trivial_lambda0.3_fine_summary.csv,
         memory_to_trivial_lambda0.3_fine_binder_vs_p.{pdf,png},
         memory_to_trivial_lambda0.3_fine_L{L}_binder_vs_p.{pdf,png}
"""

using JSON
using Statistics
using DataFrames
using Printf
using Plots
using LaTeXStrings
using DelimitedFiles

const PREFIX = "memory_to_trivial_lambda0.3_fine"

function load_results(results_dir)
    """Load all JSON files from the results directory - memory_to_trivial_lambda0.3_fine_L*.json files."""
    all_data = []

    for L_dir in sort(filter(isdir, readdir(results_dir, join=true)))
        # Only load files matching the pattern: memory_to_trivial_lambda0.3_fine_L*.json
        json_files = filter(f -> startswith(basename(f), PREFIX * "_L") && endswith(f, ".json") &&
                                 !occursin("FAILED", f), readdir(L_dir, join=true))

        for file in json_files
            try
                data = JSON.parsefile(file)
                if data isa Array
                    for entry in data
                        push!(all_data, entry)
                    end
                else
                    push!(all_data, data)
                end
            catch e
                @warn "Skipping unreadable file" file
            end
        end
    end

    df = DataFrame(all_data)
    if nrow(df) > 0
        df.P_x = round.(df.P_x; digits=6)   # remove float noise from the 0.005 grid
    end
    return df
end

function compute_statistics(df)
    """Compute mean and std of B for each (L, lambda_x, lambda_zz, P_x) combination."""
    gdf = groupby(df, [:L, :lambda_x, :lambda_zz, :P_x])

    stats = combine(gdf) do group
        (
            B_mean = mean(group.B),
            B_std = std(group.B),
            B_sem = std(group.B) / sqrt(nrow(group)),
            M2_mean = mean(group.S2_bar),
            M4_mean = mean(group.S4_bar),
            n_samples = nrow(group),
            time_mean = mean(group.time_seconds)
        )
    end

    sort!(stats, [:lambda_x, :lambda_zz, :L, :P_x])
    return stats
end

function plot_binder_vs_px(stats)
    """Create plot of Binder parameter vs P for lambda_x = 0.21, lambda_zz = 0.49."""

    p1 = plot(xlabel=L"P_x = P_{zz}" * " (dephasing probability)", ylabel="Binder Parameter",
              title=L"\lambda = 0.3:\ \lambda_x = 0.21,\ \lambda_{zz} = 0.49",
              legend=:topright, grid=true, size=(800, 600), dpi=300)

    colors = [:blue, :red, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle]

    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)

        plot!(p1, L_data.P_x, L_data.B_mean,
              yerr=L_data.B_sem,
              label="L = $L",
              color=colors[mod1(idx, length(colors))],
              marker=markers[mod1(idx, length(markers))],
              markersize=6,
              linewidth=2)
    end

    savefig(p1, "$(PREFIX)_binder_vs_p.pdf")
    savefig(p1, "$(PREFIX)_binder_vs_p.png")
end

function plot_binder_by_L(stats)
    """Create separate plots for each L."""
    colors = [:steelblue, :darkred, :darkgreen, :darkorchid, :darkorange]

    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)

        p = plot(xlabel=L"P_x = P_{zz}" * " (dephasing probability)",
                ylabel="Binder Parameter",
                title=latexstring("\\lambda = 0.3:\\ \\lambda_x = 0.21,\\ \\lambda_{zz} = 0.49,\\ L = $L"),
                legend=false, grid=true, size=(800, 600), dpi=300)

        plot!(p, L_data.P_x, L_data.B_mean,
              yerr=L_data.B_sem,
              color=colors[mod1(idx, length(colors))],
              marker=:circle,
              markersize=8,
              linewidth=2.5,
              markerstrokewidth=1.5)

        savefig(p, "$(PREFIX)_L$(L)_binder_vs_p.pdf")
        savefig(p, "$(PREFIX)_L$(L)_binder_vs_p.png")
    end
end

# First sign change of d(P) on the grid, linearly interpolated; NaN if none.
function first_crossing(ps, d)
    for i in 1:(length(ps) - 1)
        if d[i] == 0.0
            return ps[i]
        elseif d[i] * d[i + 1] < 0
            return ps[i] - d[i] * (ps[i + 1] - ps[i]) / (d[i + 1] - d[i])
        end
    end
    return NaN
end

function print_summary(stats)
    """Print summary statistics, including the pairwise crossings of the mean Binder curves."""
    println("\n" * "="^80)
    println("MEMORY-TO-TRIVIAL EA BINDER ANALYSIS SUMMARY (lambda = 0.3: λ_x = 0.21, λ_zz = 0.49)")
    println("="^80)
    println("\nSystem Parameters:")
    println("  λ_x = 0.21 (X measurement strength), λ_zz = 0.49 (ZZ measurement strength)")
    println("  P_x = P_zz scanning from $(minimum(stats.P_x)) to $(maximum(stats.P_x))")

    println("\nAnalysis Results:")

    for L in sort(unique(stats.L))
        L_data = filter(row -> row.L == L, stats)

        min_B = minimum(L_data.B_mean)
        max_B = maximum(L_data.B_mean)
        P_min = L_data[argmin(L_data.B_mean), :P_x]
        P_max = L_data[argmax(L_data.B_mean), :P_x]

        println("\n  L = $L:")
        println("    Sample points: $(nrow(L_data))")
        println("    Mean # of jobs (samples) per point: $(round(mean(L_data.n_samples), digits=1))")
        println("    Binder range: [$(round(min_B, digits=4)), $(round(max_B, digits=4))]")
        println("    Min B at P = $(round(P_min, digits=3))")
        println("    Max B at P = $(round(P_max, digits=3))")
    end

    println("\nPairwise crossings of the mean Binder curves (linear interpolation):")
    Ls = sort(unique(stats.L))
    for (i, L1) in enumerate(Ls), L2 in Ls[(i + 1):end]
        d1 = filter(row -> row.L == L1, stats)
        d2 = filter(row -> row.L == L2, stats)
        common = sort(collect(intersect(Set(d1.P_x), Set(d2.P_x))))
        length(common) < 2 && continue
        y1 = Dict(zip(d1.P_x, d1.B_mean))
        y2 = Dict(zip(d2.P_x, d2.B_mean))
        ps = first_crossing(common, [y1[p] - y2[p] for p in common])
        if isfinite(ps)
            @printf("  L = %2d vs %2d: P* = %.4f\n", L1, L2, ps)
        else
            @printf("  L = %2d vs %2d: no crossing on the grid\n", L1, L2)
        end
    end

    println("\n" * "="^80)
    println("Generated plots:")
    println("  - $(PREFIX)_binder_vs_p.{pdf,png}")
    println("  - $(PREFIX)_L{L}_binder_vs_p.{pdf,png} (for each L)")
    println("  - $(PREFIX)_summary.csv")
    println("="^80 * "\n")
end

function main()
    # Find the most recent results directory (or take the one given)
    if !isempty(ARGS)
        results_dir = ARGS[1]
    else
        results_dirs = filter(d -> startswith(d, PREFIX * "_results_") && isdir(d), readdir())
        if isempty(results_dirs)
            println("ERROR: No results directory found matching '$(PREFIX)_results_*'")
            return
        end
        results_dir = last(sort(results_dirs))  # Get most recent by timestamp
    end
    println("Loading results from: $results_dir")

    df = load_results(results_dir)

    if nrow(df) == 0
        println("ERROR: No data loaded. Check that the results directory contains $(PREFIX)_L*_lx0.21_lzz0.49_*.json files")
        return
    end

    println("Loaded $(nrow(df)) data points")
    stats = compute_statistics(df)

    writedlm("$(PREFIX)_summary.csv",
             vcat(permutedims(names(stats)), Matrix(stats)), ',')

    # Generate plots
    println("\nGenerating plots...")
    plot_binder_vs_px(stats)
    plot_binder_by_L(stats)

    # Print summary
    print_summary(stats)
end

main()
