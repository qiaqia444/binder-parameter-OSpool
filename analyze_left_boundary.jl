#!/usr/bin/env julia

"""
Analyze Left Boundary Scan Results - Plot Generation Only

Load all JSON files, compute statistics, and generate plots.
Properly groups by (L, lambda_x, lambda_zz, P_x) to avoid mixing different parameter sets.
"""

using JSON
using Statistics
using DataFrames
using Printf
using Plots
using LaTeXStrings

function load_results(results_dir)
    """Load all JSON files from the results directory - only lx0.00_lzz0.70_P*.json files."""
    all_data = []
    
    for L in [8, 16, 24, 32]
        L_dir = joinpath(results_dir, "L$L")
        if !isdir(L_dir)
            continue
        end
        
        # Only load files matching the pattern: left_boundary_L*_lx0.00_lzz0.70_P*.json
        json_files = filter(f -> contains(f, "lx0.00_lzz0.70_P") && endswith(f, ".json"), 
                           readdir(L_dir, join=true))
        
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
                # Skip failed files silently
            end
        end
    end
    
    return DataFrame(all_data)
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
    """Create plot of Binder parameter vs P for λ_x=0, λ_zz=0.7."""
    
    # Should only have one parameter set now (lambda_x=0, lambda_zz=0.7)
    p1 = plot(xlabel=L"P_x = P_{zz}" * " (dephasing probability)", ylabel="Binder Parameter", 
              title=L"\lambda_x = 0.0, \lambda_{zz} = 0.7",
              legend=:topright, grid=true, size=(800, 600), dpi=300)
    
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
        
        B_mixed = 2.0 / (3.0 * L)
        hline!(p1, [B_mixed], linestyle=:dot, color=colors[idx], 
               label=nothing, linewidth=1, alpha=0.3)
    end
    
    hline!(p1, [2/3], linestyle=:dash, color=:black, label="B = 2/3 (pure)", linewidth=2, alpha=0.7)
    
    savefig(p1, "1.pdf")
    savefig(p1, "1.png")
end

function main()
    # Find the most recent results directory
    results_dirs = filter(d -> startswith(d, "left_boundary_results_") && isdir(d), readdir())
    if isempty(results_dirs)
        println("ERROR: No results directory found matching 'left_boundary_results_*'")
        return
    end
    
    results_dir = last(sort(results_dirs))  # Get most recent by timestamp
    println("Loading results from: $results_dir")
    
    df = load_results(results_dir)
    
    if nrow(df) == 0
        println("ERROR: No data loaded. Check that results directory contains lx0.00_lzz0.70_P*.json files")
        return
    end
    
    println("Loaded $(nrow(df)) data points")
    stats = compute_statistics(df)
    plot_binder_vs_px(stats)
end

main()
