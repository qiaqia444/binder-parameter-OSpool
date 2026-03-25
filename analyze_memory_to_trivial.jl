#!/usr/bin/env julia

"""
Analyze Memory-to-Trivial Transition Scan Results - Plot Generation Only

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
    """Load all JSON files from the results directory - only lx0.10_lzz0.70_P*.json files."""
    all_data = []
    
    for L in [8, 10, 12, 14, 16]
        L_dir = joinpath(results_dir, "L$L")
        if !isdir(L_dir)
            continue
        end
        
        # Only load files matching the pattern: memory_to_trivial_L*_lx0.10_lzz0.70_P*.json
        json_files = filter(f -> contains(f, "lx0.10_lzz0.70_P") && endswith(f, ".json"), 
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
    """Create plot of Binder parameter vs P for λ_x=0.1, λ_zz=0.7."""
    
    # Should only have one parameter set now (lambda_x=0.1, lambda_zz=0.7)
    p1 = plot(xlabel=L"P_x = P_{zz}" * " (dephasing probability)", ylabel="Binder Parameter", 
              title=L"\lambda_x = 0.1, \lambda_{zz} = 0.7",
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
    
    savefig(p1, "memory_to_trivial_lx0.1_lzz0.7_binder_vs_p.pdf")
    savefig(p1, "memory_to_trivial_lx0.1_lzz0.7_binder_vs_p.png")
end

function plot_binder_by_L(stats)
    """Create separate plots for each L."""
    colors = [:steelblue, :darkred, :darkgreen, :darkorchid, :darkorange]
    
    for (idx, L) in enumerate(sort(unique(stats.L)))
        L_data = filter(row -> row.L == L, stats)
        
        p = plot(xlabel=L"P_x = P_{zz}" * " (dephasing probability)", 
                ylabel="Binder Parameter",
                title=L"\lambda_x = 0.1, \lambda_{zz} = 0.7, $L = " * "$L",
                legend=false, grid=true, size=(800, 600), dpi=300)
        
        plot!(p, L_data.P_x, L_data.B_mean,
              yerr=L_data.B_sem,
              color=colors[idx],
              marker=:circle,
              markersize=8,
              linewidth=2.5,
              markerstrokewidth=1.5)
        
        hline!(p, [2/3], linestyle=:dash, color=:black, label="Pure state (B=2/3)", 
              linewidth=2, alpha=0.7)
        hline!(p, [2.0/(3.0*L)], linestyle=:dot, color=colors[idx], 
              label="Maximally mixed (B=$(round(2.0/(3.0*L), digits=3)))", linewidth=2, alpha=0.7)
        
        savefig(p, "memory_to_trivial_L$(L)_binder_vs_p.pdf")
        savefig(p, "memory_to_trivial_L$(L)_binder_vs_p.png")
    end
end

function print_summary(stats)
    """Print summary statistics."""
    println("\n" * "="^80)
    println("MEMORY-TO-TRIVIAL TRANSITION ANALYSIS SUMMARY")
    println("="^80)
    println("\nSystem Parameters:")
    println("  λ_x = 0.1 (X measurement strength)")
    println("  λ_zz = 0.7 (ZZ measurement strength)")
    println("  P_x = P_zz scanning from 0 to 0.5")
    
    println("\nAnalysis Results:")
    
    for L in sort(unique(stats.L))
        L_data = filter(row -> row.L == L, stats)
        
        # Find critical-like features
        min_B = minimum(L_data.B_mean)
        max_B = maximum(L_data.B_mean)
        P_min = L_data[argmin(L_data.B_mean), :P_x]
        P_max = L_data[argmax(L_data.B_mean), :P_x]
        
        println("\n  L = $L:")
        println("    Sample points: $(nrow(L_data))")
        println("    Mean # of MC trajectories per point: $(Int(round(mean(L_data.n_samples))))")
        println("    Binder range: [$(round(min_B, digits=4)), $(round(max_B, digits=4))]")
        println("    Min B at P = $(round(P_min, digits=3))")
        println("    Max B at P = $(round(P_max, digits=3))")
        println("    Pure state (P=0): B ≈ $(round(L_data[findfirst(x -> x ≈ 0.0, L_data.P_x), :B_mean], digits=4))")
    end
    
    println("\n" * "="^80)
    println("Generated plots:")
    println("  - memory_to_trivial_lx0.1_lzz0.7_binder_vs_p.{pdf,png}")
    println("  - memory_to_trivial_L{L}_binder_vs_p.{pdf,png} (for each L)")
    println("="^80 * "\n")
end

function main()
    # Find the most recent results directory
    results_dirs = filter(d -> startswith(d, "memory_to_trivial_results_") && isdir(d), readdir())
    if isempty(results_dirs)
        println("ERROR: No results directory found matching 'memory_to_trivial_results_*'")
        return
    end
    
    results_dir = last(sort(results_dirs))  # Get most recent by timestamp
    println("Loading results from: $results_dir")
    
    df = load_results(results_dir)
    
    if nrow(df) == 0
        println("ERROR: No data loaded. Check that results directory contains lx0.10_lzz0.70_P*.json files")
        return
    end
    
    println("Loaded $(nrow(df)) data points")
    stats = compute_statistics(df)
    
    # Generate plots
    println("\nGenerating plots...")
    plot_binder_vs_px(stats)
    plot_binder_by_L(stats)
    
    # Print summary
    print_summary(stats)
end

main()
