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
    """Load all JSON files from the results directory."""
    all_data = []
    
    for L in [8, 10, 12, 14, 16]
        L_dir = joinpath(results_dir, "L$L")
        if !isdir(L_dir)
            continue
        end
        
        json_files = filter(f -> endswith(f, ".json"), readdir(L_dir, join=true))
        
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
    """Create plots of Binder parameter vs P_x for different system sizes."""
    
    param_sets = unique(select(stats, [:lambda_x, :lambda_zz]))
    
    for param in eachrow(param_sets)
        param_data = filter(row -> row.lambda_x == param.lambda_x && row.lambda_zz == param.lambda_zz, stats)
        
        p1 = plot(xlabel=L"P_x = P_{zz}" * " (dephasing probability)", ylabel="Binder Parameter", 
                  title=L"\lambda_x=" * "$(param.lambda_x)" * L", \lambda_{zz}=" * "$(param.lambda_zz)",
                  legend=:topright, grid=true, size=(800, 600), dpi=300)
        
        colors = [:blue, :red, :green, :purple, :orange]
        markers = [:circle, :square, :diamond, :utriangle, :dtriangle]
        
        for (idx, L) in enumerate(sort(unique(param_data.L)))
            L_data = filter(row -> row.L == L, param_data)
            
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
        
        fname_base = "left_boundary_lx$(param.lambda_x)_lzz$(param.lambda_zz)"
        savefig(p1, "$(fname_base)_binder_vs_px.pdf")
        savefig(p1, "$(fname_base)_binder_vs_px.png")
    end
end

function main()
    results_dir = "left_boundary_results_20260208_2312"
    
    if !isdir(results_dir)
        return
    end
    
    df = load_results(results_dir)
    stats = compute_statistics(df)
    plot_binder_vs_px(stats)
end

main()
