#!/usr/bin/env julia

using JSON, DataFrames, Statistics, Plots, CSV
using Printf

"""
Analysis of corrected MIPT results with fixed weak measurements.
This should show proper crossing behavior at λ = 0.5.
"""

function load_corrected_results(output_dir="output")
    results = []
    
    # Find all result files
    files = filter(f -> startswith(f, "L") && endswith(f, ".json"), readdir(output_dir))
    
    println("Found $(length(files)) result files")
    
    for file in files
        filepath = joinpath(output_dir, file)
        try
            data = JSON.parsefile(filepath)
            
            # Parse filename: L8_lam0.1_s1.json
            parts = split(replace(file, ".json" => ""), "_")
            L = parse(Int, replace(parts[1], "L" => ""))
            lambda_str = replace(parts[2], "lam" => "")
            lambda_val = parse(Float64, lambda_str)
            sample = parse(Int, replace(parts[3], "s" => ""))
            
            push!(results, (
                L = L,
                lambda = lambda_val,
                sample = sample,
                B = data["binder"],
                B_mean = get(data, "binder_mean_of_trials", data["binder"]),
                B_std = get(data, "binder_std_of_trials", 0.0),
                ntrials = get(data, "ntrials", 2000),
                filename = file
            ))
        catch e
            println("Warning: Could not parse $file: $e")
        end
    end
    
    return DataFrame(results)
end

function analyze_crossing_point(df)
    println("\n=== Crossing Point Analysis ===")
    
    # Group by lambda and compute means across samples
    grouped = combine(groupby(df, [:L, :lambda]), 
                     :B_mean => mean => :B_mean,
                     :B_std => mean => :B_std,
                     :sample => length => :n_samples)
    
    # Check crossing behavior around λ = 0.5
    critical_region = filter(row -> 0.45 <= row.lambda <= 0.55, grouped)
    
    println("\nBinder parameters in critical region (λ = 0.45 to 0.55):")
    for row in eachrow(sort(critical_region, [:lambda, :L]))
        println("λ = $(row.lambda), L = $(row.L): B = $(round(row.B_mean, digits=4)) ± $(round(row.B_std, digits=4))")
    end
    
    # Check for crossing at λ = 0.5
    lambda_05_data = filter(row -> row.lambda == 0.5, grouped)
    if nrow(lambda_05_data) == 3  # Should have L = 8, 12, 16
        println("\n*** Crossing Point Analysis at λ = 0.5 ***")
        sort!(lambda_05_data, :L)
        B_values = lambda_05_data.B_mean
        spread = maximum(B_values) - minimum(B_values)
        
        for row in eachrow(lambda_05_data)
            println("L = $(row.L): B = $(round(row.B_mean, digits=4)) ± $(round(row.B_std, digits=4))")
        end
        
        println("Spread at λ = 0.5: $(round(spread, digits=4))")
        
        if spread < 0.05
            println("✓ EXCELLENT crossing! Very small spread.")
        elseif spread < 0.1
            println("✓ Good crossing behavior.")
        elseif spread < 0.2
            println("~ Moderate crossing behavior.")
        else
            println("✗ No clear crossing - large spread.")
        end
    end
    
    return grouped
end

function plot_corrected_results(df)
    println("\nCreating plots...")
    
    # Group by lambda and compute means
    grouped = combine(groupby(df, [:L, :lambda]), 
                     :B_mean => mean => :B_mean,
                     :B_std => mean => :B_std)
    
    # Main crossing plot
    p1 = plot(xlabel="λ", ylabel="Binder Parameter B", 
              title="Corrected MIPT: Binder Parameter vs λ",
              legend=:topright, dpi=300, size=(800, 600))
    
    colors = [:red, :blue, :green]
    markers = [:circle, :square, :diamond]
    
    Ls = sort(unique(grouped.L))
    for (i, L) in enumerate(Ls)
        L_data = filter(row -> row.L == L, grouped)
        sort!(L_data, :lambda)
        
        plot!(p1, L_data.lambda, L_data.B_mean,
              label="L = $L",
              color=colors[i],
              marker=markers[i],
              markersize=4,
              linewidth=2,
              yerror=L_data.B_std)
    end
    
    # Add vertical line at λ = 0.5
    vline!(p1, [0.5], linestyle=:dash, color=:gray, alpha=0.7, 
           label="λ = 0.5", linewidth=2)
    
    savefig(p1, "analysis_results/corrected_crossing_point_analysis.png")
    println("Saved: analysis_results/corrected_crossing_point_analysis.png")
    
    # Zoom in on critical region
    p2 = plot(xlabel="λ", ylabel="Binder Parameter B", 
              title="Critical Region Detail (λ = 0.45-0.55)",
              legend=:topright, dpi=300, size=(800, 600))
    
    for (i, L) in enumerate(Ls)
        L_data = filter(row -> row.L == L && 0.45 <= row.lambda <= 0.55, grouped)
        sort!(L_data, :lambda)
        
        plot!(p2, L_data.lambda, L_data.B_mean,
              label="L = $L",
              color=colors[i],
              marker=markers[i],
              markersize=6,
              linewidth=2,
              yerror=L_data.B_std)
    end
    
    vline!(p2, [0.5], linestyle=:dash, color=:gray, alpha=0.7, 
           label="λ = 0.5", linewidth=2)
    
    savefig(p2, "analysis_results/corrected_critical_region_detail.png")
    println("Saved: analysis_results/corrected_critical_region_detail.png")
    
    return p1, p2
end

function main()
    println("=== Corrected MIPT Results Analysis ===")
    println("Loading results with fixed weak measurement operators...")
    
    # Load data
    df = load_corrected_results()
    println("Loaded $(nrow(df)) individual results")
    
    # Analyze crossing point
    grouped = analyze_crossing_point(df)
    
    # Create plots
    plot_corrected_results(df)
    
    # Save summary
    CSV.write("analysis_results/corrected_results_summary.csv", grouped)
    println("\nSaved summary to: analysis_results/corrected_results_summary.csv")
    
    println("\n=== Analysis Complete ===")
    println("Check the plots to see if the crossing point at λ = 0.5 is now visible!")
end

main()
