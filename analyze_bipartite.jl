#!/usr/bin/env julia

"""
Bipartite Entropy Analysis Script
Analyzes results from weak measurement simulations and creates plots
"""

using JSON
using Plots
using Statistics
using DataFrames
using Printf
using Glob
using Dates

# Set plot backend for publication quality
gr()  # Use GR backend for better static plots that match matplotlib style

function load_bipartite_results(results_dir::String)
    """
    Load all bipartite entropy results from JSON files
    """
    println("Loading bipartite entropy results from: $results_dir")
    
    # Find all result JSON files (exclude summary and error files)
    json_files = glob("bipartite_L*.json", results_dir)
    json_files = filter(f -> !contains(f, "FAILED") && !contains(f, "summary"), json_files)
    
    println("Found $(length(json_files)) result files")
    
    results = []
    
    for file in json_files
        try
            data = JSON.parsefile(file)
            
            # Skip failed runs
            if !get(data, "success", false)
                println("  Skipping failed run: $(basename(file))")
                continue
            end
            
            push!(results, data)
            
        catch e
            println("  Warning: Could not load $file: $e")
        end
    end
    
    println("Successfully loaded $(length(results)) results")
    return results
end

function extract_central_entropies(results)
    """
    Extract central bipartite entropies organized by L and lambda
    """
    data_by_L = Dict()
    
    for result in results
        L = result["L"]
        lambda = result["lambda"]
        central_entropy = result["central_entropy"]
        
        # Skip NaN results
        if isnan(central_entropy)
            continue
        end
        
        if !haskey(data_by_L, L)
            data_by_L[L] = Dict()
        end
        
        data_by_L[L][lambda] = central_entropy
    end
    
    return data_by_L
end

function plot_entropy_vs_lambda(data_by_L, T_max=400; save_plot=true, plot_title="")
    """
    Create plot of bipartite entropy vs lambda parameter
    Exact styling to match the reference plot
    """
    
    # Sort system sizes
    L_values = sort(collect(keys(data_by_L)))
    
    # Create plot with exact styling from reference
    p = plot(
        title = "S vs λ, T_max=$T_max, loop=True",
        xlabel = "λ (X weak measurement strength)",
        ylabel = "Final bipartite entropy S(mid cut)",
        legend = :topright,
        grid = true,
        gridwidth = 0.8,
        gridcolor = :lightgray,
        gridalpha = 0.8,
        size = (800, 600),
        dpi = 300,
        background_color = :white,
        framestyle = :box,
        titlefontsize = 14,
        guidefontsize = 12,
        legendfontsize = 11,
        tickfontsize = 10,
        margin = 5Plots.mm
    )
    
    # Exact color scheme from reference plot
    # L=8: blue, L=10: orange, L=12: green
    plot_colors = Dict(
        8 => :blue,
        10 => :orange, 
        12 => :green,
        16 => :red,     # In case we have L=16
        20 => :purple,  # Additional sizes
        24 => :brown
    )
    
    for L in L_values
        lambda_vals = sort(collect(keys(data_by_L[L])))
        entropy_vals = [data_by_L[L][λ] for λ in lambda_vals]
        
        # Get color for this L value
        color = get(plot_colors, L, :black)
        
        # Plot with lines and circular markers (matching reference exactly)
        plot!(p, lambda_vals, entropy_vals,
              label = "L=$L",
              color = color,
              marker = :circle,
              markersize = 7,
              markerstrokewidth = 0.5,
              markerstrokecolor = color,
              markercolor = color,
              linewidth = 3,
              linestyle = :solid,
              alpha = 0.9)
    end
    
    # Set axis limits exactly like reference
    xlims!(p, (0.0, 1.0))
    ylims!(p, (0.0, 1.5))
    
    # Grid styling to match reference
    plot!(p, 
          xticks = 0.0:0.2:1.0,
          yticks = 0.0:0.2:1.4,
          minorticks = true)
    
    if save_plot
        timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
        filename = "bipartite_entropy_analysis_$timestamp.png"
        savefig(p, filename)
        println("Plot saved as: $filename")
    end
    
    return p
end

function analyze_entropy_scaling(data_by_L)
    """
    Analyze how entropy scales with system size L at different lambda values
    """
    println("\n=== Entropy Scaling Analysis ===")
    
    # Get all lambda values that appear in all system sizes
    all_lambdas = Set()
    for L in keys(data_by_L)
        union!(all_lambdas, keys(data_by_L[L]))
    end
    
    L_values = sort(collect(keys(data_by_L)))
    
    println("System sizes: $L_values")
    println("Lambda values analyzed: $(sort(collect(all_lambdas)))")
    
    # Create scaling analysis plot
    p_scaling = plot(
        title = "Entropy Scaling with System Size",
        xlabel = "System Size L",
        ylabel = "Central Bipartite Entropy S",
        legend = :topright,
        grid = true,
        size = (800, 600)
    )
    
    # Plot entropy vs L for different lambda values
    key_lambdas = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]  # Key lambda points
    colors = [:blue, :orange, :green, :red, :purple, :brown, :pink]
    
    for (i, λ) in enumerate(key_lambdas)
        if λ in all_lambdas
            L_vals = []
            S_vals = []
            
            for L in L_values
                if haskey(data_by_L[L], λ)
                    push!(L_vals, L)
                    push!(S_vals, data_by_L[L][λ])
                end
            end
            
            if length(L_vals) > 1
                plot!(p_scaling, L_vals, S_vals,
                      label = "λ = $λ",
                      color = colors[mod(i-1, length(colors)) + 1],
                      marker = :circle,
                      markersize = 4,
                      linewidth = 2)
            end
        end
    end
    
    # Save scaling plot
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
    savefig(p_scaling, "entropy_scaling_$timestamp.png")
    println("Scaling plot saved as: entropy_scaling_$timestamp.png")
    
    return p_scaling
end

function create_summary_report(results, data_by_L)
    """
    Create a summary report of the analysis
    """
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
    report_file = "bipartite_analysis_report_$timestamp.txt"
    
    open(report_file, "w") do f
        println(f, "="^60)
        println(f, "BIPARTITE ENTROPY ANALYSIS REPORT")
        println(f, "="^60)
        println(f, "Generated: $(Dates.now())")
        println(f, "")
        
        println(f, "SIMULATION PARAMETERS:")
        println(f, "-"^30)
        if !isempty(results)
            sample_result = results[1]
            measurement_type = get(sample_result, "measurement_type", "weak_bipartite")
            ntrials = get(sample_result, "ntrials", "unknown")
            T_max = get(sample_result, "T_max", "unknown")
            maxdim = get(sample_result, "maxdim", "unknown")
            println(f, "Measurement type: $measurement_type")
            println(f, "Trials per job: $ntrials")
            println(f, "T_max: $T_max")
            println(f, "Bond dimension: $maxdim")
        end
        println(f, "")
        
        println(f, "RESULTS SUMMARY:")
        println(f, "-"^30)
        println(f, "Total successful runs: $(length(results))")
        
        L_values = sort(collect(keys(data_by_L)))
        println(f, "System sizes: $L_values")
        
        for L in L_values
            lambdas = sort(collect(keys(data_by_L[L])))
            entropies = [data_by_L[L][λ] for λ in lambdas]
            
            println(f, "")
            println(f, "L = $L:")
            println(f, "  Lambda range: $(minimum(lambdas)) to $(maximum(lambdas))")
            println(f, "  Entropy range: $(round(minimum(entropies), digits=3)) to $(round(maximum(entropies), digits=3))")
            println(f, "  Peak entropy: $(round(maximum(entropies), digits=3)) at λ = $(lambdas[argmax(entropies)])")
        end
        
        println(f, "")
        println(f, "FILES GENERATED:")
        println(f, "-"^30)
        println(f, "- bipartite_entropy_analysis_$timestamp.png")
        println(f, "- entropy_scaling_$timestamp.png") 
        println(f, "- $report_file")
        println(f, "")
        println(f, "="^60)
    end
    
    println("Analysis report saved as: $report_file")
    return report_file
end

function main()
    if length(ARGS) < 1
        println("Usage: julia analyze_bipartite.jl <results_directory>")
        println("Example: julia analyze_bipartite.jl bipartite_results_20250930_1200")
        return
    end
    
    results_dir = ARGS[1]
    
    if !isdir(results_dir)
        println("Error: Directory '$results_dir' not found!")
        return
    end
    
    println("="^60)
    println("BIPARTITE ENTROPY ANALYSIS")
    println("="^60)
    
    # Load results
    results = load_bipartite_results(results_dir)
    
    if isempty(results)
        println("No valid results found!")
        return
    end
    
    # Extract and organize data
    data_by_L = extract_central_entropies(results)
    
    println("\nData organization:")
    for L in sort(collect(keys(data_by_L)))
        lambdas = sort(collect(keys(data_by_L[L])))
        println("  L=$L: $(length(lambdas)) lambda points ($(minimum(lambdas)) to $(maximum(lambdas)))")
    end
    
    # Create main entropy vs lambda plot
    println("\nCreating entropy vs lambda plot...")
    p_main = plot_entropy_vs_lambda(data_by_L)
    display(p_main)
    
    # Create scaling analysis
    println("\nCreating scaling analysis...")
    p_scaling = analyze_entropy_scaling(data_by_L)
    display(p_scaling)
    
    # Generate summary report
    println("\nGenerating summary report...")
    create_summary_report(results, data_by_L)
    
    println("\n" * "="^60)
    println("Analysis complete!")
    println("Check the generated PNG files and report for results.")
    println("="^60)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end