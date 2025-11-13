#!/usr/bin/env julia

"""
Plot Binder parameter for P=0.5 dephasing data
"""

using JSON
using Statistics
using Plots

function load_p05_results()
    """Load P=0.5 results and compute statistics."""
    
    # Try both possible locations
    result_dirs = [
        "dephasing_p05_results_20251111_1727",
        "jobs/results_dephasing_p05"
    ]
    
    results = []
    for dir in result_dirs
        if isdir(dir)
            println("Loading from: $dir")
            json_files = filter(f -> endswith(f, ".json"), readdir(dir, join=true))
            println("Found $(length(json_files)) files")
            
            for file in json_files
                try
                    data = JSON.parsefile(file)
                    push!(results, data)
                catch e
                    @warn "Failed to parse $file"
                end
            end
            
            if !isempty(results)
                break
            end
        end
    end
    
    if isempty(results)
        error("No P=0.5 results found!")
    end
    
    # Group by L and lambda
    grouped = Dict()
    for r in results
        L = r["L"]
        lambda = haskey(r, "lambda") ? r["lambda"] : r["lambda_x"]
        B = r["binder_parameter"]
        
        key = (L, lambda)
        if !haskey(grouped, key)
            grouped[key] = Float64[]
        end
        push!(grouped[key], B)
    end
    
    # Compute statistics
    data_points = []
    for ((L, lambda), B_vals) in grouped
        push!(data_points, (
            L = L,
            lambda = lambda,
            B_mean = mean(B_vals),
            B_std = std(B_vals),
            n = length(B_vals)
        ))
    end
    
    return data_points
end

function plot_p05()
    """Create plot for P=0.5 data."""
    
    println("="^70)
    println("Plotting P=0.5 Dephasing Results")
    println("="^70)
    
    # Load data
    data = load_p05_results()
    println("\nLoaded $(length(data)) data points")
    
    # Create plot
    p = plot(
        xlabel = "λ",
        ylabel = "EA Binder Parameter",
        title = "Dephasing P=0.5",
        legend = :topright,
        size = (800, 600),
        dpi = 300,
        grid = true
    )
    
    # Get unique L values
    L_values = sort(unique([d.L for d in data]))
    colors = [:blue, :red, :green]
    
    # Plot each system size
    for (i, L) in enumerate(L_values)
        # Filter data for this L
        L_data = filter(d -> d.L == L, data)
        
        # Sort by lambda
        sort!(L_data, by = d -> d.lambda)
        
        # Extract arrays
        lambdas = [d.lambda for d in L_data]
        B_means = [d.B_mean for d in L_data]
        B_stds = [d.B_std for d in L_data]
        
        # Plot
        plot!(p, lambdas, B_means,
            ribbon = B_stds,
            fillalpha = 0.2,
            label = "L=$L",
            marker = :circle,
            markersize = 4,
            linewidth = 2,
            color = colors[i]
        )
    end
    
    # Save
    savefig(p, "plot_dephasing_p05.png")
    println("\n✓ Plot saved to: plot_dephasing_p05.png")
    
    # Print statistics
    println("\nStatistics:")
    for L in L_values
        L_data = filter(d -> d.L == L, data)
        B_vals = [d.B_mean for d in L_data]
        println("  L=$L: B ∈ [$(round(minimum(B_vals), digits=3)), $(round(maximum(B_vals), digits=3))]")
    end
    
    return p
end

# Run
if abspath(PROGRAM_FILE) == @__FILE__
    plot_p05()
end
