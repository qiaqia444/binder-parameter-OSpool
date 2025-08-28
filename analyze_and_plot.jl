#!/usr/bin/env julia
"""
Analysis and plotting script for Binder parameter simulation results
Creates plots similar to the original binder_parameter.ipynb notebook
"""

using Pkg; Pkg.activate(".")
using JSON, Statistics, DelimitedFiles
using Plots

function load_simulation_results(output_dir="binder-simulation-results/output")
    """Load all JSON result files from the simulation"""
    
    if !isdir(output_dir)
        println("Error: Output directory $output_dir not found")
        return nothing
    end
    
    result_files = filter(x -> endswith(x, ".json"), readdir(output_dir, join=true))
    
    if isempty(result_files)
        println("No result files found in $output_dir")
        return nothing
    end
    
    # Read all results
    results = []
    failed_files = []
    
    for file in result_files
        try
            data = JSON.parsefile(file)
            push!(results, data)
        catch e
            push!(failed_files, file)
        end
    end
    
    println("Loaded $(length(results)) simulation results")
    return results
end

function process_results(results)
    """Process raw results into organized data for plotting"""
    
    # Group by L and lambda
    grouped = Dict()
    
    for r in results
        L = r["L"]
        λ = r["lambda"]
        key = (L, λ)
        
        if !haskey(grouped, key)
            grouped[key] = []
        end
        
        push!(grouped[key], r)
    end
    
    # Calculate statistics for each (L, lambda) pair
    processed_data = Dict()
    
    for ((L, λ), group) in grouped
        binders = [r["binder"] for r in group]
        
        if !isempty(binders)
            processed_data[(L, λ)] = (
                L = L,
                lambda = λ,
                binder_mean = mean(binders),
                binder_std = std(binders),
                binder_values = binders,
                n_samples = length(binders)
            )
        end
    end
    
    return processed_data
end

function create_binder_plots(processed_data)
    """Create plots similar to the original notebook"""
    
    # Get unique L values and sort them
    L_values = sort(unique([data.L for data in values(processed_data)]))
    
    # Create the main plot
    plt_combined = plot(title="Binder Parameter vs λ", 
                       xlabel="λ (λₓ=λ, λ_zz=1-λ)", 
                       ylabel="Binder Parameter",
                       grid=true, 
                       size=(800, 600),
                       legend=:topright)
    
    # Colors for different L values
    colors = [:red, :green, :purple, :orange, :blue]
    
    # Plot each L value
    for (i, L) in enumerate(L_values)
        # Extract data for this L
        L_data = [data for data in values(processed_data) if data.L == L]
        
        if isempty(L_data)
            continue
        end
        
        # Sort by lambda
        sort!(L_data, by = x -> x.lambda)
        
        λs = [data.lambda for data in L_data]
        binders = [data.binder_mean for data in L_data]
        errors = [data.binder_std for data in L_data]
        
        # Filter out any NaN values
        valid_idx = findall(x -> !isnan(x), binders)
        
        if !isempty(valid_idx)
            color = colors[mod(i-1, length(colors)) + 1]
            
            # Plot with error bars
            plot!(plt_combined, λs[valid_idx], binders[valid_idx],
                  marker=:circle, markersize=4, linewidth=2,
                  label="L = $L", color=color,
                  yerror=errors[valid_idx])
        end
    end
    
    # Save the plot
    savefig(plt_combined, "binder_vs_lambda_combined.png")
    
    # Create individual plots for each L
    for L in L_values
        L_data = [data for data in values(processed_data) if data.L == L]
        
        if isempty(L_data)
            continue
        end
        
        sort!(L_data, by = x -> x.lambda)
        
        λs = [data.lambda for data in L_data]
        binders = [data.binder_mean for data in L_data]
        errors = [data.binder_std for data in L_data]
        
        plt_L = plot(title="Binder Parameter vs λ (L = $L)", 
                     xlabel="λ (λₓ=λ, λ_zz=1-λ)", 
                     ylabel="Binder Parameter",
                     grid=true, 
                     size=(600, 450))
        
        plot!(plt_L, λs, binders,
              marker=:circle, markersize=5, linewidth=2,
              color=:blue, label="L = $L",
              yerror=errors)
        
        savefig(plt_L, "binder_vs_lambda_L$(L).png")
    end
    
    return plt_combined
end

function print_results_summary(processed_data)
    """Print detailed results summary like the original notebook"""
    
    L_values = sort(unique([data.L for data in values(processed_data)]))
    
    println("\nBINDER PARAMETER RESULTS SUMMARY")
    
    for L in L_values
        L_data = [data for data in values(processed_data) if data.L == L]
        sort!(L_data, by = x -> x.lambda)
        
        println("\nL = $L:")
        
        for data in L_data
            λ = data.lambda
            B_mean = data.binder_mean
            B_std = data.binder_std
            n_samples = data.n_samples
            
            println("  λ = $(round(λ, digits=3)): B = $(round(B_mean, digits=6)) ± $(round(B_std, digits=6)) (n=$n_samples)")
        end
    end
    
    println("\nTRANSITION REGION ANALYSIS (λ ∈ [0.4, 0.6]):")
    
    for L in L_values
        transition_data = [data for data in values(processed_data) 
                          if data.L == L && 0.4 <= data.lambda <= 0.6]
        
        if !isempty(transition_data)
            sort!(transition_data, by = x -> x.lambda)
            println("\nL = $L (transition region):")
            
            for data in transition_data
                λ = data.lambda
                B_mean = data.binder_mean
                println("  λ = $(round(λ, digits=3)): B = $(round(B_mean, digits=6))")
            end
        end
    end
end

function save_results_csv(processed_data, filename="binder_results_summary.csv")
    """Save results to CSV file for further analysis"""
    
    # Convert to matrix format
    data_rows = []
    
    for data in values(processed_data)
        push!(data_rows, [
            data.L,
            data.lambda,
            data.binder_mean,
            data.binder_std,
            data.n_samples
        ])
    end
    
    # Sort by L, then by lambda
    sort!(data_rows, by = x -> (x[1], x[2]))
    
    # Create header
    header = ["L", "lambda", "binder_mean", "binder_std", "n_samples"]
    
    # Write to CSV
    open(filename, "w") do io
        # Write header
        println(io, join(header, ","))
        
        # Write data
        for row in data_rows
            println(io, join(row, ","))
        end
    end
end

function main()
    println("Binder Parameter Results Analysis")
    
    # Load results
    results = load_simulation_results()
    
    if results === nothing
        println("No results to analyze. Make sure to download results first.")
        return
    end
    
    # Process results
    processed_data = process_results(results)
    
    if isempty(processed_data)
        println("No valid data found in results")
        return
    end
    
    println("Processed $(length(processed_data)) (L,λ) parameter points")
    
    # Create plots
    plt = create_binder_plots(processed_data)
    display(plt)
    
    # Print summary
    print_results_summary(processed_data)
    
    # Save CSV
    save_results_csv(processed_data)
    
    println("Analysis complete.")
end

if !isinteractive()
    main()
end
