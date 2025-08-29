#!/usr/bin/env julia
"""
Analysis and plotting script for L=12 Binder parameter simulation results
"""

using JSON, Statistics, DelimitedFiles
using Plots

function load_L12_results(output_dir="output")
    """Load all L=12 JSON result files from the simulation"""
    
    if !isdir(output_dir)
        println("Error: Output directory $output_dir not found")
        return nothing
    end
    
    # Get only L12 files
    result_files = filter(x -> startswith(x, "L12_") && endswith(x, ".json"), readdir(output_dir))
    result_files = [joinpath(output_dir, f) for f in result_files]
    
    if isempty(result_files)
        println("No L12 result files found in $output_dir")
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
            println("Failed to read $file: $e")
            push!(failed_files, file)
        end
    end
    
    println("Loaded $(length(results)) L=12 simulation results")
    if !isempty(failed_files)
        println("Failed to load $(length(failed_files)) files")
    end
    
    return results
end

function process_L12_results(results)
    """Process L=12 results into organized data for plotting"""
    
    # Group by lambda
    grouped = Dict()
    
    for r in results
        λ = r["lambda"]
        
        if !haskey(grouped, λ)
            grouped[λ] = []
        end
        
        push!(grouped[λ], r)
    end
    
    # Calculate statistics for each lambda value
    processed_data = []
    
    for (λ, group) in grouped
        binders = [r["binder"] for r in group]
        
        if !isempty(binders)
            push!(processed_data, (
                lambda = λ,
                binder_mean = mean(binders),
                binder_std = std(binders),
                binder_values = binders,
                n_samples = length(binders)
            ))
        end
    end
    
    # Sort by lambda
    sort!(processed_data, by = x -> x.lambda)
    
    return processed_data
end

function create_L12_binder_plot(processed_data)
    """Create Binder parameter plot for L=12"""
    
    lambdas = [d.lambda for d in processed_data]
    binder_means = [d.binder_mean for d in processed_data]
    binder_stds = [d.binder_std for d in processed_data]
    
    # Create the plot
    plt = plot(title="Edwards-Anderson Binder Parameter (L=12)", 
               xlabel="λ (λₓ=λ, λ_zz=1-λ)", 
               ylabel="Binder Parameter",
               grid=true, 
               size=(800, 600),
               legend=:topright)
    
    # Plot with error bars
    plot!(plt, lambdas, binder_means, yerror=binder_stds,
          marker=:circle, markersize=6, linewidth=2,
          label="L=12", color=:blue)
    
    return plt
end

function save_results_csv(processed_data, filename="analysis_results/binder_results_L12.csv")
    """Save results to CSV file"""
    
    # Create analysis_results directory if it doesn't exist
    mkpath("analysis_results")
    
    # Create CSV data
    csv_data = []
    push!(csv_data, ["lambda", "binder_mean", "binder_std", "n_samples"])
    
    for d in processed_data
        push!(csv_data, [d.lambda, d.binder_mean, d.binder_std, d.n_samples])
    end
    
    # Write to file
    writedlm(filename, csv_data, ',')
    println("Results saved to $filename")
end

function main()
    println("=== Edwards-Anderson Binder Parameter Analysis (L=12) ===")
    println()
    
    # Load results
    println("Loading simulation results...")
    results = load_L12_results()
    
    if results === nothing
        println("No results found. Exiting.")
        return
    end
    
    # Process results
    println("Processing results...")
    processed_data = process_L12_results(results)
    
    println("Found $(length(processed_data)) lambda values:")
    for d in processed_data
        println("  λ=$(d.lambda): Binder=$(round(d.binder_mean, digits=4)) ± $(round(d.binder_std, digits=4)) (n=$(d.n_samples))")
    end
    println()
    
    # Create plot
    println("Creating plot...")
    plt = create_L12_binder_plot(processed_data)
    
    # Create analysis_results directory if it doesn't exist
    mkpath("analysis_results")
    
    # Save plot
    savefig(plt, "analysis_results/binder_vs_lambda_L12.png")
    println("Plot saved as analysis_results/binder_vs_lambda_L12.png")
    
    # Save CSV
    save_results_csv(processed_data)
    
    # Display plot (if in interactive mode)
    display(plt)
    
    println()
    println("=== Analysis Complete ===")
    println("Files created:")
    println("  - analysis_results/binder_vs_lambda_L12.png")
    println("  - analysis_results/binder_results_L12.csv")
end

# Run the analysis
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
