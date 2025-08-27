#!/usr/bin/env julia
"""
Script to collect and analyze results from cluster runs
"""

using JSON, Statistics, DelimitedFiles
using Pkg; Pkg.activate(".")

function collect_results(output_dir="output")
    # Find all JSON result files
    result_files = filter(x -> endswith(x, ".json"), readdir(output_dir, join=true))
    
    if isempty(result_files)
        println("No result files found in $output_dir")
        return nothing
    end
    
    println("Found $(length(result_files)) result files")
    
    # Read all results
    results = []
    for file in result_files
        try
            data = JSON.parsefile(file)
            push!(results, data)
        catch e
            println("Error reading $file: $e")
        end
    end
    
    return results
end

function analyze_results(results)
    if isempty(results)
        println("No results to analyze")
        return
    end
    
    # Group by L and lambda
    grouped = Dict()
    for r in results
        key = (r["L"], r["lambda"])
        if !haskey(grouped, key)
            grouped[key] = []
        end
        push!(grouped[key], r)
    end
    
    # Calculate statistics for each (L, lambda) pair
    println("\nResults Summary:")
    println("="^80)
    println("L\tλ\tBinder\t\tStd\t\tN_samples")
    println("="^80)
    
    summary_data = []
    for ((L, λ), group) in sort(collect(grouped))
        binders = [r["binder"] for r in group]
        mean_binder = mean(binders)
        std_binder = std(binders)
        n_samples = length(binders)
        
        println("$L\t$(round(λ, digits=3))\t$(round(mean_binder, digits=6))\t$(round(std_binder, digits=6))\t$n_samples")
        
        push!(summary_data, (L, λ, mean_binder, std_binder, n_samples))
    end
    
    return summary_data
end

function save_summary(summary_data, output_file="results_summary.csv")
    if isempty(summary_data)
        println("No data to save")
        return
    end
    
    # Convert to matrix for writing
    data_matrix = hcat([d[1] for d in summary_data],  # L
                      [d[2] for d in summary_data],  # lambda
                      [d[3] for d in summary_data],  # mean_binder
                      [d[4] for d in summary_data],  # std_binder
                      [d[5] for d in summary_data])  # n_samples
    
    header = ["L", "lambda", "binder_mean", "binder_std", "n_samples"]
    
    writedlm(output_file, vcat(reshape(header, 1, :), data_matrix), ',')
    println("\nSummary saved to $output_file")
end

function create_plots_script(summary_data)
    # Create a Julia script for plotting
    plot_script = """
using DelimitedFiles, Plots

# Read the summary data
data = readdlm("results_summary.csv", ',', header=true)[1]
Ls = unique(data[:, 1])
lambdas = unique(data[:, 2])

# Create plot
plt = plot(title="Binder Parameter vs λ", 
          xlabel="λ (λₓ=λ, λ_zz=1-λ)", 
          ylabel="Binder Parameter",
          grid=true, size=(800, 600))

colors = [:red, :green, :purple, :orange, :blue]

for (i, L) in enumerate(Ls)
    # Filter data for this L
    mask = data[:, 1] .== L
    λs = data[mask, 2]
    binders = data[mask, 3]
    stds = data[mask, 4]
    
    # Sort by lambda
    sorted_idx = sortperm(λs)
    λs_sorted = λs[sorted_idx]
    binders_sorted = binders[sorted_idx]
    stds_sorted = stds[sorted_idx]
    
    color = colors[mod(i-1, length(colors)) + 1]
    
    plot!(plt, λs_sorted, binders_sorted,
          marker=:circle, markersize=4, linewidth=2,
          label="L = \$L", color=color,
          yerror=stds_sorted)
end

savefig(plt, "binder_vs_lambda.png")
display(plt)
println("Plot saved as binder_vs_lambda.png")
"""
    
    open("create_plots.jl", "w") do io
        write(io, plot_script)
    end
    
    println("Plot script created: create_plots.jl")
    println("Run with: julia create_plots.jl")
end

function main()
    println("Collecting results from cluster runs...")
    results = collect_results()
    
    if results !== nothing
        summary_data = analyze_results(results)
        save_summary(summary_data)
        create_plots_script(summary_data)
        
        println("\nAnalysis complete!")
        println("- Summary saved to results_summary.csv")
        println("- Plot script created: create_plots.jl")
    end
end

if !isinteractive()
    main()
end
