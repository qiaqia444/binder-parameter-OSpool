#!/usr/bin/env julia

"""
Plot Binder parameter for P=0.5 dephasing data (CORRECTED VERSION)
All system sizes: L=8,12,16,20,24,28,32,36
Corrected quantum channel implementation: ρ → (1-p)ρ + p·X·ρ·X
"""

using JSON
using Plots
using Statistics
using CSV
using DataFrames

# Load and process results from multiple directories
function load_results(directories)
    all_results = Dict{Int, Vector{Dict{String, Any}}}()
    
    for dir in directories
        if !isdir(dir)
            continue
        end
        
        files = filter(f -> endswith(f, ".json"), readdir(dir, join=true))
        if !isempty(files)
            println("Found $(length(files)) files in $dir")
        end
        
        for file in files
            try
                data = JSON.parsefile(file)
                L = data["L"]
                if !haskey(all_results, L)
                    all_results[L] = []
                end
                push!(all_results[L], data)
            catch e
                println("Error reading $file: $e")
            end
        end
    end
    
    return all_results
end

# Calculate average Binder parameter for each (L, lambda)
function calculate_binder_averages(results_dict)
    binder_data = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
    
    for L in sort(collect(keys(results_dict)))
        results = results_dict[L]
        # Group by lambda
        lambda_groups = Dict{Float64, Vector{Float64}}()
        
        for data in results
            # Use lambda_x (for P=0.5, lambda_x = lambda_zz = lambda)
            lambda = data["lambda_x"]
            binder = data["binder_parameter"]
            
            if !haskey(lambda_groups, lambda)
                lambda_groups[lambda] = []
            end
            push!(lambda_groups[lambda], binder)
        end
        
        # Calculate averages and standard errors
        lambdas = Float64[]
        binders = Float64[]
        errors = Float64[]
        
        for (lambda, binder_list) in sort(collect(lambda_groups))
            push!(lambdas, lambda)
            push!(binders, mean(binder_list))
            push!(errors, std(binder_list) / sqrt(length(binder_list)))
        end
        
        binder_data[L] = (lambdas, binders, errors)
        println("L=$L: $(length(lambdas)) lambda points, samples per point: $(length(lambda_groups[lambdas[1]]))")
    end
    
    return binder_data
end

# Main analysis
println("="^60)
println("Dephasing P=0.5 Analysis (CORRECTED VERSION)")
println("Proper quantum channel: ρ → (1-p)ρ + p·X·ρ·X")
println("="^60)

# Search for result directories (corrected version)
directories = [
    "dephasing_p05_results_20251118_0044",
    "large_dephasing_p05_results_20251120_2236",
    "large_dephasing_p05_results_20251121_1142"
]

# Also check jobs/results_* directories and any timestamped directories
jobs_dirs = filter(d -> startswith(basename(d), "results_"), 
                   filter(isdir, readdir("jobs", join=true)))
append!(directories, jobs_dirs)

# Also check for any new timestamped directories in current folder
for entry in readdir(".", join=false)
    if isdir(entry) && (startswith(entry, "large_dephasing_p05_results_") || 
                        startswith(entry, "dephasing_p05_results_"))
        if !(entry in directories)
            push!(directories, entry)
        end
    end
end

println("\nSearching directories:")
for dir in directories
    if isdir(dir)
        println("  - $dir")
    end
end
println()

# Load all results
all_results = load_results(directories)
println("\nLoaded results for L values: ", sort(collect(keys(all_results))))

# Calculate Binder parameter averages
binder_data = calculate_binder_averages(all_results)

# Create plot
println("\nGenerating plot...")
p = plot(
    xlabel="λ",
    ylabel="EA Binder Parameter",
    title="Dephasing P=0.5",
    legend=:topright,
    grid=true,
    size=(800, 600),
    dpi=300,
    framestyle=:box
)

# Colors and markers for different L
colors = [:blue, :red, :green, :purple, :orange, :brown, :pink, :gray]
markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5, :hexagon, :cross]

for (idx, L) in enumerate(sort(collect(keys(binder_data))))
    lambdas, binders, errors = binder_data[L]
    color = colors[mod1(idx, length(colors))]
    marker = markers[mod1(idx, length(markers))]
    
    plot!(p, lambdas, binders, 
          yerr=errors,
          marker=marker,
          markersize=5,
          color=color,
          linewidth=2,
          label="L=$L")
end

# Save plot
output_file = "plot_dephasing_p05_corrected.png"
savefig(p, output_file)
println("\nPlot saved to: $output_file")

# Save summary statistics to CSV
println("\nSaving summary statistics...")
df_rows = []
for L in sort(collect(keys(binder_data)))
    lambdas, binders, errors = binder_data[L]
    for i in 1:length(lambdas)
        push!(df_rows, (L=L, lambda=lambdas[i], binder=binders[i], error=errors[i]))
    end
end

df = DataFrame(df_rows)
csv_file = "dephasing_p05_corrected_summary.csv"
CSV.write(csv_file, df)
println("Summary saved to: $csv_file")

println("\n" * "="^60)
println("Analysis complete! (CORRECTED VERSION)")
println("="^60)
