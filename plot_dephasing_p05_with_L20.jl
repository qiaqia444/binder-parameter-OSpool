#!/usr/bin/env julia

"""
Plot Binder parameter for P=0.5 dephasing data including L=20
"""

using JSON
using Statistics
using Plots

println("="^70)
println("Loading Dephasing P=0.5 Results (L=8,12,16,20)")
println("="^70)

# Search for result files in multiple locations
result_files = String[]

# Small systems (L=8,12,16)
small_dirs = [
    "dephasing_p05_results_20251118_0044",
    "jobs/results_dephasing_p05"
]

# Large systems (L=20)
large_dirs = [
    "large_dephasing_p05_results_20251120_2236",
    "jobs/results_large_dephasing_p05"
]

all_dirs = vcat(small_dirs, large_dirs)

for dir in all_dirs
    if isdir(dir)
        files = filter(f -> endswith(f, ".json"), readdir(dir, join=true))
        if !isempty(files)
            println("Found $(length(files)) files in $dir")
            append!(result_files, files)
        end
    end
end

println("\nTotal files found: $(length(result_files))")

if isempty(result_files)
    error("No result files found!")
end

# Parse results
results = Dict()

for file in result_files
    try
        data = JSON.parsefile(file)
        L = data["L"]
        lambda = haskey(data, "lambda") ? data["lambda"] : data["lambda_x"]
        P = data["P_x"]
        binder = data["binder_parameter"]
        
        # Only include P=0.5 data
        if P ≈ 0.5
            key = (L, lambda)
            if !haskey(results, key)
                results[key] = Float64[]
            end
            push!(results[key], binder)
        end
    catch e
        println("Warning: Could not parse $file: $e")
    end
end

println("Parsed $(length(results)) unique (L, λ) combinations")

# Get L values
L_values = sort(unique([k[1] for k in keys(results)]))
println("\nSystem sizes: $L_values")

# Aggregate data for each L
plot_data = Dict()

for L in L_values
    lambda_vals = Float64[]
    binder_means = Float64[]
    binder_stds = Float64[]
    
    for ((l, lam), binders) in sort(collect(results), by=x->x[1][2])
        if l == L
            push!(lambda_vals, lam)
            push!(binder_means, mean(binders))
            push!(binder_stds, std(binders) / sqrt(length(binders)))
        end
    end
    
    plot_data[L] = (lambda_vals, binder_means, binder_stds)
    println("  L=$L: $(length(lambda_vals)) λ points, $(length(results[(L, lambda_vals[1])])) samples/point")
end

# Create plot
println("\nCreating plot...")

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

colors = [:blue, :red, :green, :purple]
markers = [:circle, :square, :diamond, :utriangle]

for (idx, L) in enumerate(L_values)
    lambda_vals, binder_means, binder_stds = plot_data[L]
    
    plot!(p, lambda_vals, binder_means,
          yerr=binder_stds,
          marker=markers[idx],
          markersize=5,
          color=colors[idx],
          linewidth=2,
          label="L=$L")
end

savefig(p, "plot_dephasing_p05_with_L20.png")
println("\n✓ Plot saved to: plot_dephasing_p05_with_L20.png")

# Print summary
println("\n" * "="^70)
println("Summary Statistics")
println("="^70)

for L in L_values
    lambda_vals, binder_means, binder_stds = plot_data[L]
    println("L=$L:")
    println("  λ points: $(length(lambda_vals))")
    println("  λ range: [$(minimum(lambda_vals)), $(maximum(lambda_vals))]")
    println("  B range: [$(round(minimum(binder_means), digits=3)), $(round(maximum(binder_means), digits=3))]")
    println("  Mean error: $(round(mean(binder_stds), digits=4))")
end

println("="^70)
