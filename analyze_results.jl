#!/usr/bin/env julia

# Binder Parameter Analysis Script for Forced Measurements
# Run this on your Mac after transferring results from OSPool

using JSON
using Statistics
using Plots
using DataFrames
using CSV

println("=== Binder Parameter Analysis (Forced Measurements) ===")

# Check if results directory exists
if !isdir("results")
    println("ERROR: No results directory found!")
    println("Make sure you've extracted the transferred archive first:")
    println("  tar -xzf forced_binder_results_*.tar.gz")
    exit(1)
end

# Load all results
function load_results(system_size::Int)
    results_dir = "results/forced_measurements/L$system_size"
    if !isdir(results_dir)
        println("WARNING: No results found for L=$system_size")
        return DataFrame()
    end
    
    files = filter(f -> endswith(f, ".json"), readdir(results_dir))
    println("Loading $(length(files)) files for L=$system_size...")
    
    data = []
    for file in files
        try
            json_data = JSON.parsefile(joinpath(results_dir, file))
            if json_data["success"]
                push!(data, json_data)
            else
                println("  Skipping failed job: $file")
            end
        catch e
            println("  Error loading $file: $e")
        end
    end
    
    if isempty(data)
        return DataFrame()
    end
    
    # Convert to DataFrame
    df = DataFrame()
    df.L = [d["L"] for d in data]
    df.lambda_x = [d["lambda_x"] for d in data] 
    df.lambda_zz = [d["lambda_zz"] for d in data]
    df.lambda = [d["lambda"] for d in data]
    df.binder_parameter = [d["binder_parameter"] for d in data]
    df.sample = [d["sample"] for d in data]
    df.seed = [d["seed"] for d in data]
    
    return df
end

# Load data for all system sizes
println("Loading simulation results...")
df_L8 = load_results(8)
df_L12 = load_results(12)
df_L16 = load_results(16)

# Combine all data
all_data = vcat(df_L8, df_L12, df_L16)

if nrow(all_data) == 0
    println("ERROR: No valid results found!")
    exit(1)
end

println("Loaded $(nrow(all_data)) successful simulations")
println("System sizes: L=8 ($(nrow(df_L8))), L=12 ($(nrow(df_L12))), L=16 ($(nrow(df_L16)))")

# Calculate statistics for each (L, lambda) combination
grouped_stats = combine(groupby(all_data, [:L, :lambda]), 
                       :binder_parameter => mean => :B_mean,
                       :binder_parameter => std => :B_std,
                       :binder_parameter => length => :n_samples)

println("\\nStatistics by system size and lambda:")
for group in groupby(grouped_stats, :L)
    L = group.L[1]
    println("L=$L: $(nrow(group)) lambda values, $(sum(group.n_samples)) total samples")
end

# Create plots
println("\\nGenerating plots...")

# 1. Binder parameter vs lambda for each L
p1 = plot(title="Binder Parameter vs Lambda (Forced Measurements)", 
          xlabel="λ", ylabel="Binder Parameter B")

for L in [8, 12, 16]
    data_L = filter(row -> row.L == L, grouped_stats)
    if nrow(data_L) > 0
        scatter!(p1, data_L.lambda, data_L.B_mean, 
                yerror=data_L.B_std, label="L=$L", markersize=4)
        plot!(p1, data_L.lambda, data_L.B_mean, label="", alpha=0.5)
    end
end

savefig(p1, "binder_vs_lambda_forced.png")
println("Saved: binder_vs_lambda_forced.png")

# 2. Data coverage heatmap
p2 = plot(title="Data Coverage (Samples per L, λ)", 
          xlabel="λ", ylabel="System Size L")

# Create coverage matrix
lambdas = sort(unique(grouped_stats.lambda))
Ls = [8, 12, 16] 
coverage_matrix = zeros(length(Ls), length(lambdas))

for (i, L) in enumerate(Ls)
    for (j, lambda) in enumerate(lambdas)
        matches = filter(row -> row.L == L && row.lambda == lambda, grouped_stats)
        coverage_matrix[i, j] = nrow(matches) > 0 ? matches.n_samples[1] : 0
    end
end

heatmap!(p2, lambdas, Ls, coverage_matrix, color=:viridis)
savefig(p2, "data_coverage_forced.png") 
println("Saved: data_coverage_forced.png")

# 3. Save processed data
CSV.write("processed_binder_results_forced.csv", grouped_stats)
println("Saved: processed_binder_results_forced.csv")

# Summary statistics
println("\\n=== Analysis Summary ===")
println("Total lambda values tested: $(length(unique(all_data.lambda)))")
println("Lambda range: $(minimum(all_data.lambda)) to $(maximum(all_data.lambda))")
println("Average samples per (L,λ): $(round(mean(grouped_stats.n_samples), digits=1))")

# Look for transition region (where Binder parameter ≈ 0.6)
transition_candidates = filter(row -> 0.55 <= row.B_mean <= 0.65, grouped_stats)
if nrow(transition_candidates) > 0
    println("\\nPotential transition region (B ≈ 0.6):")
    for row in eachrow(transition_candidates)
        println("  L=$(row.L), λ=$(row.lambda): B = $(round(row.B_mean, digits=3)) ± $(round(row.B_std, digits=3))")
    end
end

println("\\nAnalysis complete! Generated files:")
println("  - binder_vs_lambda_forced.png")  
println("  - data_coverage_forced.png")
println("  - processed_binder_results_forced.csv")
println("\\nNext steps:")
println("  1. Compare with standard quantum trajectory results")
println("  2. Perform finite-size scaling analysis") 
println("  3. Identify critical lambda_c from crossing points")