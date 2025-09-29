#!/usr/bin/env julia

# Standard Trajectory Binder Parameter Analysis Script
# Analyzes standard quantum trajectory results with L*_lam* naming

using JSON
using Statistics
using Plots
using DataFrames
using CSV

# Set plot defaults for better appearance
default(fontfamily="Times", 
        linewidth=2,
        gridwidth=1,
        framestyle=:box,
        size=(800, 600),
        dpi=300)

println("=== Standard Trajectory Binder Parameter Analysis ===")

# Check if standard results directory exists
if !isdir("my_standard_results")
    println("ERROR: No my_standard_results directory found!")
    println("Make sure you've extracted the standard results archive first:")
    println("  tar -xzf my_standard_*.tar.gz")
    exit(1)
end

# Load standard trajectory results
function load_standard_results(system_size::Int)
    results_dir = "my_standard_results/L$system_size"
    if !isdir(results_dir)
        println("WARNING: No results found for L=$system_size")
        return DataFrame()
    end
    
    files = filter(f -> endswith(f, ".json"), readdir(results_dir))
    println("Loading $(length(files)) standard trajectory files for L=$system_size...")
    
    data = []
    for file in files
        try
            json_data = JSON.parsefile(joinpath(results_dir, file))
            if haskey(json_data, "binder")  # Check for binder field (standard trajectory format)
                push!(data, json_data)
            else
                println("  Skipping file without binder data: $file")
            end
        catch e
            println("  Error loading $file: $e")
        end
    end
    
    if isempty(data)
        return DataFrame()
    end
    
    # Convert to DataFrame - handle the actual JSON structure
    df = DataFrame()
    df.L = [get(d, "L", system_size) for d in data]
    df.lambda_x = [get(d, "lambda_x", NaN) for d in data] 
    df.lambda_zz = [get(d, "lambda_zz", NaN) for d in data]
    df.lambda = [get(d, "lambda", NaN) for d in data]
    # Use "binder" field instead of "binder_parameter"
    df.binder_parameter = [get(d, "binder", NaN) for d in data]
    df.sample = [get(d, "sample", 1) for d in data]
    df.seed = [get(d, "seed", 0) for d in data]
    
    return df
end

# Load data for all system sizes
println("Loading standard trajectory simulation results...")
df_L8 = load_standard_results(8)
df_L12 = load_standard_results(12)
df_L16 = load_standard_results(16)

# Combine all data
all_data = vcat(df_L8, df_L12, df_L16)

if nrow(all_data) == 0
    println("ERROR: No valid results found!")
    exit(1)
end

println("Loaded $(nrow(all_data)) successful standard trajectory simulations")
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

# Create enhanced plots using notebook style
function plot_binder_analysis_standard(lambdas::Vector{Float64}, 
                                      all_binder_data::Dict{Int, Vector{Float64}})

    p_binder = plot(title="Standard Trajectory Binder Parameter vs λ\\nLook for crossing points indicating critical transitions",
                   xlabel="λ", ylabel="Binder Parameter", 
                   grid=true, gridwidth=1, gridcolor=:lightgray, gridalpha=0.5,
                   size=(800, 600), margin=5Plots.mm)
    
    colors = [:blue, :red, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :triangle, :star5]
    
    all_x_data = Float64[]
    all_y_data = Float64[]
    
    for (idx, L) in enumerate(sort(collect(keys(all_binder_data))))
        color = colors[idx % length(colors) + 1]
        marker = markers[idx % length(markers) + 1]
        
        binder_data = all_binder_data[L]
        valid_idx = findall(!isnan, binder_data)
        if !isempty(valid_idx)
            x_vals = lambdas[valid_idx]
            y_vals = binder_data[valid_idx]
            
            plot!(p_binder, x_vals, y_vals,
                  marker=marker, markersize=6, linewidth=3,
                  label="L = $L", color=color, alpha=0.8)
            
            append!(all_x_data, x_vals)
            append!(all_y_data, y_vals)
        end
    end
    
    # Set dynamic axis limits with minimal padding for close data
    if !isempty(all_x_data) && !isempty(all_y_data)
        x_min, x_max = minimum(all_x_data), maximum(all_x_data)
        y_min, y_max = minimum(all_y_data), maximum(all_y_data)
        
        # For closely spaced data, use minimal padding
        x_range = x_max - x_min
        y_range = y_max - y_min
        
        x_padding = max(0.02 * x_range, 0.005)  # 2% padding, min 0.005
        y_padding = max(0.02 * y_range, 0.002)  # 2% padding, min 0.002
        
        if y_range < 0.01
            y_padding = max(y_padding, 0.001)
        end
        
        plot!(xlims=(x_min - x_padding, x_max + x_padding),
              ylims=(y_min - y_padding, y_max + y_padding))
    end
    
    # Add a vertical line at λ = 0.5 as a reference (if it's within the data range)
    if !isempty(all_x_data)
        x_min_data, x_max_data = minimum(all_x_data), maximum(all_x_data)
        if 0.5 >= x_min_data && 0.5 <= x_max_data
            y_min_plot, y_max_plot = ylims(p_binder)
            plot!(p_binder, [0.5, 0.5], [y_min_plot, y_max_plot], 
                  linestyle=:dot, color=:gray, alpha=0.7,
                  label="λ = 0.5 (balanced)", linewidth=1)
        end
    end
    
    plot!(legendfontsize=10, guidefontsize=12, titlefontsize=11)

    return p_binder
end

# Convert grouped_stats to the format expected by plotting function
binder_data_dict = Dict{Int, Vector{Float64}}()
lambdas_vector = sort(unique(grouped_stats.lambda))

for L in [8, 12, 16]
    L_data = filter(row -> row.L == L, grouped_stats)
    sort!(L_data, :lambda)
    
    # Create vector aligned with lambdas_vector
    binder_values = Float64[]
    for λ in lambdas_vector
        matching_row = filter(row -> row.lambda ≈ λ, L_data)
        if nrow(matching_row) > 0
            push!(binder_values, matching_row.B_mean[1])
        else
            push!(binder_values, NaN)
        end
    end
    binder_data_dict[L] = binder_values
end

# Generate plots
println("\\nGenerating standard trajectory plots...")
p1 = plot_binder_analysis_standard(lambdas_vector, binder_data_dict)
savefig(p1, "standard_binder_plot.png")
println("Saved: standard_binder_plot.png")

# Save processed data
CSV.write("standard_binder_results.csv", grouped_stats)
println("Saved: standard_binder_results.csv")

# Enhanced summary with critical point analysis
println("\\n=== Standard Trajectory Analysis Summary ===")
println("Total lambda values tested: $(length(unique(all_data.lambda)))")
println("Lambda range: $(minimum(all_data.lambda)) to $(maximum(all_data.lambda))")
println("Average samples per (L,λ): $(round(mean(grouped_stats.n_samples), digits=1))")

# Critical point analysis
println("\\n=== Critical Point Analysis ===")
critical_candidates = filter(row -> 0.55 <= row.B_mean <= 0.65, grouped_stats)
if nrow(critical_candidates) > 0
    println("Binder parameter ≈ 0.6 crossings (potential critical points):")
    for row in eachrow(critical_candidates)
        error = row.B_std / sqrt(row.n_samples)
        println("  L=$(row.L), λ=$(row.lambda): B = $(round(row.B_mean, digits=3)) ± $(round(error, digits=3))")
    end
    
    # Estimate critical lambda for each L
    println("\\nEstimated critical λc for each system size:")
    for L in [8, 12, 16]
        L_data = filter(row -> row.L == L, critical_candidates)
        if nrow(L_data) > 0
            # Find lambda closest to B = 0.6
            closest_idx = argmin(abs.(L_data.B_mean .- 0.6))
            λc_est = L_data.lambda[closest_idx]
            B_est = L_data.B_mean[closest_idx]
            println("  L=$L: λc ≈ $λc_est (B = $(round(B_est, digits=3)))")
        end
    end
end

println("\\n=== Generated Files ===")
println("  - standard_binder_plot.png (Main physics plot)")
println("  - standard_binder_results.csv (Statistical data)")

println("\\n=== Physics Conclusion ===")
println("Standard trajectory Binder parameter analysis complete!")
println("Ready for comparison with forced measurement results!")