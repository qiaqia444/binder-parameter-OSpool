#!/usr/bin/env julia

# Enhanced Binder Parameter Analysis Script for Forced Measurements
# Generates publication-quality plots with improved styling

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


function plot_binder_analysis_enhanced(lambdas::Vector{Float64}, 
                                      all_binder_data::Dict{Int, Vector{Float64}})

    p_binder = plot(title="Binder Parameter vs λ (Forced Measurements)",
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
    
    plot!(legendfontsize=10, guidefontsize=12, titlefontsize=11)

    return p_binder
end

# Convert grouped_stats to the format expected by your plotting function
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

# Generate plot using your style
p1 = plot_binder_analysis_enhanced(lambdas_vector, binder_data_dict)
savefig(p1, "notebook_style_binder_plot.png")
println("Saved: notebook_style_binder_plot.png")