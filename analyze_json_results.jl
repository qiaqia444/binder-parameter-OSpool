#!/usr/bin/env julia

using JSON, DataFrames, Statistics, Plots, CSV
using Printf

"""
Analyze JSON format results from Edwards-Anderson Binder parameter simulations
Handles the JSON output format from your cluster runs
"""

function load_json_results(output_dir="output")
    println("Loading JSON results from $output_dir...")
    
    all_data = DataFrame(L=Int[], lambda=Float64[], lambda_x=Float64[], lambda_zz=Float64[], 
                        B=Float64[], B_mean=Float64[], B_std=Float64[], 
                        S2_bar=Float64[], S4_bar=Float64[], ntrials=Int[], 
                        seed=Int[], filename=String[])
    
    json_files = filter(f -> endswith(f, ".json"), readdir(output_dir))
    println("Found $(length(json_files)) JSON files")
    
    for file in json_files
        filepath = joinpath(output_dir, file)
        try
            # Read JSON file
            data = JSON.parsefile(filepath)
            
            # Parse filename to extract parameters
            # Expected format: std_L8_lam0.1_s1.json
            parts = split(replace(file, ".json" => ""), "_")
            if length(parts) >= 4 && parts[1] == "std"
                L = parse(Int, parts[2][2:end])  # Remove 'L' prefix from L8
                lambda_str = parts[3][4:end]     # Remove 'lam' prefix
                lambda = parse(Float64, lambda_str)
                seed_str = parts[4][2:end]       # Remove 's' prefix
                seed = parse(Int, seed_str)
                
                # Extract data from JSON  
                B = get(data, "binder", NaN)
                B_mean = get(data, "binder_mean_of_trials", NaN)
                B_std = get(data, "binder_std_of_trials", NaN)
                S2_bar = get(data, "S2_bar", NaN)
                S4_bar = get(data, "S4_bar", NaN)
                ntrials = get(data, "ntrials", 0)
                
                # Calculate lambda_x and lambda_zz
                lambda_x = lambda
                lambda_zz = 1.0 - lambda
                
                # Add to dataframe
                push!(all_data, (L, lambda, lambda_x, lambda_zz, B, B_mean, B_std, 
                                S2_bar, S4_bar, ntrials, seed, file))
            end
        catch e
            println("Warning: Could not read $file: $e")
        end
    end
    
    println("Loaded $(nrow(all_data)) data points")
    if nrow(all_data) > 0
        println("System sizes: $(sort(unique(all_data.L)))")
        println("Lambda range: $(round(minimum(all_data.lambda), digits=3)) - $(round(maximum(all_data.lambda), digits=3))")
        println("Number of trials per job: $(unique(all_data.ntrials))")
    end
    
    return all_data
end

function compute_statistics(data)
    """Compute mean and std for each (L, lambda) combination"""
    if nrow(data) == 0
        return DataFrame()
    end
    
    grouped = groupby(data, [:L, :lambda])
    
    stats = combine(grouped) do df
        DataFrame(
            lambda_x = df.lambda_x[1],
            lambda_zz = df.lambda_zz[1],
            B_mean = mean(df.B),
            B_std = std(df.B),
            B_sem = length(df.B) > 1 ? std(df.B) / sqrt(length(df.B)) : 0.0,
            n_seeds = length(df.B),
            S2_mean = mean(df.S2_bar),
            S4_mean = mean(df.S4_bar),
            ntrials_mean = mean(df.ntrials)
        )
    end
    
    return sort(stats, [:L, :lambda])
end

function plot_available_data(stats_data)
    """Create plots for whatever data is available"""
    
    if nrow(stats_data) == 0
        println("No data to plot!")
        return nothing
    end
    
    L_values = sort(unique(stats_data.L))
    println("Creating plots for L = $L_values")
    
    # Main Binder parameter plot
    p1 = plot(title="Edwards-Anderson Binder Parameter vs λ", 
              xlabel="λ (λₓ=λ, λ_zz=1-λ)", ylabel="Binder Parameter B",
              grid=true, size=(800, 600), dpi=300)
    
    colors = [:red, :blue, :green, :purple, :orange]
    
    for (i, L) in enumerate(L_values)
        L_data = filter(row -> row.L == L, stats_data)
        color = colors[(i-1) % length(colors) + 1]
        
        # Plot with error bars if we have multiple seeds
        if maximum(L_data.n_seeds) > 1
            plot!(p1, L_data.lambda, L_data.B_mean,
                  yerror=L_data.B_sem,
                  marker=:circle, markersize=4, linewidth=2,
                  label="L = $L", color=color)
        else
            plot!(p1, L_data.lambda, L_data.B_mean,
                  marker=:circle, markersize=4, linewidth=2,
                  label="L = $L", color=color)
        end
    end
    
    # If we have data around the critical region, zoom in
    lambda_min = minimum(stats_data.lambda)
    lambda_max = maximum(stats_data.lambda)
    
    if lambda_min <= 0.4 && lambda_max >= 0.6
        p2 = plot(title="Critical Region Detail", 
                  xlabel="λ", ylabel="Binder Parameter B",
                  grid=true, size=(600, 400), dpi=300,
                  xlims=(max(0.3, lambda_min), min(0.7, lambda_max)))
        
        for (i, L) in enumerate(L_values)
            L_data = filter(row -> row.L == L, stats_data)
            color = colors[(i-1) % length(colors) + 1]
            
            if maximum(L_data.n_seeds) > 1
                plot!(p2, L_data.lambda, L_data.B_mean,
                      yerror=L_data.B_sem,
                      marker=:circle, markersize=6, linewidth=3,
                      label="L = $L", color=color)
            else
                plot!(p2, L_data.lambda, L_data.B_mean,
                      marker=:circle, markersize=6, linewidth=3,
                      label="L = $L", color=color)
            end
        end
        
        return p1, p2
    end
    
    return p1, nothing
end

function print_data_summary(data, stats_data)
    """Print summary of what data we have"""
    println("\n" * "="^60)
    println("DATA SUMMARY")
    println("="^60)
    
    if nrow(data) == 0
        println("No data found!")
        return
    end
    
    println("Total data points: $(nrow(data))")
    println("System sizes: $(sort(unique(data.L)))")
    println("Lambda values: $(sort(unique(data.lambda)))")
    println("Seeds per (L,λ): $(sort(unique(stats_data.n_seeds)))")
    println("Trials per job: $(sort(unique(data.ntrials)))")
    
    println("\n" * "="^60)
    println("BINDER PARAMETER RESULTS")
    println("="^60)
    println(@sprintf("%-3s %-6s %-10s %-10s %-10s %-6s", "L", "λ", "B_mean", "B_std", "B_sem", "n"))
    println("-"^60)
    
    for row in eachrow(sort(stats_data, [:L, :lambda]))
        println(@sprintf("%-3d %-6.3f %-10.6f %-10.6f %-10.6f %-6d", 
                row.L, row.lambda, row.B_mean, row.B_std, row.B_sem, row.n_seeds))
    end
end

function main()
    println("="^60)
    println("EDWARDS-ANDERSON BINDER PARAMETER ANALYSIS")
    println("JSON Format Results")
    println("="^60)
    
    # Load and process data
    all_data = load_json_results()
    stats_data = compute_statistics(all_data)
    
    # Print summary
    print_data_summary(all_data, stats_data)
    
    # Create plots if we have data
    if nrow(stats_data) > 0
        println("\nCreating plots...")
        plots = plot_available_data(stats_data)
        
        if plots isa Tuple
            p1, p2 = plots
            savefig(p1, "binder_parameter_results.png")
            if p2 !== nothing
                savefig(p2, "binder_critical_region.png")
                println("Plots saved: binder_parameter_results.png, binder_critical_region.png")
            else
                println("Plot saved: binder_parameter_results.png")
            end
        else
            savefig(plots, "binder_parameter_results.png")
            println("Plot saved: binder_parameter_results.png")
        end
        
        # Save processed data
        CSV.write("binder_statistics.csv", stats_data)
        println("Statistics saved to: binder_statistics.csv")
    end
    
    println("\n" * "="^60)
    println("ANALYSIS COMPLETE!")
    println("="^60)
    
    return all_data, stats_data
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
