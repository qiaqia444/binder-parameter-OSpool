#!/usr/bin/env julia

using DataFrames, CSV, Statistics, Plots, StatsPlots
using Printf

"""
Comprehensive analysis of L=8,10,12 Binder parameter results
Analyzes phase transition with 2000 trials per job and 20 seeds per (L,λ)
"""

function load_all_results(output_dir="output")
    println("Loading results from $output_dir...")
    
    all_data = DataFrame(L=Int[], lambda=Float64[], lambda_x=Float64[], lambda_zz=Float64[], 
                        B=Float64[], B_mean=Float64[], B_std=Float64[], 
                        S2_bar=Float64[], S4_bar=Float64[], ntrials=Int[], 
                        seed=Int[], job_type=String[])
    
    csv_files = filter(f -> endswith(f, ".csv"), readdir(output_dir))
    println("Found $(length(csv_files)) CSV files")
    
    for file in csv_files
        filepath = joinpath(output_dir, file)
        try
            data = CSV.read(filepath, DataFrame)
            if nrow(data) > 0
                # Extract job type from filename
                job_type = startswith(file, "test_") ? "test" : "standard"
                data.job_type .= job_type
                append!(all_data, data)
            end
        catch e
            println("Warning: Could not read $file: $e")
        end
    end
    
    println("Loaded $(nrow(all_data)) data points")
    println("System sizes: $(sort(unique(all_data.L)))")
    println("Lambda range: $(round(minimum(all_data.lambda), digits=3)) - $(round(maximum(all_data.lambda), digits=3))")
    
    return all_data
end

function compute_statistics(data)
    """Compute mean and std for each (L, lambda) combination"""
    grouped = groupby(data, [:L, :lambda])
    
    stats = combine(grouped) do df
        DataFrame(
            lambda_x = df.lambda_x[1],
            lambda_zz = df.lambda_zz[1],
            B_mean = mean(df.B),
            B_std = std(df.B),
            B_sem = std(df.B) / sqrt(length(df.B)),  # Standard error
            n_seeds = length(df.B),
            S2_mean = mean(df.S2_bar),
            S4_mean = mean(df.S4_bar)
        )
    end
    
    return sort(stats, [:L, :lambda])
end

function plot_phase_transition(stats_data)
    """Create comprehensive phase transition plots"""
    
    # Main Binder parameter plot
    p1 = plot(title="Edwards-Anderson Binder Parameter vs λ", 
              xlabel="λ (λₓ=λ, λ_zz=1-λ)", ylabel="Binder Parameter B",
              grid=true, size=(800, 600), dpi=300)
    
    colors = [:red, :blue, :green, :purple]
    L_values = sort(unique(stats_data.L))
    
    for (i, L) in enumerate(L_values)
        L_data = filter(row -> row.L == L, stats_data)
        color = colors[i]
        
        # Plot with error bars
        plot!(p1, L_data.lambda, L_data.B_mean,
              yerror=L_data.B_sem,
              marker=:circle, markersize=4, linewidth=2,
              label="L = $L", color=color)
    end
    
    # Add critical region shading
    vspan!(p1, [0.45, 0.55], alpha=0.2, color=:gray, label="Critical Region")
    
    # Zoom plot around critical region
    p2 = plot(title="Critical Region Detail (λ = 0.4-0.6)", 
              xlabel="λ", ylabel="Binder Parameter B",
              grid=true, size=(600, 400), dpi=300,
              xlims=(0.4, 0.6))
    
    for (i, L) in enumerate(L_values)
        L_data = filter(row -> row.L == L && 0.4 <= row.lambda <= 0.6, stats_data)
        if nrow(L_data) > 0
            color = colors[i]
            plot!(p2, L_data.lambda, L_data.B_mean,
                  yerror=L_data.B_sem,
                  marker=:circle, markersize=6, linewidth=3,
                  label="L = $L", color=color)
        end
    end
    
    # Finite-size scaling plot
    p3 = plot(title="Finite-Size Scaling at Critical Region", 
              xlabel="1/L", ylabel="Binder Parameter B",
              grid=true, size=(600, 400), dpi=300)
    
    critical_lambdas = [0.47, 0.48, 0.49, 0.50, 0.51]
    for (i, λ) in enumerate(critical_lambdas)
        λ_data = filter(row -> abs(row.lambda - λ) < 0.005, stats_data)
        if nrow(λ_data) > 0
            inv_L = 1 ./ λ_data.L
            plot!(p3, inv_L, λ_data.B_mean,
                  yerror=λ_data.B_sem,
                  marker=:circle, markersize=5, linewidth=2,
                  label="λ = $λ")
        end
    end
    
    return p1, p2, p3
end

function find_critical_point(stats_data)
    """Estimate critical point from steepest slope"""
    println("\n" * "="^60)
    println("CRITICAL POINT ANALYSIS")
    println("="^60)
    
    for L in sort(unique(stats_data.L))
        L_data = filter(row -> row.L == L, stats_data)
        sort!(L_data, :lambda)
        
        # Find steepest slope
        max_slope = 0.0
        critical_lambda = 0.0
        
        for i in 1:(nrow(L_data)-1)
            slope = abs(L_data.B_mean[i+1] - L_data.B_mean[i]) / 
                   (L_data.lambda[i+1] - L_data.lambda[i])
            if slope > max_slope
                max_slope = slope
                critical_lambda = (L_data.lambda[i] + L_data.lambda[i+1]) / 2
            end
        end
        
        println("L = $L: Critical λ ≈ $(round(critical_lambda, digits=4)) (max slope = $(round(max_slope, digits=3)))")
    end
end

function print_summary_table(stats_data)
    """Print a nice summary table"""
    println("\n" * "="^80)
    println("BINDER PARAMETER SUMMARY TABLE")
    println("="^80)
    println(@sprintf("%-3s %-6s %-10s %-10s %-10s %-6s", "L", "λ", "B_mean", "B_std", "B_sem", "n"))
    println("-"^80)
    
    for L in sort(unique(stats_data.L))
        L_data = filter(row -> row.L == L, stats_data)
        sort!(L_data, :lambda)
        
        for row in eachrow(L_data)
            println(@sprintf("%-3d %-6.3f %-10.6f %-10.6f %-10.6f %-6d", 
                    row.L, row.lambda, row.B_mean, row.B_std, row.B_sem, row.n_seeds))
        end
        println("-"^80)
    end
end

function main()
    println("="^60)
    println("EDWARDS-ANDERSON BINDER PARAMETER ANALYSIS")
    println("L = [8, 10, 12] with 2000 trials per job")
    println("="^60)
    
    # Load and process data
    all_data = load_all_results()
    stats_data = compute_statistics(all_data)
    
    # Create plots
    println("\nCreating plots...")
    p1, p2, p3 = plot_phase_transition(stats_data)
    
    # Create analysis_results directory if it doesn't exist
    mkpath("analysis_results")
    
    # Save plots
    savefig(p1, "analysis_results/binder_parameter_L8_10_12.png")
    savefig(p2, "analysis_results/binder_critical_region.png")
    savefig(p3, "analysis_results/finite_size_scaling.png")
    println("Plots saved: analysis_results/binder_parameter_L8_10_12.png, analysis_results/binder_critical_region.png, analysis_results/finite_size_scaling.png")
    
    # Analysis
    find_critical_point(stats_data)
    print_summary_table(stats_data)
    
    # Save processed data
    CSV.write("analysis_results/binder_statistics_L8_10_12.csv", stats_data)
    println("\nProcessed statistics saved to: analysis_results/binder_statistics_L8_10_12.csv")
    
    println("\n" * "="^60)
    println("ANALYSIS COMPLETE!")
    println("Check the generated plots and CSV file for detailed results.")
    println("="^60)
    
    return stats_data
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
