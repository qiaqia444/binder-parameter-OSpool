#!/usr/bin/env julia

"""
Analyze Left Boundary Scan Results

Load all JSON files, compute statistics, and prepare data for plotting.
"""

using JSON
using Statistics
using DataFrames
using Printf
using CSV

function load_results(results_dir)
    """Load all JSON files from the results directory."""
    all_data = []
    
    for L in [8, 10, 12, 14, 16]
        L_dir = joinpath(results_dir, "L$L")
        if !isdir(L_dir)
            println("Warning: Directory $L_dir not found")
            continue
        end
        
        json_files = filter(f -> endswith(f, ".json"), readdir(L_dir, join=true))
        println("Loading L=$L: $(length(json_files)) files")
        
        for file in json_files
            try
                data = JSON.parsefile(file)
                # Handle both single dict and array format
                if data isa Array
                    for entry in data
                        push!(all_data, entry)
                    end
                else
                    push!(all_data, data)
                end
            catch e
                println("Warning: Failed to load $file: $e")
            end
        end
    end
    
    return DataFrame(all_data)
end

function compute_statistics(df)
    """Compute mean and std of B for each (L, P_x) combination."""
    # Group by L and P_x
    gdf = groupby(df, [:L, :P_x])
    
    stats = combine(gdf) do group
        (
            B_mean = mean(group.B),
            B_std = std(group.B),
            B_sem = std(group.B) / sqrt(nrow(group)),
            M2_mean = mean(group.S2_bar),
            M4_mean = mean(group.S4_bar),
            n_samples = nrow(group),
            time_mean = mean(group.time_seconds)
        )
    end
    
    sort!(stats, [:L, :P_x])
    return stats
end

function print_summary(stats)
    """Print a nice summary table."""
    println("\n" * "="^80)
    println("LEFT BOUNDARY SCAN RESULTS SUMMARY")
    println("="^80)
    
    for L in unique(stats.L)
        L_data = filter(row -> row.L == L, stats)
        println("\nL = $L ($(nrow(L_data)) P_x values, $(L_data.n_samples[1]) samples each)")
        println("-"^80)
        println("  P_x    |   B_mean  |  B_std   |  B_sem   |   M₂²    |   M₄²    | Time(s)")
        println("-"^80)
        
        for row in eachrow(L_data)
            @printf("  %.2f   |  %.4f   |  %.4f  |  %.4f  | %.6f | %.6f |  %.1f\n",
                    row.P_x, row.B_mean, row.B_std, row.B_sem, 
                    row.M2_mean, row.M4_mean, row.time_mean)
        end
    end
    
    println("\n" * "="^80)
end

function check_phase_transition(stats)
    """Identify potential phase transition region."""
    println("\n" * "="^80)
    println("PHASE TRANSITION ANALYSIS")
    println("="^80)
    
    for L in unique(stats.L)
        L_data = filter(row -> row.L == L, stats)
        
        # Find where B drops from ~2/3 to ~0
        high_B = filter(row -> row.B_mean > 0.5, L_data)
        low_B = filter(row -> row.B_mean < 0.3, L_data)
        
        if nrow(high_B) > 0 && nrow(low_B) > 0
            P_high = maximum(high_B.P_x)
            P_low = minimum(low_B.P_x)
            println("\nL = $L:")
            println("  Volume-law phase (B ≈ 2/3): P_x ≤ $(P_high)")
            println("  Area-law phase (B ≈ 0):     P_x ≥ $(P_low)")
            println("  Transition region:          $(P_high) < P_x < $(P_low)")
            
            # Find steepest drop
            if nrow(L_data) > 1
                diffs = diff(L_data.B_mean)
                max_drop_idx = argmin(diffs)
                println("  Steepest drop at:           P_x ≈ $(L_data.P_x[max_drop_idx])")
            end
        else
            println("\nL = $L: No clear phase transition in this range")
        end
    end
    
    println("\n" * "="^80)
end

# Main execution
function main()
    # Find results directory
    results_dirs = filter(d -> startswith(basename(d), "left_boundary_results_") && isdir(d), 
                          readdir(".", join=true))
    
    if isempty(results_dirs)
        println("ERROR: No results directory found!")
        println("Expected directory name: left_boundary_results_<timestamp>")
        return
    end
    
    results_dir = results_dirs[end]  # Use most recent
    println("Analyzing results from: $results_dir")
    
    # Load data
    println("\nLoading JSON files...")
    df = load_results(results_dir)
    println("Total data points loaded: $(nrow(df))")
    
    # Compute statistics
    println("\nComputing statistics...")
    stats = compute_statistics(df)
    
    # Save processed data
    output_file = "left_boundary_analysis.csv"
    CSV.write(output_file, stats)
    println("\nProcessed data saved to: $output_file")
    
    # Print summary
    print_summary(stats)
    
    # Check for phase transition
    check_phase_transition(stats)
    
    println("\nAnalysis complete!")
    println("Next step: julia plot_left_boundary.jl")
end

main()
