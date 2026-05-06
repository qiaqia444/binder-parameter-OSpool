#!/usr/bin/env julia

"""
Analyze Learning-to-Trivial Transition Results (lambda=0.7 family)

Parameters: λ_x = 0.49, λ_zz = 0.21
"""

using JSON
using DataFrames
using Statistics
using Printf

# ====================================================================
# Configuration
# ====================================================================
const DEFAULT_RESULTS_DIR = "learning_to_trivial_lambda0.7_results"

function load_results(results_dir::String=DEFAULT_RESULTS_DIR)
    """Load all JSON results from directory"""
    
    if !isdir(results_dir)
        @error "Results directory not found: $results_dir"
        return DataFrame()
    end
    
    # Find all JSON files
    json_files = filter(f -> endswith(f, ".json"), readdir(results_dir; join=true))
    
    if isempty(json_files)
        @error "No JSON files found in $results_dir"
        return DataFrame()
    end
    
    # Load data
    data = []
    for file in json_files
        try
            content = read(file, String)
            json_data = JSON.parse(content)
            
            # Handle both single results and arrays
            results = isa(json_data, Array) ? json_data : [json_data]
            
            for result in results
                push!(data, result)
            end
        catch e
            @warn "Failed to load $file: $e"
        end
    end
    
    # Convert to DataFrame
    if isempty(data)
        return DataFrame()
    end
    
    df = DataFrame()
    for (i, result) in enumerate(data)
        for (key, value) in result
            if !haskey(df, key)
                df[!, Symbol(key)] = Any[]
            end
        end
    end
    
    # Fill data
    for result in data
        row = Dict()
        for col in names(df)
            row[col] = get(result, String(col), missing)
        end
        push!(df, row)
    end
    
    return df
end

function analyze_by_parameter(df::DataFrame)
    """Group results by parameters and compute statistics"""
    
    if nrow(df) == 0
        println("No data to analyze")
        return
    end
    
    # Group by (L, P_x, P_zz)
    grouped = combine(
        groupby(df, [:L, :P_x]),
        :B => mean => :B_mean,
        :B => std => :B_std,
        :B_mean_of_trials => mean => :B_trial_mean,
        :B_std_of_trials => mean => :B_trial_std_mean,
        :M2_bar => mean => :M2_mean,
        :M4_bar => mean => :M4_mean,
        :purity_bar => mean => :purity_mean,
        :purity_bar => std => :purity_std,
        nrow => :n_results
    )
    
    return sort(grouped, [:L, :P_x])
end

function print_analysis(df::DataFrame)
    """Print analysis summary"""
    
    println("="^80)
    println("LEARNING-TO-TRIVIAL TRANSITION ANALYSIS")
    println("lambda=0.7 family: λ_x = 0.49, λ_zz = 0.21")
    println("="^80)
    println()
    
    # Overall statistics
    println("Overall Statistics:")
    println("  Total results: $(nrow(df))")
    println("  System sizes: $(sort(unique(df.L)))")
    println("  P_x values: $(sort(unique(df.P_x)))")
    println()
    
    # Analyze by system size
    for L in sort(unique(df.L))
        df_L = filter(row -> row.L == L, df)
        println("L = $L:")
        
        for P in sort(unique(df_L.P_x))
            df_LP = filter(row -> row.P_x == P, df_L)
            B_mean = mean(df_LP.B)
            B_std = std(df_LP.B)
            purity_mean = mean(df_LP.purity_bar)
            
            @printf("  P_x = %.2f: B = %.4f ± %.4f, purity = %.4f\n", 
                    P, B_mean, B_std, purity_mean)
        end
        println()
    end
    
    # Phase transition indicators
    println("Phase Transition Indicators:")
    println()
    
    # Check Binder crossing at critical P for different L
    for L in sort(unique(df.L))
        df_L = filter(row -> row.L == L, df)
        Ps = sort(unique(df_L.P_x))
        Bs = [mean(filter(row -> row.P_x == P, df_L).B) for P in Ps]
        
        @printf("L = %2d: B varies from %.4f (P=0) to %.4f (P=0.5)\n", 
                L, Bs[1], Bs[end])
    end
    println()
    
    println("="^80)
end

function main()
    # Load results
    df = load_results()
    
    if nrow(df) == 0
        println("ERROR: Could not load results. Check directory: $DEFAULT_RESULTS_DIR")
        return
    end
    
    # Analyze
    print_analysis(df)
    
    # Save summary
    summary_df = analyze_by_parameter(df)
    summary_file = "learning_to_trivial_lambda0.7_analysis_summary.csv"
    CSV.write(summary_file, summary_df)
    println("Summary saved to: $summary_file")
end

main()
