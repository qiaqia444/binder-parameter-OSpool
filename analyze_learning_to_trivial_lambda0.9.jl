#!/usr/bin/env julia

"""
Analyze Learning-to-Trivial Transition Results (lambda=0.9 family)

Parameters: λ_x = 0.63, λ_zz = 0.07
"""

using JSON
using CSV
using DataFrames
using Statistics
using Printf
using Plots
using LaTeXStrings

# ====================================================================
# Configuration
# ====================================================================
const DEFAULT_RESULTS_DIR = "learning_to_trivial_lambda0.9_results"

function resolve_results_dir(results_dir::String=DEFAULT_RESULTS_DIR)
    if isdir(results_dir)
        return results_dir
    end

    candidates = filter(
        name -> startswith(name, "learning_to_trivial_lambda0.9_results_") && isdir(name),
        readdir(pwd())
    )
    if !isempty(candidates)
        return joinpath(pwd(), sort(candidates)[end])
    end

    return results_dir
end

function load_results(results_dir::String=DEFAULT_RESULTS_DIR)
    """Load all JSON results from directory"""

    results_dir = resolve_results_dir(results_dir)
    
    if !isdir(results_dir)
        @error "Results directory not found: $results_dir"
        return DataFrame()
    end
    
    # Find all JSON files recursively (results are stored in L-specific subdirectories)
    json_files = String[]
    for (root, _, files) in walkdir(results_dir)
        for file in files
            if endswith(file, ".json")
                push!(json_files, joinpath(root, file))
            end
        end
    end
    
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
            if !(Symbol(key) in names(df))
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
        :B => (x -> std(x) / sqrt(length(x))) => :B_sem,
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
    println("lambda=0.9 family: λ_x = 0.63, λ_zz = 0.07")
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

function save_plot(summary_df::DataFrame)
    if nrow(summary_df) == 0
        return
    end

    gr()
    p = plot(
        xlabel = L"P_x = P_{zz}" * " (dephasing probability)",
        ylabel = "Binder Parameter",
        title = L"\lambda_x = 0.63, \lambda_{zz} = 0.07",
        legend = :topright,
        grid = true,
        size = (800, 600),
        dpi = 300,
        framestyle = :box,
    )

    colors = [:blue, :red, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle]

    for L in sort(unique(summary_df.L))
        df_L = filter(row -> row.L == L, summary_df)
        sort!(df_L, :P_x)

        color = colors[(mod1(findfirst(==(L), sort(unique(summary_df.L))), length(colors)))]
        marker = markers[(mod1(findfirst(==(L), sort(unique(summary_df.L))), length(markers)))]

        plot!(p, df_L.P_x, df_L.B_mean,
              yerr = df_L.B_sem,
              label = "L = $L",
              color = color,
              marker = marker,
              markersize = 6,
              linewidth = 2)

    end

    hline!(p, [2/3], linestyle = :dash, color = :black, label = "B = 2/3 (pure)", linewidth = 2, alpha = 0.7)

    png_file = "learning_to_trivial_lambda0.9_binder_vs_px.png"
    pdf_file = "learning_to_trivial_lambda0.9_binder_vs_px.pdf"
    savefig(p, png_file)
    savefig(p, pdf_file)
    println("Plot saved to: $png_file")
    println("Plot saved to: $pdf_file")
end

function main()
    results_dir = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_RESULTS_DIR

    # Load results
    df = load_results(results_dir)
    
    if nrow(df) == 0
        println("ERROR: Could not load results. Check directory: $(resolve_results_dir(results_dir))")
        return
    end
    
    # Analyze
    print_analysis(df)
    
    # Save summary
    summary_df = analyze_by_parameter(df)
    summary_file = "learning_to_trivial_lambda0.9_analysis_summary.csv"
    CSV.write(summary_file, summary_df)
    println("Summary saved to: $summary_file")

    save_plot(summary_df)
end

main()
