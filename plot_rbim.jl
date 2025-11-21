#!/usr/bin/env julia

"""
Plot RBIM phase diagram
Shows the Random Bond Ising Model transition at different dephasing strengths
X-axis: λ_zz (ZZ measurement strength)
Y-axis: Binder parameter
Different curves for different P_x (dephasing strengths)
"""

using JSON
using Plots
using Statistics
using CSV
using DataFrames

println("="^60)
println("RBIM Phase Diagram Analysis")
println("Strong dephasing limit: λ_x = 0")
println("Order parameter: EA correlations [⟨Z_i Z_j⟩²]_J")
println("="^60)

# Load results
function load_rbim_results(result_dirs)
    all_results = Dict{Tuple{Int,Float64}, Vector{Dict{String, Any}}}()
    
    for dir in result_dirs
        if !isdir(dir)
            continue
        end
        
        files = filter(f -> endswith(f, ".json") && startswith(basename(f), "rbim_"), 
                      readdir(dir, join=true))
        
        if !isempty(files)
            println("Found $(length(files)) files in $dir")
        end
        
        for file in files
            try
                data = JSON.parsefile(file)
                L = data["L"]
                P_x = data["P_x"]
                key = (L, P_x)
                
                if !haskey(all_results, key)
                    all_results[key] = []
                end
                push!(all_results[key], data)
            catch e
                println("Error reading $file: $e")
            end
        end
    end
    
    return all_results
end

# Calculate averages
function calculate_rbim_averages(results_dict)
    binder_data = Dict{Tuple{Int,Float64}, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
    
    for ((L, P_x), results) in results_dict
        # Group by lambda_zz
        lambda_groups = Dict{Float64, Vector{Float64}}()
        
        for data in results
            lambda_zz = data["lambda_zz"]
            binder = data["binder_parameter"]
            
            if !haskey(lambda_groups, lambda_zz)
                lambda_groups[lambda_zz] = []
            end
            push!(lambda_groups[lambda_zz], binder)
        end
        
        # Calculate averages
        lambdas = Float64[]
        binders = Float64[]
        errors = Float64[]
        
        for (lambda_zz, binder_list) in sort(collect(lambda_groups))
            push!(lambdas, lambda_zz)
            push!(binders, mean(binder_list))
            push!(errors, std(binder_list) / sqrt(length(binder_list)))
        end
        
        binder_data[(L, P_x)] = (lambdas, binders, errors)
        println("L=$L, P=$P_x: $(length(lambdas)) λ_zz points, $(sum(length(v) for v in values(lambda_groups))) total samples")
    end
    
    return binder_data
end

# Search for result directories
result_dirs = ["rbim_results", "jobs/results_rbim"]
for entry in readdir(".", join=false)
    if isdir(entry) && startswith(entry, "rbim_results_")
        push!(result_dirs, entry)
    end
end

println("\nSearching directories:")
for dir in result_dirs
    if isdir(dir)
        println("  - $dir")
    end
end
println()

# Load and process
all_results = load_rbim_results(result_dirs)
println("\nLoaded results for (L, P_x) combinations: ", sort(collect(keys(all_results))))

binder_data = calculate_rbim_averages(all_results)

# Get unique L and P values
L_values = sort(unique([k[1] for k in keys(binder_data)]))
P_values = sort(unique([k[2] for k in keys(binder_data)]))

println("\nSystem sizes: $L_values")
println("Dephasing strengths: $P_values")

# Create plots: Binder vs P at each λ_zz (looking for Nishimori point)
lambda_zz_values = sort(unique([k[2] for (L, P) in keys(binder_data) for k in [(P, get(binder_data, (L, P), ([],))[1])] if !isempty(k[2])]))

for lambda_zz in lambda_zz_values
    println("\nGenerating plot for λ_zz=$lambda_zz...")
    
    p = plot(
        xlabel="P (Dephasing Probability)",
        ylabel="EA Binder Parameter",
        title="RBIM Nishimori Point Search (λ_zz=$lambda_zz, λ_x=0)",
        legend=:best,
        grid=true,
        size=(800, 600),
        dpi=300,
        framestyle=:box
    )
    
    colors = [:blue, :red, :green, :purple, :orange, :brown, :pink, :gray]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5, :hexagon, :cross]
    
    for (idx, L) in enumerate(L_values)
        P_vals = Float64[]
        B_vals = Float64[]
        E_vals = Float64[]
        
        for P_x in P_values
            if haskey(binder_data, (L, P_x))
                lambdas, binders, errors = binder_data[(L, P_x)]
                # Find the index for this lambda_zz
                lambda_idx = findfirst(x -> isapprox(x, lambda_zz, atol=0.01), lambdas)
                if !isnothing(lambda_idx)
                    push!(P_vals, P_x)
                    push!(B_vals, binders[lambda_idx])
                    push!(E_vals, errors[lambda_idx])
                end
            end
        end
        
        if !isempty(P_vals)
            plot!(p, P_vals, B_vals,
                  yerr=E_vals,
                  marker=markers[mod1(idx, length(markers))],
                  markersize=5,
                  color=colors[mod1(idx, length(colors))],
                  linewidth=2,
                  label="L=$L")
        end
    end
    
    output_file = "plot_rbim_lzz$(replace(string(lambda_zz), "." => "p")).png"
    savefig(p, output_file)
    println("  Saved: $output_file")
end

# Save summary CSV
println("\nSaving summary statistics...")
df_rows = []
for ((L, P_x), (lambdas, binders, errors)) in binder_data
    for i in 1:length(lambdas)
        push!(df_rows, (L=L, P_x=P_x, lambda_zz=lambdas[i], binder=binders[i], error=errors[i]))
    end
end

df = DataFrame(df_rows)
csv_file = "rbim_summary.csv"
CSV.write(csv_file, df)
println("Summary saved to: $csv_file")

println("\n" * "="^60)
println("RBIM analysis complete!")
println("="^60)
