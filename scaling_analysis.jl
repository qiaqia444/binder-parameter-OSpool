using JSON
using Statistics
using Plots
using LsqFit
using Printf

default(fontfamily="Times", 
        linewidth=2,
        gridwidth=1,
        framestyle=:box,
        size=(800, 600),
        dpi=300)


function load_scaling_data_lambda05()
    results_dir = "standard_results_20251003_1508"
    scaling_data = Dict{Int, Any}()
    
    # Load L=8,12,16 data from existing results directory (if available)
    for L in [8, 12, 16]
        if isdir(results_dir)
            files = filter(f -> contains(f, "standard_L$(L)_lam0.5_s"), readdir(results_dir))
            
            x2_values = Float64[]
            x4_values = Float64[]
            
            for file in files
                try
                    filepath = joinpath(results_dir, file)
                    data = JSON.parsefile(filepath)
                    
                    if haskey(data, "S2_bar") && haskey(data, "S4_bar")
                        # x₂ = Σᵢⱼ |⟨ZᵢZⱼ⟩|² = S2_bar × L²
                        # x₄ = Σᵢⱼₖₗ |⟨ZᵢZⱼZₖZₗ⟩|² = S4_bar × L⁴
                        x2 = data["S2_bar"] * L^2
                        x4 = data["S4_bar"] * L^4
                        
                        push!(x2_values, x2)
                        push!(x4_values, x4)
                    end
                catch e
                    continue
                end
            end
            
            if !isempty(x2_values)
                scaling_data[L] = Dict(
                    "x2_mean" => mean(x2_values),
                    "x2_std" => std(x2_values),
                    "x4_mean" => mean(x4_values),
                    "x4_std" => std(x4_values),
                    "n_samples" => length(x2_values)
                )
            end
        end
    end
    
    # Add L=20 data directly from cluster output
    L20_data = [
        (0.2997898911606759, 0.17709772525408954),  # (S2_bar, S4_bar) sample 1
        (0.3063767706105352, 0.18326162981229593),  # sample 2  
        (0.31127483886816226, 0.18721760505236631)  # sample 3
    ]
    
    L = 20
    x2_values_20 = [s2 * L^2 for (s2, s4) in L20_data]
    x4_values_20 = [s4 * L^4 for (s2, s4) in L20_data]
    
    scaling_data[L] = Dict(
        "x2_mean" => mean(x2_values_20),
        "x2_std" => std(x2_values_20),
        "x4_mean" => mean(x4_values_20), 
        "x4_std" => std(x4_values_20),
        "n_samples" => length(x2_values_20)
    )
    
    println("Data loaded:")
    for size in sort(collect(keys(scaling_data)))
        d = scaling_data[size]
        println("L=$size: x2=$(round(d["x2_mean"], digits=2))±$(round(d["x2_std"], digits=2)), x4=$(round(d["x4_mean"], digits=2))±$(round(d["x4_std"], digits=2)) (n=$(d["n_samples"]))")
    end
    
    return scaling_data
end

scaling_data = load_scaling_data_lambda05()

system_sizes = sort(collect(keys(scaling_data)))
x2_means = [scaling_data[L]["x2_mean"] for L in system_sizes]
x4_means = [scaling_data[L]["x4_mean"] for L in system_sizes]
x2_stds = [scaling_data[L]["x2_std"] for L in system_sizes]
x4_stds = [scaling_data[L]["x4_std"] for L in system_sizes]

# Initialize result variables
α2 = A2 = r2_x2 = α4 = A4 = r2_x4 = ratio = NaN
hypothesis_consistent = false

# Perform power law fits
if length(system_sizes) >= 2
    # Convert to log scale for linear fitting
    log_L = log.(system_sizes)
    log_x2 = log.(x2_means)
    log_x4 = log.(x4_means)
    
    # Linear model: log(y) = log(A) + α*log(L)
    linear_model(x, p) = p[1] .+ p[2] .* x
    
    # Fit x₂ scaling: x₂ = A₂ * L^α₂
    try
        fit_x2 = curve_fit(linear_model, log_L, log_x2, [0.0, 1.0])
        log_A2 = fit_x2.param[1]
        α2 = fit_x2.param[2]
        A2 = exp(log_A2)
        
        # Calculate R²
        y_pred_x2 = log_A2 .+ α2 .* log_L
        ss_res_x2 = sum((log_x2 .- y_pred_x2).^2)
        ss_tot_x2 = sum((log_x2 .- mean(log_x2)).^2)
        r2_x2 = 1 - ss_res_x2 / ss_tot_x2
    catch e
        α2 = A2 = r2_x2 = NaN
    end
    
    # Fit x₄ scaling: x₄ = A₄ * L^α₄
    try
        fit_x4 = curve_fit(linear_model, log_L, log_x4, [0.0, 2.0])
        log_A4 = fit_x4.param[1]
        α4 = fit_x4.param[2]
        A4 = exp(log_A4)
        
        # Calculate R²
        y_pred_x4 = log_A4 .+ α4 .* log_L
        ss_res_x4 = sum((log_x4 .- y_pred_x4).^2)
        ss_tot_x4 = sum((log_x4 .- mean(log_x4)).^2)
        r2_x4 = 1 - ss_res_x4 / ss_tot_x4
    catch e
        α4 = A4 = r2_x4 = NaN
    end
    
    # Calculate ratio and test hypothesis
    if !isnan(α2) && !isnan(α4)
        ratio = α4 / α2
        hypothesis_consistent = abs(ratio - 2.0) < 0.15
    end
    
    # Create improved log-log plot with extended range including L=20
    p = plot(title="Scaling at Critical Point λ = 0.5\\nPower Law Analysis (Including L=20 Data)",
             xlabel="System Size L", ylabel="Correlator Sums x2, x4",
             xscale=:log10, yscale=:log10,
             grid=true, gridcolor=:lightgray, gridalpha=0.5,
             margin=5Plots.mm, legendfontsize=9)
    
    # Plot data points with error bars using consistent colors
    scatter!(p, system_sizes, x2_means, yerror=x2_stds,
             label="x₂ = S₂×L² (2nd moment)", color=:blue, markersize=7, 
             markerstrokewidth=2, alpha=0.9, markershape=:circle)
    
    scatter!(p, system_sizes, x4_means, yerror=x4_stds,
             label="x₄ = S₄×L⁴ (4th moment)", color=:red, markersize=7,
             markerstrokewidth=2, alpha=0.9, markershape=:square)
    
    # Plot fitted lines with extended range
    if !isnan(α2) && !isnan(A2)
        L_min = minimum(system_sizes) * 0.8
        L_max = maximum(system_sizes) * 1.3
        L_fit = exp.(range(log(L_min), log(L_max), length=100))
        plot!(p, L_fit, A2 .* L_fit.^α2, color=:blue, linestyle=:dash, linewidth=2,
              label="x₂ ~ L^$(round(α2, digits=2)) (R²=$(round(r2_x2, digits=3)))", alpha=0.8)
    end
    
    if !isnan(α4) && !isnan(A4)
        L_min = minimum(system_sizes) * 0.8
        L_max = maximum(system_sizes) * 1.3
        L_fit = exp.(range(log(L_min), log(L_max), length=100))
        plot!(p, L_fit, A4 .* L_fit.^α4, color=:red, linestyle=:dash, linewidth=2,
              label="x₄ ~ L^$(round(α4, digits=2)) (R²=$(round(r2_x4, digits=3)))", alpha=0.8)
    end
    
    # Add improved annotation showing theoretical expectation
    if !isnan(ratio)
        y_pos = exp(0.7 * log(maximum(x4_means)) + 0.3 * log(minimum(x4_means)))
        x_pos = exp(0.3 * log(maximum(system_sizes)) + 0.7 * log(minimum(system_sizes)))
        annotate!(p, [(x_pos, y_pos, 
                      text("Ratio α₄/α₂ = $(round(ratio, digits=3))\\nTheory: 2.0", 9, :black, :center))])
    end
    
    # Set consistent legend and font sizes
    plot!(legendfontsize=10, guidefontsize=12, titlefontsize=11)
    
    savefig(p, "scaling_analysis_lambda05_with_L20.png")
    
    # Print comprehensive analysis summary
    println("\\n" * "="^60)
    println("SCALING ANALYSIS SUMMARY (λ = 0.5 Critical Point)")
    println("="^60)
    println("System sizes analyzed: $(system_sizes)")
    println("Total data points: $(sum([scaling_data[L]["n_samples"] for L in system_sizes]))")
    println()
    
    if !isnan(α2) && !isnan(α4)
        println("POWER LAW SCALING RESULTS:")
        println("• x₂ scaling exponent: α₂ = $(round(α2, digits=3)) ± $(round(sqrt(abs(r2_x2)), digits=3))")
        println("• x₄ scaling exponent: α₄ = $(round(α4, digits=3)) ± $(round(sqrt(abs(r2_x4)), digits=3))")
        println("• Ratio α₄/α₂ = $(round(ratio, digits=3))")
        println("• Theoretical expectation: α₄/α₂ = 2.0")
        println("• Deviation: $(round(abs(ratio - 2.0), digits=3)) $(abs(ratio - 2.0) < 0.1 ? "✓ Excellent!" : abs(ratio - 2.0) < 0.2 ? "✓ Good" : "⚠ Significant")")
        println("• Fit quality: R²(x₂) = $(round(r2_x2, digits=3)), R²(x₄) = $(round(r2_x4, digits=3))")
        
        println("\\nDATA QUALITY ASSESSMENT:")
        if length(system_sizes) >= 4
            println("• Excellent: ≥4 system sizes for robust scaling")
        elseif length(system_sizes) >= 3
            println("• Good: 3 system sizes provide reasonable scaling estimate")
        else
            println("• Limited: Only $(length(system_sizes)) sizes - need more data for robust analysis")
        end
        
        if abs(ratio - 2.0) < 0.15
            println("• Theory consistency: PASSED (ratio within 15% of expectation)")
        else
            println("• Theory consistency: MARGINAL (ratio deviates >15% from theory)")
        end
    else
        println("SCALING FIT FAILED - insufficient data or numerical issues")
    end
    
    println("\\nOUTPUT FILES:")
    println("• Plot: scaling_analysis_lambda05_with_L20.png")
    println("\\nWith L=24,28,32,36 data from your critical jobs, the scaling analysis will be much more robust!")
    
else
    println("ERROR: Need at least 2 system sizes for power law fitting")
end