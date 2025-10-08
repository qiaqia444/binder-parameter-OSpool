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
    
    for L in [8, 12, 16]
        # Find all λ = 0.5 files for this system size
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
    
    # Create log-log plot with ASCII labels to avoid font issues
    p = plot(title="Scaling at Critical Point λ = 0.5\\nPower Law Analysis",
             xlabel="System Size L", ylabel="Correlator Sums x2, x4",
             xscale=:log10, yscale=:log10,
             grid=true, gridcolor=:lightgray, gridalpha=0.5,
             margin=5Plots.mm)
    
    # Plot data points with error bars using consistent colors
    scatter!(p, system_sizes, x2_means, yerror=x2_stds,
             label="x2 = sum_ij <Zi*Zj>^2", color=:blue, markersize=6, 
             markerstrokewidth=2, alpha=0.8)
    
    scatter!(p, system_sizes, x4_means, yerror=x4_stds,
             label="x4 = sum_ijkl <Zi*Zj*Zk*Zl>^2", color=:red, markersize=6,
             markerstrokewidth=2, alpha=0.8)
    
    # Plot fitted lines if successful
    if !isnan(α2) && !isnan(A2)
        L_fit = 7.5:0.05:16.5
        plot!(p, L_fit, A2 .* L_fit.^α2, color=:blue, linestyle=:dash, linewidth=3,
              label="x2 ~ L^$(round(α2, digits=2)) (R²=$(round(r2_x2, digits=3)))", alpha=0.8)
    end
    
    if !isnan(α4) && !isnan(A4)
        L_fit = 7.5:0.05:16.5
        plot!(p, L_fit, A4 .* L_fit.^α4, color=:red, linestyle=:dash, linewidth=3,
              label="x4 ~ L^$(round(α4, digits=2)) (R²=$(round(r2_x4, digits=3)))", alpha=0.8)
    end
    
    # Add ratio annotation with ASCII
    if !isnan(ratio)
        annotate!(p, [(14, maximum(x4_means)*0.3, 
                      text("α4/α2 = $(round(ratio, digits=2))", 12, :black, :center))])
    end
    
    # Set consistent legend and font sizes
    plot!(legendfontsize=10, guidefontsize=12, titlefontsize=11)
    
    savefig(p, "scaling_analysis_lambda05_clean.png")
    
else
    println("ERROR: Need at least 2 system sizes for power law fitting")
end