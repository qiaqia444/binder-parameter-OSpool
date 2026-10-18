#!/usr/bin/env julia

"""
Analyze Right Boundary dual EA Binder Scan Results with a GHZ initial state

Usage:
    julia analyze_dual_ea_right_boundary_ghz.jl [ghz_results_dir] [allup_dual_ea_results_dir] [min_L]

Loads every `dual_ea_right_boundary_ghz_L*_...json` under <ghz_results_dir>/L*/ (default: the
newest `dual_ea_right_boundary_ghz_results_*`), de-duplicates jobs with the same
(L, q, point_seed), pools the per-trajectory M2/M4/B arrays of all jobs at each (L, q), and
reports BOTH Binder conventions separately:

    B_mean_of_trials = mean over trajectories of B_m,      SE = std(B_m)/sqrt(n)
    B_ratio_of_means = 1 - mean(M4_m)/(3 mean(M2_m)^2),    SE = bootstrap over trajectories

If a results directory of the ALL-UP dual-EA scan (from collect_dual_ea_right_boundary_results.sh)
is given, its B_mean_of_trials is added as a comparison column / dotted curves. The GHZ jobs
use seeds from 120001 and the all-up jobs seeds from 9001, so the comparison is statistical,
not trajectory by trajectory. `min_L` (default 16) drops smaller sizes at load time.
No crossing is assumed or fitted.

Outputs: dual_ea_right_boundary_ghz_summary.csv and dual_ea_right_boundary_ghz_*_vs_q.{pdf,png}
"""

using JSON
using Statistics
using Printf
using Random
using DelimitedFiles
using Plots

binder(M2, M4) = (isfinite(M2) && M2 > 0) ? 1.0 - M4 / (3.0 * M2^2) : NaN

# Pool per-trajectory arrays at each (L, q) from `<dir>/L*/<prefix>*.json`.
function load_pooled(results_dir, prefix; min_L::Int=16)
    pooled = Dict{Tuple{Int,Float64},Dict{String,Any}}()
    seen = Set{Tuple{Int,Float64,Int}}()
    n_dup = 0

    for L_dir in filter(isdir, readdir(results_dir, join=true))
        for file in filter(f -> startswith(basename(f), prefix) && endswith(f, ".json") &&
                                !occursin("FAILED", f), readdir(L_dir, join=true))
            try
                parsed = JSON.parsefile(file)
                for entry in (parsed isa AbstractVector ? parsed : [parsed])
                    L = Int(entry["L"])
                    L < min_L && continue
                    q = round(Float64(entry["P_x"]); digits=6)
                    key = (L, q, Int(entry["point_seed"]))
                    if key in seen
                        n_dup += 1
                        continue
                    end
                    push!(seen, key)
                    g = get!(pooled, (L, q), Dict{String,Any}(
                        "M2" => Float64[], "M4" => Float64[], "B" => Float64[],
                        "n_files" => 0, "settings" => Set{Tuple}(),
                    ))
                    append!(g["M2"], Float64.(entry["M2_per_trajectory"]))
                    append!(g["M4"], Float64.(entry["M4_per_trajectory"]))
                    append!(g["B"], Float64.(entry["B_per_trajectory"]))
                    g["n_files"] += 1
                    push!(g["settings"], (entry["T_max"], entry["maxdim"], entry["cutoff"]))
                end
            catch e
                @warn "Could not read JSON file" file exception=(e, catch_backtrace())
            end
        end
    end

    n_dup > 0 && println("  ($results_dir) skipped $n_dup duplicate jobs (same L, q, point_seed)")
    return pooled
end

function summarize(pooled; nboot=1000)
    rows = []
    for ((L, q), g) in sort(collect(pooled); by=first)
        length(g["settings"]) > 1 &&
            @warn "Mixed (T_max, maxdim, cutoff) settings pooled at this point" L q settings=g["settings"]

        M2, M4, B = g["M2"], g["M4"], g["B"]
        n = length(B)
        rng = MersenneTwister(1000L + round(Int, 1000q))
        boot = [begin
                    idx = rand(rng, 1:n, n)
                    binder(mean(M2[idx]), mean(M4[idx]))
                end for _ in 1:nboot]

        push!(rows, (
            L=L, q=q, n_traj=n, n_files=g["n_files"],
            M2_mean=mean(M2), M4_mean=mean(M4),
            B_mean_of_trials=mean(B), B_mean_of_trials_se=std(B) / sqrt(n),
            B_ratio_of_means=binder(mean(M2), mean(M4)),
            B_ratio_of_means_se=std(boot),
        ))
    end
    return rows
end

function main()
    dirs = filter(d -> startswith(d, "dual_ea_right_boundary_ghz_results_") && isdir(d), readdir())
    results_dir = !isempty(ARGS) ? ARGS[1] : (isempty(dirs) ? "" : last(sort(dirs)))
    min_L = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 16
    if isempty(results_dir) || !isdir(results_dir)
        println("ERROR: No results directory found (pass one, or run collect_dual_ea_right_boundary_ghz_results.sh).")
        return
    end
    println("Loading GHZ-start dual-EA results from: $results_dir (sizes L >= $min_L)")

    pooled = load_pooled(results_dir, "dual_ea_right_boundary_ghz_"; min_L=min_L)
    if isempty(pooled)
        println("ERROR: No dual_ea_right_boundary_ghz_*.json files found under $results_dir/L*/")
        return
    end

    rows = summarize(pooled)

    allup = Dict{Tuple{Int,Float64},Float64}()
    if length(ARGS) >= 2 && isdir(ARGS[2])
        println("Loading all-up dual-EA results for comparison from: $(ARGS[2])")
        up_pooled = load_pooled(ARGS[2], "dual_ea_right_boundary_L"; min_L=min_L)   # excludes *_fine_*, *_ghz_*
        for ((L, q), g) in up_pooled
            allup[(L, q)] = mean(g["B"])
        end
    end

    header = ["L" "q" "n_traj" "n_files" "M2_mean" "M4_mean" "B_mean_of_trials" "B_mean_of_trials_se" "B_ratio_of_means" "B_ratio_of_means_se" "B_mean_allup"]
    table = Any[header]
    println()
    @printf("%4s %6s %7s %10s %10s %18s %18s %11s\n", "L", "q", "n_traj", "M2_mean", "M4_mean",
            "B_mean_of_trials", "B_ratio_of_means", "B_allup")
    for r in rows
        bu = get(allup, (r.L, r.q), NaN)
        push!(table, [r.L r.q r.n_traj r.n_files r.M2_mean r.M4_mean r.B_mean_of_trials r.B_mean_of_trials_se r.B_ratio_of_means r.B_ratio_of_means_se bu])
        @printf("%4d %6.3f %7d %10.5f %10.5f %10.5f ± %.5f %10.5f ± %.5f %10.5f\n", r.L, r.q, r.n_traj,
                r.M2_mean, r.M4_mean, r.B_mean_of_trials, r.B_mean_of_trials_se,
                r.B_ratio_of_means, r.B_ratio_of_means_se, bu)
    end
    writedlm("dual_ea_right_boundary_ghz_summary.csv", vcat(table...), ',')

    Ls = sort(unique(r.L for r in rows))
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5, :hexagon]

    for (name, field, se_field) in (
        ("B_mean_of_trials", :B_mean_of_trials, :B_mean_of_trials_se),
        ("B_ratio_of_means", :B_ratio_of_means, :B_ratio_of_means_se),
    )
        p = plot(xlabel="q = P_x = P_zz", ylabel="dual EA Binder ($name)",
                 title="Right boundary, GHZ start: λ_x = 0.7, λ_zz = 0", legend=:outertopright,
                 grid=true, size=(850, 600), dpi=300)
        for (idx, L) in enumerate(Ls)
            sub = filter(r -> r.L == L, rows)
            plot!(p, [r.q for r in sub], [getfield(r, field) for r in sub];
                  yerr=[getfield(r, se_field) for r in sub], label="L = $L (GHZ)",
                  marker=markers[mod1(idx, length(markers))], markersize=5, linewidth=2)
        end
        hline!(p, [2 / 3], linestyle=:dash, color=:black, alpha=0.6, label="2/3")
        if !isempty(allup)
            for L in Ls
                pts = sort([(q, b) for ((l, q), b) in allup if l == L])
                isempty(pts) && continue
                plot!(p, first.(pts), last.(pts); linestyle=:dot, label="all-up, L = $L", linewidth=1.5)
            end
        end
        savefig(p, "dual_ea_right_boundary_ghz_$(name)_vs_q.pdf")
        savefig(p, "dual_ea_right_boundary_ghz_$(name)_vs_q.png")
    end

    println("\n✓ Wrote dual_ea_right_boundary_ghz_summary.csv and dual_ea_right_boundary_ghz_*_vs_q.{pdf,png}")
end

main()
