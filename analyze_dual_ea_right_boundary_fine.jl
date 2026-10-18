#!/usr/bin/env julia

"""
Analyze Right Boundary dual EA Binder FINE-grid Scan Results (q = 0.25:0.01:0.35)

Usage:
    julia analyze_dual_ea_right_boundary_fine.jl [results_dir] [min_L]

Loads every `dual_ea_right_boundary_fine_L*_...json` under <results_dir>/L*/ (default: the
newest `dual_ea_right_boundary_fine_results_*`), de-duplicates jobs with the same
(L, q, point_seed), pools the per-trajectory M2/M4/B arrays of all jobs at each (L, q) and
reports BOTH Binder conventions separately:

    B_mean_of_trials = mean over trajectories of B_m,      SE = std(B_m)/sqrt(n)
    B_ratio_of_means = 1 - mean(M4_m)/(3 mean(M2_m)^2),    SE = bootstrap over trajectories

`min_L` (default 16) drops smaller sizes at load time.

Because the q grid is fine (step 0.01) the pairwise crossings q*(L1, L2) of B_mean_of_trials
are located by linear interpolation between grid points, with a whole-trajectory bootstrap
(16-84%). A pair with no sign change on the grid is reported as such; nothing is forced.

Outputs: dual_ea_right_boundary_fine_summary.csv, dual_ea_right_boundary_fine_crossings.csv,
         dual_ea_right_boundary_fine_*_vs_q.{pdf,png}, dual_ea_right_boundary_fine_crossings.png
"""

using JSON
using Statistics
using Printf
using Random
using DelimitedFiles
using Plots

const NBOOT = 400

binder(M2, M4) = (isfinite(M2) && M2 > 0) ? 1.0 - M4 / (3.0 * M2^2) : NaN

function load_dual_ea_fine(results_dir; min_L::Int=16)
    pooled = Dict{Tuple{Int,Float64},Dict{String,Any}}()
    seen = Set{Tuple{Int,Float64,Int}}()
    n_dup = 0

    for L_dir in filter(isdir, readdir(results_dir, join=true))
        for file in filter(f -> startswith(basename(f), "dual_ea_right_boundary_fine_") &&
                                endswith(f, ".json") && !occursin("FAILED", f),
                           readdir(L_dir, join=true))
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

    n_dup > 0 && println("  skipped $n_dup duplicate jobs (same L, q, point_seed)")
    return pooled
end

function summarize(pooled; nboot=NBOOT)
    rows = []
    boots = Dict{Tuple{Int,Float64},Vector{Float64}}()   # bootstrap draws of B_mean_of_trials
    for ((L, q), g) in sort(collect(pooled); by=first)
        length(g["settings"]) > 1 &&
            @warn "Mixed (T_max, maxdim, cutoff) settings pooled at this point" L q settings=g["settings"]

        M2, M4, B = g["M2"], g["M4"], g["B"]
        n = length(B)
        rng = MersenneTwister(1000L + round(Int, 1000q))
        boot_ratio = Float64[]
        boot_mean = Float64[]
        for _ in 1:nboot
            idx = rand(rng, 1:n, n)
            push!(boot_ratio, binder(mean(M2[idx]), mean(M4[idx])))
            push!(boot_mean, mean(B[idx]))
        end
        boots[(L, q)] = boot_mean

        push!(rows, (
            L=L, q=q, n_traj=n, n_files=g["n_files"],
            M2_mean=mean(M2), M4_mean=mean(M4),
            B_mean_of_trials=mean(B), B_mean_of_trials_se=std(B) / sqrt(n),
            B_ratio_of_means=binder(mean(M2), mean(M4)),
            B_ratio_of_means_se=std(boot_ratio),
        ))
    end
    return rows, boots
end

# First sign change of d(q) on the grid, linearly interpolated; NaN if none.
function first_crossing(qs, d)
    for i in 1:(length(qs) - 1)
        if d[i] == 0.0
            return qs[i]
        elseif d[i] * d[i + 1] < 0
            return qs[i] - d[i] * (qs[i + 1] - qs[i]) / (d[i + 1] - d[i])
        end
    end
    return NaN
end

function crossing_table(rows, boots)
    Ls = sort(unique(r.L for r in rows))
    byL = Dict(L => sort(filter(r -> r.L == L, rows); by=r -> r.q) for L in Ls)
    out = []
    for (i, L1) in enumerate(Ls), L2 in Ls[(i + 1):end]
        common = sort(collect(intersect(Set(r.q for r in byL[L1]), Set(r.q for r in byL[L2]))))
        length(common) < 2 && continue
        y(L, q) = only(r.B_mean_of_trials for r in byL[L] if r.q == q)
        qstar = first_crossing(common, [y(L1, q) - y(L2, q) for q in common])
        nb = length(boots[(L1, common[1])])
        draws = [first_crossing(common, [boots[(L1, q)][k] - boots[(L2, q)][k] for q in common]) for k in 1:nb]
        ok = filter(isfinite, draws)
        push!(out, (L1=L1, L2=L2, adjacent=(findfirst(==(L2), Ls) == findfirst(==(L1), Ls) + 1),
                    qstar=qstar, lo=isempty(ok) ? NaN : quantile(ok, 0.16),
                    hi=isempty(ok) ? NaN : quantile(ok, 0.84), frac=length(ok) / nb,
                    Lbar=sqrt(L1 * L2)))
    end
    return out
end

function main()
    dirs = filter(d -> startswith(d, "dual_ea_right_boundary_fine_results_") && isdir(d), readdir())
    results_dir = !isempty(ARGS) ? ARGS[1] : (isempty(dirs) ? "" : last(sort(dirs)))
    min_L = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 16
    if isempty(results_dir) || !isdir(results_dir)
        println("ERROR: No results directory found (pass one, or run collect_dual_ea_right_boundary_fine_results.sh).")
        return
    end
    println("Loading fine-grid dual-EA results from: $results_dir (sizes L >= $min_L)")

    pooled = load_dual_ea_fine(results_dir; min_L=min_L)
    if isempty(pooled)
        println("ERROR: No dual_ea_right_boundary_fine_*.json files found under $results_dir/L*/")
        return
    end

    rows, boots = summarize(pooled)

    header = ["L" "q" "n_traj" "n_files" "M2_mean" "M4_mean" "B_mean_of_trials" "B_mean_of_trials_se" "B_ratio_of_means" "B_ratio_of_means_se"]
    table = Any[header]
    println()
    @printf("%4s %6s %7s %10s %10s %18s %18s\n", "L", "q", "n_traj", "M2_mean", "M4_mean",
            "B_mean_of_trials", "B_ratio_of_means")
    for r in rows
        push!(table, [r.L r.q r.n_traj r.n_files r.M2_mean r.M4_mean r.B_mean_of_trials r.B_mean_of_trials_se r.B_ratio_of_means r.B_ratio_of_means_se])
        @printf("%4d %6.3f %7d %10.5f %10.5f %10.5f ± %.5f %10.5f ± %.5f\n", r.L, r.q, r.n_traj,
                r.M2_mean, r.M4_mean, r.B_mean_of_trials, r.B_mean_of_trials_se,
                r.B_ratio_of_means, r.B_ratio_of_means_se)
    end
    writedlm("dual_ea_right_boundary_fine_summary.csv", vcat(table...), ',')

    # ---- crossings -------------------------------------------------------------
    cross = crossing_table(rows, boots)
    println("\nPairwise crossings of B_mean_of_trials (linear interpolation on the fine grid):")
    @printf("%9s %9s %20s %11s\n", "pair", "q*", "16-84% (bootstrap)", "boot found")
    ctable = Any[["L1" "L2" "adjacent" "qstar" "q16" "q84" "boot_found_frac"]]
    for c in cross
        if isfinite(c.qstar)
            @printf("%4d-%-4d %9.4f   [%.4f, %.4f] %9.0f%%  %s\n", c.L1, c.L2, c.qstar, c.lo, c.hi, 100 * c.frac, c.adjacent ? "adj" : "")
        else
            @printf("%4d-%-4d   no crossing on the grid (%.0f%% of bootstrap draws have one)\n", c.L1, c.L2, 100 * c.frac)
        end
        push!(ctable, [c.L1 c.L2 c.adjacent c.qstar c.lo c.hi c.frac])
    end
    writedlm("dual_ea_right_boundary_fine_crossings.csv", vcat(ctable...), ',')

    Ls = sort(unique(r.L for r in rows))
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5, :hexagon]

    for (name, field, se_field) in (
        ("B_mean_of_trials", :B_mean_of_trials, :B_mean_of_trials_se),
        ("B_ratio_of_means", :B_ratio_of_means, :B_ratio_of_means_se),
    )
        p = plot(xlabel="q = P_x = P_zz", ylabel="dual EA Binder ($name)",
                 title="Right boundary (fine grid): λ_x = 0.7, λ_zz = 0", legend=:outertopright,
                 grid=true, size=(850, 600), dpi=300)
        for (idx, L) in enumerate(Ls)
            sub = filter(r -> r.L == L, rows)
            plot!(p, [r.q for r in sub], [getfield(r, field) for r in sub];
                  yerr=[getfield(r, se_field) for r in sub], label="L = $L",
                  marker=markers[mod1(idx, length(markers))], markersize=5, linewidth=2)
        end
        savefig(p, "dual_ea_right_boundary_fine_$(name)_vs_q.pdf")
        savefig(p, "dual_ea_right_boundary_fine_$(name)_vs_q.png")
    end

    found = filter(c -> isfinite(c.qstar), cross)
    if !isempty(found)
        pc = plot(xlabel="geometric mean of the two sizes", ylabel="crossing q*",
                  title="pairwise crossings (bootstrap 16-84%)", legend=:bottomright,
                  grid=true, size=(700, 520), dpi=300)
        for adj in (true, false)
            sel = filter(c -> c.adjacent == adj, found)
            isempty(sel) && continue
            scatter!(pc, [c.Lbar for c in sel], [c.qstar for c in sel];
                     yerr=([c.qstar - c.lo for c in sel], [c.hi - c.qstar for c in sel]),
                     label=adj ? "adjacent sizes" : "non-adjacent", marker=adj ? :circle : :x)
        end
        savefig(pc, "dual_ea_right_boundary_fine_crossings.png")
    end

    println("\n✓ Wrote dual_ea_right_boundary_fine_summary.csv, dual_ea_right_boundary_fine_crossings.csv, " *
            "dual_ea_right_boundary_fine_*_vs_q.{pdf,png} and dual_ea_right_boundary_fine_crossings.png")
end

main()
