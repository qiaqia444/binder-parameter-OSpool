#!/usr/bin/env julia

"""
Analyze Right Boundary Renyi-2 Binder measured AFTER THE LAST X DEPHASING

Usage:
    julia --project=. analyze_renyi2_right_boundary_xnoise.jl [xnoise_results_dir] [production_results_dir]

* xnoise_results_dir     : output of collect_renyi2_right_boundary_xnoise_results.sh
                           (default: newest renyi2_xnoise_right_boundary_results_*)
* production_results_dir : OPTIONAL output of collect_right_boundary_results.sh (production Renyi-2,
                           measured after the LAST ZZ DEPHASING; default: newest right_boundary_results_*
                           if present, pass "none" to skip the comparison)

Per (L, q): the primary estimate is the mean over jobs of the per-job trajectory-averaged Binder
`B` (jobs weighted by their number of valid trajectories); its standard error is the job-to-job
scatter / sqrt(n_jobs). The ratio of mean moments 1 - M4bar/(3 M2bar^2) is reported separately as a
diagnostic. Because both families use the SAME seeds (job k of one family evolves the same Born
trajectories as job k of the other), jobs are PAIRED by (L, q, sample); the paired difference
B(after last ZZ) - B(after last X) and its job-level standard error are reported. Missing partner
jobs are simply not paired. No crossing is assumed.

Outputs: renyi2_xnoise_right_boundary_summary.csv and renyi2_xnoise_right_boundary_*.{pdf,png}
"""

using JSON
using Statistics
using Printf
using DelimitedFiles
using Plots

const JOB_RE = r"_Px([0-9.]+)_s(\d+)\.json$"

# Returns Dict{(L,q,sample)} => (B, M2, M4, n_valid, n_invalid)
function load_jobs(results_dir; prefix)
    jobs = Dict{Tuple{Int,Float64,Int},NamedTuple}()
    for L_dir in filter(isdir, readdir(results_dir, join=true))
        for file in sort(readdir(L_dir, join=true))
            name = basename(file)
            (startswith(name, prefix) && endswith(name, ".json") && !occursin("FAILED", name)) || continue
            m = match(JOB_RE, name)
            m === nothing && continue
            try
                parsed = JSON.parsefile(file)
                for e in (parsed isa AbstractVector ? parsed : [parsed])
                    e["B"] === nothing && continue
                    n = Int(get(e, "n_valid", get(e, "ntrials", 1)))
                    jobs[(Int(e["L"]), round(Float64(e["P_x"]); digits=6), parse(Int, m.captures[2]))] = (
                        B=Float64(e["B"]), M2=Float64(e["M2_bar"]), M4=Float64(e["M4_bar"]),
                        n=n, n_invalid=Int(get(e, "n_invalid", 0)),
                    )
                end
            catch err
                @warn "Could not read JSON file" file exception=(err, catch_backtrace())
            end
        end
    end
    return jobs
end

binder(M2, M4) = (isfinite(M2) && M2 > 0) ? 1.0 - M4 / (3.0 * M2^2) : NaN

function summarize(jobs)
    groups = Dict{Tuple{Int,Float64},Vector{NamedTuple}}()
    for ((L, q, _), j) in jobs
        push!(get!(groups, (L, q), NamedTuple[]), j)
    end
    S = Dict{Tuple{Int,Float64},NamedTuple}()
    for (k, js) in groups
        w = [j.n for j in js]
        W = sum(w)
        B = sum(w .* [j.B for j in js]) / W
        M2 = sum(w .* [j.M2 for j in js]) / W
        M4 = sum(w .* [j.M4 for j in js]) / W
        se = length(js) > 1 ? std([j.B for j in js]) / sqrt(length(js)) : NaN
        S[k] = (B=B, B_se=se, B_ratio=binder(M2, M4), M2=M2, M4=M4, njobs=length(js), ntraj=W,
                n_invalid=sum(j.n_invalid for j in js))
    end
    return S
end

function paired_difference(jx, jz)
    diffs = Dict{Tuple{Int,Float64},Vector{Float64}}()
    for (key, a) in jx
        haskey(jz, key) || continue
        push!(get!(diffs, (key[1], key[2]), Float64[]), jz[key].B - a.B)
    end
    return Dict(k => (d=mean(v), se=length(v) > 1 ? std(v) / sqrt(length(v)) : NaN, n=length(v)) for (k, v) in diffs)
end

function main()
    xdirs = filter(d -> startswith(d, "renyi2_xnoise_right_boundary_results_") && isdir(d), readdir())
    xdir = length(ARGS) >= 1 ? ARGS[1] : (isempty(xdirs) ? "" : last(sort(xdirs)))
    if isempty(xdir) || !isdir(xdir)
        println("ERROR: no X-noise results directory (pass one, or run collect_renyi2_right_boundary_xnoise_results.sh)")
        return
    end
    pdirs = filter(d -> startswith(d, "right_boundary_results_") && isdir(d), readdir())
    pdir = length(ARGS) >= 2 ? ARGS[2] : (isempty(pdirs) ? "none" : last(sort(pdirs)))
    println("X-noise results     : $xdir")
    println("production results  : $pdir")

    jx = load_jobs(xdir; prefix="renyi2_xnoise_right_boundary_")
    isempty(jx) && (println("ERROR: no renyi2_xnoise_right_boundary_*.json files under $xdir/L*/"); return)
    Sx = summarize(jx)

    have_prod = pdir != "none" && isdir(pdir)
    jz = have_prod ? load_jobs(pdir; prefix="right_boundary_") : Dict{Tuple{Int,Float64,Int},NamedTuple}()
    Sz = isempty(jz) ? Dict{Tuple{Int,Float64},NamedTuple}() : summarize(jz)
    PD = isempty(jz) ? Dict{Tuple{Int,Float64},NamedTuple}() : paired_difference(jx, jz)

    keys_sorted = sort(collect(keys(Sx)))
    @printf("\n%4s %6s %6s | %-19s %-9s | %-19s | %-19s\n", "L", "q", "jobs", "B after last X (se)", "B ratio", "B after last ZZ (se)", "paired dB = ZZ - X")
    table = Any["L" "q" "n_jobs" "n_traj" "B_after_X" "B_after_X_se" "B_ratio_after_X" "M2_after_X" "M4_after_X" "B_after_ZZ" "B_after_ZZ_se" "dB_paired" "dB_paired_se" "n_paired_jobs"]
    for (L, q) in keys_sorted
        x = Sx[(L, q)]
        z = get(Sz, (L, q), nothing)
        d = get(PD, (L, q), nothing)
        @printf("%4d %6.3f %6d | %8.5f ± %-8.5f %9.5f | %-19s | %s\n", L, q, x.njobs, x.B, x.B_se, x.B_ratio,
                z === nothing ? "-" : @sprintf("%8.5f ± %-8.5f", z.B, z.B_se),
                d === nothing ? "-" : @sprintf("%8.5f ± %-8.5f (%d)", d.d, d.se, d.n))
        table = vcat(table, [L q x.njobs x.ntraj x.B x.B_se x.B_ratio x.M2 x.M4 (z === nothing ? NaN : z.B) (z === nothing ? NaN : z.B_se) (d === nothing ? NaN : d.d) (d === nothing ? NaN : d.se) (d === nothing ? 0 : d.n)])
    end
    writedlm("renyi2_xnoise_right_boundary_summary.csv", table, ',')

    Ls = sort(unique(k[1] for k in keys_sorted))
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5, :hexagon]
    colors = [:blue, :red, :green, :purple, :orange, :brown, :black]
    rowsof(S, L) = sort([(q, S[(l, q)]) for (l, q) in keys(S) if l == L]; by=first)

    p1 = plot(xlabel="P_x = P_zz", ylabel="Binder Parameter (Rényi-2)", title="λ_x = 0.7, λ_zz = 0.0: after the last X dephasing",
              legend=:topright, grid=true, size=(900, 650), dpi=200)
    for (i, L) in enumerate(Ls)
        r = rowsof(Sx, L)
        plot!(p1, first.(r), [x.B for (_, x) in r]; yerr=[x.B_se for (_, x) in r], color=colors[mod1(i, 7)], marker=markers[mod1(i, 7)],
              markersize=6, linewidth=2, markerstrokecolor=:black, label="L = $L")
    end
    hline!(p1, [2 / 3]; linestyle=:dash, color=:black, alpha=0.7, label="B = 2/3")
    savefig(p1, "renyi2_xnoise_right_boundary_binder_vs_q.pdf"); savefig(p1, "renyi2_xnoise_right_boundary_binder_vs_q.png")

    if !isempty(Sz)
        p2 = plot(xlabel="P_x = P_zz", ylabel="Binder Parameter (Rényi-2)", title="observation time: after last X (solid) vs after last ZZ (dashed)",
                  legend=:topright, grid=true, size=(900, 650), dpi=200)
        for (i, L) in enumerate(Ls)
            r, rz = rowsof(Sx, L), rowsof(Sz, L)
            plot!(p2, first.(r), [x.B for (_, x) in r]; color=colors[mod1(i, 7)], marker=markers[mod1(i, 7)], linewidth=2, label="L = $L, after X")
            isempty(rz) || plot!(p2, first.(rz), [x.B for (_, x) in rz]; color=colors[mod1(i, 7)], markershape=:x, linestyle=:dash, linewidth=1.5, label="L = $L, after ZZ")
        end
        savefig(p2, "renyi2_xnoise_right_boundary_vs_after_zz.pdf"); savefig(p2, "renyi2_xnoise_right_boundary_vs_after_zz.png")

        p3 = plot(xlabel="P_x = P_zz", ylabel="B(after last ZZ) - B(after last X)", title="paired difference (same Born trajectories)", legend=:topleft, grid=true, size=(900, 650), dpi=200)
        for (i, L) in enumerate(Ls)
            r = sort([(q, d) for ((l, q), d) in PD if l == L]; by=first)
            isempty(r) || plot!(p3, first.(r), [d.d for (_, d) in r]; yerr=[d.se for (_, d) in r], color=colors[mod1(i, 7)], marker=markers[mod1(i, 7)], linewidth=2, label="L = $L")
        end
        hline!(p3, [0.0]; color=:black, linestyle=:dash, label="")
        savefig(p3, "renyi2_xnoise_right_boundary_paired_difference.pdf"); savefig(p3, "renyi2_xnoise_right_boundary_paired_difference.png")
    end

    println("\n✓ Wrote renyi2_xnoise_right_boundary_summary.csv and renyi2_xnoise_right_boundary_*.{pdf,png}")
end

main()
