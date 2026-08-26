#!/usr/bin/env julia

"""
Aggregate and plot the reference-entropy diagnostic (see
reference_entropy_core.jl for the derivation and reference_entropy_run.jl
for the raw CSV schema).

Independent trajectories are the statistical unit; times within a
trajectory are correlated and are NEVER treated as independent samples.
All bootstrap resampling below resamples TRIAL IDs (preserving each
trial's full time series), never individual (trial, time) rows.
"""

using Statistics
using Printf
using Plots

# ============================================================
# Raw CSV loading (manual parser matching reference_entropy_run.jl's
# manually-written CSV).
# ============================================================
struct RawRow
    L::Int
    q::Float64
    lambda_x::Float64
    lambda_zz::Float64
    trial::Int
    seed::Int
    time::Int
    T_max::Int
    entropy::Float64
    trace_error::Float64
    hermiticity_error::Float64
    minimum_eigenvalue::Float64
    max_linkdim::Int
    record_checksum::String
end

function parse_float_field(s::AbstractString)
    s == "NaN" && return NaN
    s == "Inf" && return Inf
    s == "-Inf" && return -Inf
    return parse(Float64, s)
end

function load_raw_rows(files)
    rows = RawRow[]
    for file in files
        lines = readlines(file)
        isempty(lines) && continue
        header = split(lines[1], ",")
        col = Dict(name => idx for (idx, name) in enumerate(header))
        for line in lines[2:end]
            isempty(line) && continue
            f = split(line, ",")
            push!(rows, RawRow(
                parse(Int, f[col["L"]]), parse_float_field(f[col["q"]]),
                parse_float_field(f[col["lambda_x"]]), parse_float_field(f[col["lambda_zz"]]),
                parse(Int, f[col["trial"]]), parse(Int, f[col["seed"]]),
                parse(Int, f[col["time"]]), parse(Int, f[col["T_max"]]),
                parse_float_field(f[col["reference_entropy"]]),
                parse_float_field(f[col["trace_error"]]),
                parse_float_field(f[col["hermiticity_error"]]),
                parse_float_field(f[col["minimum_eigenvalue"]]),
                parse(Int, f[col["max_linkdim"]]),
                f[col["record_checksum"]],
            ))
        end
    end
    return rows
end

function bootstrap_ci(values::Vector{Float64}; nboot::Int=2000)
    n = length(values)
    n < 2 && return (lo=NaN, hi=NaN)
    samples = Vector{Float64}(undef, nboot)
    for b in 1:nboot
        idx = rand(1:n, n)
        samples[b] = mean(values[idx])
    end
    sort!(samples)
    lo = samples[max(1, round(Int, 0.025 * nboot))]
    hi = samples[min(nboot, round(Int, 0.975 * nboot))]
    return (lo=lo, hi=hi)
end

# ============================================================
# Aggregate per (L, q, time).
# ============================================================
function compute_statistics(rows::Vector{RawRow})
    groups = Dict{Tuple{Int,Float64,Int},Vector{RawRow}}()
    for row in rows
        push!(get!(groups, (row.L, row.q, row.time), RawRow[]), row)
    end

    stats = NamedTuple[]
    for ((L, q, t), grouprows) in groups
        entropies = [r.entropy for r in grouprows]
        n = length(entropies)
        ci = bootstrap_ci(entropies)
        push!(stats, (
            L=L, q=q, time=t, ntrials=n,
            mean_reference_entropy=mean(entropies),
            standard_error=n > 1 ? std(entropies) / sqrt(n) : NaN,
            bootstrap_ci_low=ci.lo, bootstrap_ci_high=ci.hi,
            mean_trace_error=mean(r.trace_error for r in grouprows),
            maximum_trace_error=maximum(r.trace_error for r in grouprows),
            minimum_reference_eigenvalue=minimum(r.minimum_eigenvalue for r in grouprows),
            maximum_bond_dimension=maximum(r.max_linkdim for r in grouprows),
        ))
    end
    sort!(stats, by=row -> (row.L, row.q, row.time))
    return stats
end

# ============================================================
# Plateau / steady-state convergence flag: compare the mean entropy at
# the last saved time to the mean entropy `n_back` snapshots earlier at
# the SAME (L,q); converged if the difference is smaller than the
# combined standard error (a documented, simple plateau criterion).
# ============================================================
function plateau_value(stats, L::Int, q::Float64; n_back::Int=2, tol_factor::Float64=2.0)
    data = sort(filter(row -> row.L == L && row.q == q, stats), by=row -> row.time)
    isempty(data) && return (value=NaN, se=NaN, converged=false, final_time=0)

    final = data[end]
    if length(data) <= n_back
        return (value=final.mean_reference_entropy, se=final.standard_error, converged=false, final_time=final.time)
    end

    earlier = data[end - n_back]
    diff = abs(final.mean_reference_entropy - earlier.mean_reference_entropy)
    combined_se = sqrt(final.standard_error^2 + earlier.standard_error^2)
    converged = isnan(combined_se) ? false : diff <= tol_factor * combined_se

    return (value=final.mean_reference_entropy, se=final.standard_error, converged=converged, final_time=final.time)
end

function compute_plateau_summary(stats)
    keys_Lq = unique((row.L, row.q) for row in stats)
    summary = NamedTuple[]
    for (L, q) in keys_Lq
        p = plateau_value(stats, L, q)
        push!(summary, (L=L, q=q, S_R=p.value, standard_error=p.se, converged=p.converged, final_time=p.final_time))
    end
    sort!(summary, by=row -> (row.L, row.q))
    return summary
end

# ============================================================
# Pairwise crossing estimates between adjacent L (linear interpolation
# of the plateau S_R(q) curves) with trajectory-level bootstrap CI.
# ============================================================
function find_crossing(qs::Vector{Float64}, y1::Vector{Float64}, y2::Vector{Float64})
    diff = y1 .- y2
    for k in 1:(length(qs) - 1)
        if diff[k] == 0
            return qs[k]
        end
        if sign(diff[k]) != sign(diff[k + 1])
            # linear interpolation for the zero crossing
            frac = diff[k] / (diff[k] - diff[k + 1])
            return qs[k] + frac * (qs[k + 1] - qs[k])
        end
    end
    return NaN
end

function compute_crossings(plateau_summary, raw_rows::Vector{RawRow}; nboot::Int=500)
    Ls = sort(unique(row.L for row in plateau_summary))
    crossings = NamedTuple[]

    for k in 1:(length(Ls) - 1)
        L1, L2 = Ls[k], Ls[k + 1]
        d1 = sort(filter(row -> row.L == L1, plateau_summary), by=row -> row.q)
        d2 = sort(filter(row -> row.L == L2, plateau_summary), by=row -> row.q)
        common_qs = sort(collect(intersect(Set(row.q for row in d1), Set(row.q for row in d2))))
        length(common_qs) < 2 && continue

        y1 = [only(row.S_R for row in d1 if row.q == q) for q in common_qs]
        y2 = [only(row.S_R for row in d2 if row.q == q) for q in common_qs]
        q_c = find_crossing(common_qs, y1, y2)

        # Trajectory-level bootstrap: resample trials WITHIN each (L,q) group
        # of the final-time raw entropies and recompute the crossing.
        boot_crossings = Float64[]
        final_rows_L1 = Dict(q => [r.entropy for r in raw_rows if r.L == L1 && r.q == q && r.time == d1[findfirst(row -> row.q == q, d1)].final_time] for q in common_qs)
        final_rows_L2 = Dict(q => [r.entropy for r in raw_rows if r.L == L2 && r.q == q && r.time == d2[findfirst(row -> row.q == q, d2)].final_time] for q in common_qs)
        for _ in 1:nboot
            yb1 = Float64[mean(v[rand(1:length(v), length(v))]) for v in (final_rows_L1[q] for q in common_qs)]
            yb2 = Float64[mean(v[rand(1:length(v), length(v))]) for v in (final_rows_L2[q] for q in common_qs)]
            qc_b = find_crossing(common_qs, yb1, yb2)
            isfinite(qc_b) && push!(boot_crossings, qc_b)
        end

        ci_lo = isempty(boot_crossings) ? NaN : simple_quantile(sort(boot_crossings), 0.025)
        ci_hi = isempty(boot_crossings) ? NaN : simple_quantile(sort(boot_crossings), 0.975)
        inv_L_mid = 2.0 / (L1 + L2)

        push!(crossings, (L1=L1, L2=L2, q_crossing=q_c, ci_low=ci_lo, ci_high=ci_hi, inv_L_mid=inv_L_mid))
    end

    return crossings
end

function simple_quantile(sorted_vals::Vector{Float64}, p::Float64)
    isempty(sorted_vals) && return NaN
    idx = clamp(round(Int, p * length(sorted_vals)), 1, length(sorted_vals))
    return sorted_vals[idx]
end

# ============================================================
# Finite-size scaling collapse: S_R(L,q) ~ f[(q - q_c) L^(1/nu)].
# Grid search over (q_c, nu) minimizing a binned cross-L variance
# collapse-quality metric; bootstrap over trials for CIs. This is a
# simple, transparent collapse metric -- NOT a full nonlinear
# regression framework.
# ============================================================
function collapse_cost(plateau_summary, q_c::Float64, nu::Float64; nbins::Int=15)
    xs = Float64[]
    ys = Float64[]
    for row in plateau_summary
        push!(xs, (row.q - q_c) * row.L^(1.0 / nu))
        push!(ys, row.S_R)
    end
    isempty(xs) && return Inf
    lo, hi = extrema(xs)
    (hi - lo) < 1e-12 && return Inf

    bin_of(x) = clamp(1 + floor(Int, (x - lo) / (hi - lo) * nbins), 1, nbins)
    bins = Dict{Int,Vector{Float64}}()
    for (x, y) in zip(xs, ys)
        push!(get!(bins, bin_of(x), Float64[]), y)
    end

    total = 0.0
    n_used = 0
    for (_, vals) in bins
        length(vals) < 2 && continue
        total += var(vals)
        n_used += 1
    end
    return n_used == 0 ? Inf : total / n_used
end

function fit_scaling_collapse(plateau_summary; q_c_range=0.0:0.01:0.5, nu_range=0.5:0.1:3.0)
    best = (q_c=NaN, nu=NaN, cost=Inf)
    for q_c in q_c_range, nu in nu_range
        cost = collapse_cost(plateau_summary, q_c, nu)
        if cost < best.cost
            best = (q_c=q_c, nu=nu, cost=cost)
        end
    end
    return best
end

function bootstrap_scaling_collapse(plateau_summary, raw_rows::Vector{RawRow}; nboot::Int=100)
    q_cs = Float64[]
    nus = Float64[]
    keys_Lq = [(row.L, row.q) for row in plateau_summary]
    final_time_of = Dict((row.L, row.q) => row.final_time for row in plateau_summary)

    for _ in 1:nboot
        resampled = NamedTuple[]
        for (L, q) in keys_Lq
            t = final_time_of[(L, q)]
            vals = [r.entropy for r in raw_rows if r.L == L && r.q == q && r.time == t]
            isempty(vals) && continue
            resampled_val = mean(vals[rand(1:length(vals), length(vals))])
            push!(resampled, (L=L, q=q, S_R=resampled_val))
        end
        fit = fit_scaling_collapse(resampled)
        if isfinite(fit.cost)
            push!(q_cs, fit.q_c)
            push!(nus, fit.nu)
        end
    end

    isempty(q_cs) && return (q_c_lo=NaN, q_c_hi=NaN, nu_lo=NaN, nu_hi=NaN)
    sort!(q_cs); sort!(nus)
    return (
        q_c_lo=simple_quantile(q_cs, 0.025), q_c_hi=simple_quantile(q_cs, 0.975),
        nu_lo=simple_quantile(nus, 0.025), nu_hi=simple_quantile(nus, 0.975),
    )
end

# ============================================================
# CSV output.
# ============================================================
function csv_field(x::AbstractFloat)
    isnan(x) && return "NaN"
    isinf(x) && return x > 0 ? "Inf" : "-Inf"
    return string(x)
end
csv_field(x) = string(x)

function save_csv(rows, header::Vector{String}, path::String)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(csv_field.(getfield.(Ref(row), Symbol.(header))), ","))
        end
    end
end

# ============================================================
# Plots.
# ============================================================
function plot_SR_vs_q(plateau_summary)
    p = plot(xlabel="q", ylabel="reference entropy S_R (converged value)",
              title="Reference entropy vs q", legend=:outertopright,
              grid=true, size=(800, 600), dpi=200)
    colors = palette(:tab10)
    for (idx, L) in enumerate(sort(unique(row.L for row in plateau_summary)))
        data = sort(filter(row -> row.L == L, plateau_summary), by=row -> row.q)
        isempty(data) && continue
        converged_mask = [row.converged for row in data]
        plot!(p, [row.q for row in data], [row.S_R for row in data];
              yerr=[row.standard_error for row in data],
              label="L = $L", color=colors[mod1(idx, length(colors))],
              marker=:circle, markersize=6, linewidth=2)
        unconverged = filter(row -> !row.converged, data)
        if !isempty(unconverged)
            scatter!(p, [row.q for row in unconverged], [row.S_R for row in unconverged];
                      marker=:xcross, markersize=10, color=:red, label=(idx == 1 ? "unconverged" : ""))
        end
    end
    savefig(p, "reference_entropy_SR_vs_q.pdf")
    savefig(p, "reference_entropy_SR_vs_q.png")
end

function plot_SR_vs_t_over_L(stats)
    Ls = sort(unique(row.L for row in stats))
    subplots = Plots.Plot[]
    colors = palette(:tab10)
    for (li, L) in enumerate(Ls)
        p = plot(xlabel="t/L", ylabel="S_R(t)", title="L = $L",
                  legend=(li == length(Ls) ? :outertopright : false))
        panel = filter(row -> row.L == L, stats)
        for (qi, q) in enumerate(sort(unique(row.q for row in panel)))
            data = sort(filter(row -> row.q == q, panel), by=row -> row.time)
            isempty(data) && continue
            plot!(p, [row.time / L for row in data], [row.mean_reference_entropy for row in data];
                  yerr=[row.standard_error for row in data],
                  label=@sprintf("q=%.2f", q), color=colors[mod1(qi, length(colors))],
                  marker=:circle, markersize=3)
        end
        push!(subplots, p)
    end
    combined = plot(subplots...; layout=(1, length(Ls)), size=(1800, 500),
                     plot_title="S_R(t) vs t/L", dpi=200)
    savefig(combined, "reference_entropy_SR_vs_t_over_L.pdf")
    savefig(combined, "reference_entropy_SR_vs_t_over_L.png")
end

function plot_crossings_vs_q(plateau_summary, crossings)
    p = plot(xlabel="q", ylabel="S_R", title="Pairwise curve crossings",
              legend=:outertopright, grid=true, size=(800, 600), dpi=200)
    colors = palette(:tab10)
    for (idx, L) in enumerate(sort(unique(row.L for row in plateau_summary)))
        data = sort(filter(row -> row.L == L, plateau_summary), by=row -> row.q)
        plot!(p, [row.q for row in data], [row.S_R for row in data];
              label="L = $L", color=colors[mod1(idx, length(colors))], marker=:circle)
    end
    for c in crossings
        isfinite(c.q_crossing) && vline!(p, [c.q_crossing]; linestyle=:dash, color=:gray, label="")
    end
    savefig(p, "reference_entropy_pairwise_crossings.pdf")
    savefig(p, "reference_entropy_pairwise_crossings.png")
end

function plot_crossing_vs_inv_L(crossings)
    isempty(crossings) && return
    p = plot(xlabel="1/L (pair midpoint)", ylabel="crossing q_c estimate",
              title="Crossing location vs inverse system size", legend=false,
              grid=true, size=(800, 600), dpi=200)
    xs = [c.inv_L_mid for c in crossings]
    ys = [c.q_crossing for c in crossings]
    yerr = [(c.q_crossing - c.ci_low, c.ci_high - c.q_crossing) for c in crossings]
    plot!(p, xs, ys; ribbon=(first.(yerr), last.(yerr)), marker=:circle, markersize=6)
    savefig(p, "reference_entropy_crossing_vs_inv_L.pdf")
    savefig(p, "reference_entropy_crossing_vs_inv_L.png")
end

function plot_scaling_collapse(plateau_summary, fit)
    isnan(fit.q_c) && return
    p = plot(xlabel="(q - q_c) L^(1/nu)", ylabel="S_R",
              title=@sprintf("Scaling collapse: q_c=%.3f, nu=%.2f", fit.q_c, fit.nu),
              legend=:outertopright, grid=true, size=(800, 600), dpi=200)
    colors = palette(:tab10)
    for (idx, L) in enumerate(sort(unique(row.L for row in plateau_summary)))
        data = filter(row -> row.L == L, plateau_summary)
        xs = [(row.q - fit.q_c) * L^(1.0 / fit.nu) for row in data]
        ys = [row.S_R for row in data]
        order = sortperm(xs)
        plot!(p, xs[order], ys[order]; label="L = $L", color=colors[mod1(idx, length(colors))], marker=:circle)
    end
    savefig(p, "reference_entropy_scaling_collapse.pdf")
    savefig(p, "reference_entropy_scaling_collapse.png")
end

function default_result_files()
    result_dirs = filter(d -> startswith(d, "reference_entropy_results_") && isdir(d), readdir())
    isempty(result_dirs) && return String[]
    results_dir = last(sort(result_dirs))
    return String[file for file in readdir(results_dir, join=true) if endswith(file, ".csv")]
end

function print_interpretation_guide(plateau_summary, crossings, fit)
    has_crossing = any(isfinite(c.q_crossing) for c in crossings)
    println("""

    ================= Interpretation guide =================
    S_R ~ 1 means the measurement record has learned little about the
    encoded logical bit (reference stays entangled with the system).
    S_R ~ 0 means the record has learned the logical information
    (conditional reference state becomes nearly pure).

    Pairwise crossings found: $(count(c -> isfinite(c.q_crossing), crossings))/$(length(crossings))
    """)
    if fit !== nothing && isfinite(fit.q_c)
        @printf("    Scaling-collapse fit: q_c ~ %.3f, nu ~ %.2f (cost=%.4g)\n", fit.q_c, fit.nu, fit.cost)
    end
    println("""

    A credible transition requires ALL of: crossings drifting toward a
    common value as L increases; adequate steady-state convergence
    (see the `converged` flag); bootstrap uncertainty smaller than the
    curve separation; a stable scaling collapse; and results converged
    with respect to MPS truncation (see reference_entropy_validate.jl's
    convergence test). This script reports the raw crossing/collapse
    numbers computed above; it does NOT itself certify a transition.
    If curves are ordered without crossing, reference entropy does not
    detect a transition for the studied parameter slice and sizes --
    report that explicitly rather than forcing an interpretation.

    Reference entropy and the Rényi-2 Binder are INDEPENDENT
    diagnostics. Do not force them to agree: reference entropy probes
    what the measurement record has learned about one encoded logical
    qubit, while the Binder cumulant probes replica-overlap order-
    parameter fluctuations. A crossing (or lack of one) in S_R should be
    reported and compared to the Binder crossing (or lack of one)
    ONLY after both have independently satisfied their own convergence
    and uncertainty requirements.
    ==========================================================
    """)
end

function main()
    files = isempty(ARGS) ? default_result_files() : collect(ARGS)
    if isempty(files)
        println("ERROR: No CSV files found. Pass files explicitly or run " *
                "reference_entropy_collect.sh first.")
        return
    end

    println("Loading $(length(files)) file(s)...")
    raw_rows = load_raw_rows(files)
    println("Loaded $(length(raw_rows)) raw (trial, time) rows")

    stats = compute_statistics(raw_rows)
    save_csv(stats, [
        "L", "q", "time", "ntrials", "mean_reference_entropy", "standard_error",
        "bootstrap_ci_low", "bootstrap_ci_high", "mean_trace_error",
        "maximum_trace_error", "minimum_reference_eigenvalue", "maximum_bond_dimension",
    ], "reference_entropy_aggregated_statistics.csv")
    println("Saved reference_entropy_aggregated_statistics.csv ($(length(stats)) rows)")

    plateau_summary = compute_plateau_summary(stats)
    save_csv(plateau_summary, ["L", "q", "S_R", "standard_error", "converged", "final_time"],
              "reference_entropy_plateau_summary.csv")

    println("\nSummary (L, q, plateau S_R, SE, converged):")
    for row in plateau_summary
        @printf("  L=%-3d q=%.2f  S_R=%.4f ± %-8.4f converged=%s\n",
                row.L, row.q, row.S_R, row.standard_error, row.converged)
    end

    crossings = compute_crossings(plateau_summary, raw_rows)
    save_csv(crossings, ["L1", "L2", "q_crossing", "ci_low", "ci_high", "inv_L_mid"],
              "reference_entropy_crossings.csv")

    fit = length(unique(row.L for row in plateau_summary)) >= 2 ?
        fit_scaling_collapse(plateau_summary) : (q_c=NaN, nu=NaN, cost=Inf)
    if isfinite(fit.q_c)
        ci = bootstrap_scaling_collapse(plateau_summary, raw_rows)
        open("reference_entropy_scaling_fit.csv", "w") do io
            println(io, "q_c,nu,cost,q_c_ci_low,q_c_ci_high,nu_ci_low,nu_ci_high")
            println(io, join(csv_field.([fit.q_c, fit.nu, fit.cost, ci.q_c_lo, ci.q_c_hi, ci.nu_lo, ci.nu_hi]), ","))
        end
    end

    println("\nGenerating plots...")
    plot_SR_vs_q(plateau_summary)
    plot_SR_vs_t_over_L(stats)
    plot_crossings_vs_q(plateau_summary, crossings)
    plot_crossing_vs_inv_L(crossings)
    plot_scaling_collapse(plateau_summary, fit)

    println("✓ Plots generated:")
    for name in (
        "reference_entropy_SR_vs_q", "reference_entropy_SR_vs_t_over_L",
        "reference_entropy_pairwise_crossings", "reference_entropy_crossing_vs_inv_L",
        "reference_entropy_scaling_collapse",
    )
        println("  - $name.pdf")
    end

    print_interpretation_guide(plateau_summary, crossings, fit)
end

main()
