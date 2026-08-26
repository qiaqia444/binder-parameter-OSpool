#!/usr/bin/env julia

"""
Aggregate and plot the right-boundary Renyi-2 disorder-string overlap
diagnostic (see defect_insertion_core.jl for the full derivation).

R_ij = zij/z0 is a SIGNED overlap, not a partition-function ratio. This
script never discards a trajectory, never clips R, and never converts
R to |R| while dropping its sign; it reports signed statistics
(mean_R, average_sign, positive_fraction, ...) alongside magnitude-only
statistics (mean_absR, string_log_magnitude = -log|R|).

Pseudoreplication: each separation r has TWO translated endpoint pairs
measured on the SAME trajectory (not independent samples). This script
averages the two pairs' R within each trajectory FIRST, and only then
computes standard errors / bootstrap CIs across independent trajectories.
"""

using Statistics
using Printf
using Plots

const ZERO_ATOL = 1e-12  # must match defect_insertion_core.jl's DEFECT_INSERTION_ZERO_ATOL

# ============================================================
# Raw CSV loading (manual parser; matches defect_insertion_run.jl's
# manually-written CSV, so no CSV.jl dependency is required).
# ============================================================
struct RawRow
    L::Int
    q::Float64
    trial::Int
    seed::Int
    i::Int
    j::Int
    r::Int
    z0::Float64
    zij::Float64
    R::Float64
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
            fields = split(line, ",")
            push!(rows, RawRow(
                parse(Int, fields[col["L"]]),
                parse_float_field(fields[col["q"]]),
                parse(Int, fields[col["trial"]]),
                parse(Int, fields[col["seed"]]),
                parse(Int, fields[col["i"]]),
                parse(Int, fields[col["j"]]),
                parse(Int, fields[col["r"]]),
                parse_float_field(fields[col["z0"]]),
                parse_float_field(fields[col["zij"]]),
                parse_float_field(fields[col["R"]]),
                fields[col["record_checksum"]],
            ))
        end
    end
    return rows
end

# ============================================================
# Step 1: average translated pairs WITHIN each trajectory (avoids
# pseudoreplication -- the two pairs at fixed (L,q,r,trial) are NOT
# independent samples).
# ============================================================
function average_translated_pairs(rows::Vector{RawRow})
    groups = Dict{Tuple{Int,Float64,Int,Int},Vector{Float64}}()
    for row in rows
        key = (row.L, row.q, row.r, row.trial)
        push!(get!(groups, key, Float64[]), row.R)
    end

    R_traj = NamedTuple[]
    for ((L, q, r, trial), Rs) in groups
        push!(R_traj, (L=L, q=q, r=r, trial=trial, R=mean(Rs), n_pairs=length(Rs)))
    end
    return R_traj
end

# ============================================================
# Step 2: aggregate the per-trajectory R over independent trials.
# ============================================================
function classify_sign(R::Float64)
    return abs(R) <= ZERO_ATOL ? 0.0 : sign(R)
end

function string_log_magnitude_of(R::Float64)
    return abs(R) <= ZERO_ATOL ? Inf : -log(abs(R))
end

function bootstrap_ci(values::Vector{Float64}; nboot::Int=2000, rng=nothing)
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

function compute_statistics(R_traj::Vector{NamedTuple})
    groups = Dict{Tuple{Int,Float64,Int},Vector{Float64}}()
    for row in R_traj
        key = (row.L, row.q, row.r)
        push!(get!(groups, key, Float64[]), row.R)
    end

    stats = NamedTuple[]
    for ((L, q, r), Rs) in groups
        n_total = length(Rs)
        signs = classify_sign.(Rs)
        n_positive = count(==(1.0), signs)
        n_negative = count(==(-1.0), signs)
        n_zero = count(==(0.0), signs)
        @assert n_positive + n_negative + n_zero == n_total  # bookkeeping identity

        absRs = abs.(Rs)
        slm = string_log_magnitude_of.(Rs)
        finite_slm = filter(isfinite, slm)

        mean_R = mean(Rs)
        se_R = n_total > 1 ? std(Rs) / sqrt(n_total) : NaN
        ci_R = bootstrap_ci(Rs)
        ci_absR = bootstrap_ci(absRs)
        ci_slm = bootstrap_ci(finite_slm)

        push!(stats, (
            L=L, q=q, r=r, n_total=n_total,
            n_positive=n_positive, n_negative=n_negative, n_zero=n_zero,
            mean_R=mean_R, standard_error_R=se_R,
            mean_R_ci_lo=ci_R.lo, mean_R_ci_hi=ci_R.hi,
            mean_absR=mean(absRs), median_absR=median(absRs),
            mean_absR_ci_lo=ci_absR.lo, mean_absR_ci_hi=ci_absR.hi,
            positive_fraction=n_positive / n_total,
            average_sign=mean(signs),
            median_string_log_magnitude=median(slm),
            mean_finite_string_log_magnitude=isempty(finite_slm) ? NaN : mean(finite_slm),
            n_finite_log_magnitude=length(finite_slm),
            n_excluded_log_magnitude=n_total - length(finite_slm),
            string_log_magnitude_ci_lo=ci_slm.lo, string_log_magnitude_ci_hi=ci_slm.hi,
        ))
    end

    sort!(stats, by=row -> (row.L, row.q, row.r))
    return stats
end

# ============================================================
# CSV output (aggregated statistics).
# ============================================================
const AGGREGATE_CSV_HEADER = [
    "L", "q", "r", "n_total", "n_positive", "n_negative", "n_zero",
    "mean_R", "standard_error_R", "mean_R_ci_lo", "mean_R_ci_hi",
    "mean_absR", "median_absR", "mean_absR_ci_lo", "mean_absR_ci_hi",
    "positive_fraction", "average_sign",
    "median_string_log_magnitude", "mean_finite_string_log_magnitude",
    "n_finite_log_magnitude", "n_excluded_log_magnitude",
    "string_log_magnitude_ci_lo", "string_log_magnitude_ci_hi",
]

function csv_field(x::AbstractFloat)
    isnan(x) && return "NaN"
    isinf(x) && return x > 0 ? "Inf" : "-Inf"
    return string(x)
end
csv_field(x) = string(x)

function save_aggregated_csv(stats, path)
    open(path, "w") do io
        println(io, join(AGGREGATE_CSV_HEADER, ","))
        for row in stats
            println(io, join(csv_field.(getfield.(Ref(row), Symbol.(AGGREGATE_CSV_HEADER))), ","))
        end
    end
end

# ============================================================
# Plots. Labels explicitly say "signed disorder-string overlap" /
# "string log-magnitude" -- never "free energy" or "partition function."
# ============================================================
function panel_by_L(stats, Ls)
    return [filter(row -> row.L == L, stats) for L in Ls]
end

function plot_metric_vs_r(stats, Ls; ykey::Symbol, ylabel::String, title::String, filename::String,
                           ci_lo_key=nothing, ci_hi_key=nothing)
    panels = panel_by_L(stats, Ls)
    subplots = Plots.Plot[]
    colors = palette(:tab10)

    for (idx, (Lval, panel)) in enumerate(zip(Ls, panels))
        p = plot(xlabel="r", ylabel=ylabel, title="L = $Lval", legend=(idx == length(Ls) ? :outertopright : false))
        for (cidx, q) in enumerate(sort(unique(panel.q for panel in panel)))
            data = sort(filter(row -> row.q == q, panel), by=row -> row.r)
            isempty(data) && continue
            y = getfield.(data, ykey)
            if ci_lo_key !== nothing
                lo = getfield.(data, ci_lo_key)
                hi = getfield.(data, ci_hi_key)
                yerr = [(y[k] - lo[k], hi[k] - y[k]) for k in eachindex(y)]
                plot!(p, [row.r for row in data], y; ribbon=(first.(yerr), last.(yerr)),
                      label=@sprintf("q=%.2f", q), color=colors[mod1(cidx, length(colors))],
                      marker=:circle, markersize=4)
            else
                plot!(p, [row.r for row in data], y; label=@sprintf("q=%.2f", q),
                      color=colors[mod1(cidx, length(colors))], marker=:circle, markersize=4)
            end
        end
        push!(subplots, p)
    end

    combined = plot(subplots...; layout=(1, length(Ls)), size=(1800, 500), plot_title=title, dpi=200)
    savefig(combined, filename * ".pdf")
    savefig(combined, filename * ".png")
end

function plot_half_system_vs_L(stats, Ls)
    half_rows = filter(row -> row.r == div(row.L, 2), stats)

    p_lin = plot(xlabel="L", ylabel="string log-magnitude (r = L/2)",
                 title="Half-system string log-magnitude vs L", legend=:outertopright,
                 grid=true, size=(800, 600), dpi=200)
    p_log = plot(xlabel="log L", ylabel="string log-magnitude (r = L/2)",
                 title="Half-system string log-magnitude vs log L", legend=:outertopright,
                 grid=true, size=(800, 600), dpi=200)
    colors = palette(:tab10)

    for (idx, q) in enumerate(sort(unique(row.q for row in half_rows)))
        data = sort(filter(row -> row.q == q, half_rows), by=row -> row.L)
        isempty(data) && continue
        y = getfield.(data, :mean_finite_string_log_magnitude)
        lo = getfield.(data, :string_log_magnitude_ci_lo)
        hi = getfield.(data, :string_log_magnitude_ci_hi)
        yerr = ([y[k] - lo[k] for k in eachindex(y)], [hi[k] - y[k] for k in eachindex(y)])
        color = colors[mod1(idx, length(colors))]
        plot!(p_lin, [row.L for row in data], y; ribbon=yerr,
              label=@sprintf("q=%.2f", q), color=color, marker=:circle, markersize=6)
        plot!(p_log, log.([row.L for row in data]), y; ribbon=yerr,
              label=@sprintf("q=%.2f", q), color=color, marker=:circle, markersize=6)
    end

    savefig(p_lin, "defect_insertion_half_system_string_log_magnitude_vs_L.pdf")
    savefig(p_lin, "defect_insertion_half_system_string_log_magnitude_vs_L.png")
    savefig(p_log, "defect_insertion_half_system_string_log_magnitude_vs_log_L.pdf")
    savefig(p_log, "defect_insertion_half_system_string_log_magnitude_vs_log_L.png")
end

function plot_string_log_magnitude_vs_log_r(stats, Ls)
    panels = panel_by_L(stats, Ls)
    subplots = Plots.Plot[]
    colors = palette(:tab10)

    for (idx, (Lval, panel)) in enumerate(zip(Ls, panels))
        p = plot(xlabel="log r", ylabel="string log-magnitude", title="L = $Lval",
                  legend=(idx == length(Ls) ? :outertopright : false))
        for (cidx, q) in enumerate(sort(unique(row.q for row in panel)))
            data = sort(filter(row -> row.q == q, panel), by=row -> row.r)
            isempty(data) && continue
            plot!(p, log.([row.r for row in data]), getfield.(data, :mean_finite_string_log_magnitude);
                  label=@sprintf("q=%.2f", q), color=colors[mod1(cidx, length(colors))],
                  marker=:circle, markersize=4)
        end
        push!(subplots, p)
    end

    combined = plot(subplots...; layout=(1, length(Ls)), size=(1800, 500),
                     plot_title="String log-magnitude vs log r", dpi=200)
    savefig(combined, "defect_insertion_string_log_magnitude_vs_log_r.pdf")
    savefig(combined, "defect_insertion_string_log_magnitude_vs_log_r.png")
end

function default_result_files()
    result_dirs = filter(d -> startswith(d, "defect_insertion_results_") && isdir(d), readdir())
    isempty(result_dirs) && return String[]
    results_dir = last(sort(result_dirs))
    return String[
        file
        for L in (16, 24, 32)
        for subdir in (joinpath(results_dir, "L$L"),)
        if isdir(subdir)
        for file in readdir(subdir, join=true)
        if endswith(file, ".csv")
    ]
end

function print_interpretation_guide()
    println("""

    ================= Interpretation guide =================
    R_ij(m) = zij(m)/z0(m) is a SIGNED Renyi-2 disorder-string overlap
    for one trajectory m. It is NOT automatically a partition-function
    ratio: showing D_ij mu_p D_ij = -mu_p on the segment does not by
    itself prove R_ij = Z_DW/Z_0 for two independently positive classical
    partition functions. That interpretation requires a separately
    derived statistical-mechanics mapping (modified bond weights or
    boundary conditions, both contractions independently proven
    positive), which has NOT been derived here.

    Comparing with the Rényi-2 Binder results requires checking, in
    order:
      1. Signed data: mean_R, positive_fraction, and average_sign across
         q. If R changes sign frequently or averages near zero even in a
         region where Binder shows strong order, the two diagnostics are
         NOT measuring the same thing in a naive sense.
      2. Distance (r) and size (L) scaling of string_log_magnitude:
           -log|R(r)| ~ sigma*r + C     (putative exponential decay)
           -log|R(r)| ~ 2*x*log(r) + C  (putative algebraic/critical decay)
           -log|R(r)| ~ C               (putative long-distance saturation)
         A transition requires a QUALITATIVE CHANGE in this r- or
         L-scaling across q -- not merely an increase in
         string_log_magnitude with q at fixed r or fixed L=r/2.
      3. Numerical-zero behavior: whether R is compatible with an
         analytically-zero overlap at a given q (as found here at
         q = 0.5; see defect_insertion_validate.jl's "q = 0.5
         investigation" test) before assigning any physical meaning to a
         large string_log_magnitude value.

    Do NOT claim agreement or disagreement with the Binder transition
    until all three of the above have been checked for the actual
    production data.
    ==========================================================
    """)
end

function main()
    files = isempty(ARGS) ? default_result_files() : collect(ARGS)
    if isempty(files)
        println("ERROR: No CSV files found. Pass files explicitly or run " *
                "defect_insertion_collect.sh first.")
        return
    end

    println("Loading $(length(files)) file(s)...")
    raw_rows = load_raw_rows(files)
    println("Loaded $(length(raw_rows)) raw pair-level rows")

    R_traj = average_translated_pairs(raw_rows)
    println("Averaged to $(length(R_traj)) (L,q,r,trial) trajectory rows " *
            "(translated pairs averaged within trajectory first)")

    stats = compute_statistics(R_traj)
    save_aggregated_csv(stats, "defect_insertion_aggregated_statistics.csv")
    println("Saved aggregated statistics to defect_insertion_aggregated_statistics.csv")

    println("\nSummary (L, q, r, n_total, mean_R, mean_absR, positive_fraction, average_sign):")
    for row in stats
        @printf(
            "  L=%-3d q=%.2f r=%-3d n=%-5d mean_R=%8.4f mean_absR=%7.4f pos_frac=%.3f avg_sign=%+.3f\n",
            row.L, row.q, row.r, row.n_total, row.mean_R, row.mean_absR,
            row.positive_fraction, row.average_sign,
        )
    end

    Ls = sort(unique(row.L for row in stats))

    println("\nGenerating plots...")
    plot_metric_vs_r(stats, Ls; ykey=:mean_R, ylabel="mean R", title="Mean signed disorder-string overlap vs r",
                      filename="defect_insertion_mean_R_vs_r", ci_lo_key=:mean_R_ci_lo, ci_hi_key=:mean_R_ci_hi)
    plot_metric_vs_r(stats, Ls; ykey=:mean_absR, ylabel="mean |R|", title="Mean |disorder-string overlap| vs r",
                      filename="defect_insertion_mean_absR_vs_r", ci_lo_key=:mean_absR_ci_lo, ci_hi_key=:mean_absR_ci_hi)
    plot_metric_vs_r(stats, Ls; ykey=:mean_finite_string_log_magnitude, ylabel="string log-magnitude",
                      title="String log-magnitude vs r", filename="defect_insertion_string_log_magnitude_vs_r",
                      ci_lo_key=:string_log_magnitude_ci_lo, ci_hi_key=:string_log_magnitude_ci_hi)
    plot_metric_vs_r(stats, Ls; ykey=:positive_fraction, ylabel="fraction with R > 0", title="Positive fraction vs r",
                      filename="defect_insertion_positive_fraction_vs_r")
    plot_metric_vs_r(stats, Ls; ykey=:average_sign, ylabel="average sign(R)", title="Average sign vs r",
                      filename="defect_insertion_average_sign_vs_r")
    plot_string_log_magnitude_vs_log_r(stats, Ls)
    plot_half_system_vs_L(stats, Ls)

    println("✓ Plots generated:")
    for name in (
        "defect_insertion_mean_R_vs_r", "defect_insertion_mean_absR_vs_r",
        "defect_insertion_string_log_magnitude_vs_r",
        "defect_insertion_string_log_magnitude_vs_log_r",
        "defect_insertion_positive_fraction_vs_r", "defect_insertion_average_sign_vs_r",
        "defect_insertion_half_system_string_log_magnitude_vs_L",
        "defect_insertion_half_system_string_log_magnitude_vs_log_L",
    )
        println("  - $name.pdf")
    end

    print_interpretation_guide()
end

main()
