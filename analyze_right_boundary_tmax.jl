#!/usr/bin/env julia

"""
Analyze right-boundary finite-time test results.

Loads JSON outputs for L in {16,24,32}, P in {0.25,0.30,0.35,0.40,0.45},
and T_max_factor in {4,8}. Generates right-boundary-style plots and a direct
comparison table for B(T_max=8L) versus B(T_max=4L).
"""

using JSON
using Statistics
using DataFrames
using CSV
using Plots
using LaTeXStrings

function infer_t_factor(file::AbstractString, row::AbstractDict)
    if haskey(row, "T_max_factor") && !ismissing(row["T_max_factor"])
        return Int(row["T_max_factor"])
    end

    m = match(r"_Tf(\d+)_", basename(file))
    if m !== nothing
        return parse(Int, m.captures[1])
    end

    if haskey(row, "T_max") && haskey(row, "L")
        L = Int(row["L"])
        T_max = Int(row["T_max"])
        if L > 0
            return Int(round(T_max / L))
        end
    end

    return -1
end

function load_results(results_dir)
    rows = Dict{String,Any}[]

    for L in [16, 24, 32]
        for Tf in [4, 8]
            subdir = joinpath(results_dir, "L$L", "Tf$Tf")
            if !isdir(subdir)
                continue
            end

            json_files = filter(
                f -> contains(f, "right_boundary_tmax_L$(L)_Tf$(Tf)_lx0.70_lzz0.00_Px") && endswith(f, ".json"),
                readdir(subdir, join=true),
            )

            for file in json_files
                try
                    parsed = JSON.parsefile(file)
                    entries = parsed isa AbstractVector ? parsed : [parsed]

                    for entry in entries
                        if !(entry isa AbstractDict)
                            continue
                        end

                        row = Dict{String,Any}(String(k) => v for (k, v) in entry)
                        row["_source_file"] = file
                        row["T_max_factor"] = infer_t_factor(file, row)
                        push!(rows, row)
                    end
                catch e
                    @warn "Could not read JSON file" file=file exception=(e, catch_backtrace())
                end
            end
        end
    end

    isempty(rows) && return DataFrame()

    all_keys = Set{String}()
    for row in rows
        union!(all_keys, keys(row))
    end

    for row in rows
        for key in all_keys
            if !haskey(row, key)
                row[key] = missing
            end
        end
    end

    return DataFrame(rows)
end

function compute_statistics(df)
    gdf = groupby(df, [:L, :lambda_x, :lambda_zz, :P_x, :T_max_factor])

    stats = combine(gdf) do group
        (
            B_mean = mean(group.B),
            B_std = std(group.B),
            B_sem = std(group.B) / sqrt(nrow(group)),
            M2_mean = mean(group.M2_bar),
            M4_mean = mean(group.M4_bar),
            purity_mean = mean(group.purity_bar),
            n_samples = nrow(group),
            time_mean = mean(group.time_seconds),
        )
    end

    sort!(stats, [:L, :P_x, :T_max_factor])
    return stats
end

function plot_binder_by_tfactor(stats)
    colors = Dict(4 => :blue, 8 => :red)
    markers = Dict(4 => :circle, 8 => :diamond)

    for L in sort(unique(stats.L))
        p = plot(
            xlabel=L"P_x = P_{zz}",
            ylabel="Binder Parameter (Renyi-2)",
            title="Right boundary finite-time test, L=$L",
            legend=:outertopright,
            grid=true,
            size=(800, 600),
            dpi=300,
        )

        for Tf in [4, 8]
            d = filter(row -> row.L == L && row.T_max_factor == Tf, stats)
            if nrow(d) == 0
                continue
            end

            plot!(
                p,
                d.P_x,
                d.B_mean,
                yerr=d.B_sem,
                label="T_max=$(Tf)L",
                color=colors[Tf],
                marker=markers[Tf],
                markersize=6,
                linewidth=2,
            )
        end

        hline!(p, [2 / 3], linestyle=:dash, color=:black, label="B = 2/3", linewidth=2, alpha=0.7)

        savefig(p, "right_boundary_tmax_L$(L)_renyi2_binder_vs_px.pdf")
        savefig(p, "right_boundary_tmax_L$(L)_renyi2_binder_vs_px.png")
    end
end

function plot_moment_comparison(stats)
    for metric in ["M2_mean", "M4_mean"]
        ylabel_txt = metric == "M2_mean" ? L"M_2" : L"M_4"

        p = plot(
            xlabel=L"P_x = P_{zz}",
            ylabel=ylabel_txt,
            title="Right boundary finite-time test: $(metric)",
            legend=:outertopright,
            grid=true,
            size=(900, 650),
            dpi=300,
        )

        for L in sort(unique(stats.L))
            for Tf in [4, 8]
                d = filter(row -> row.L == L && row.T_max_factor == Tf, stats)
                if nrow(d) == 0
                    continue
                end

                plot!(
                    p,
                    d.P_x,
                    d[!, Symbol(metric)],
                    label="L=$L, T=$(Tf)L",
                    linewidth=2,
                )
            end
        end

        if metric == "M2_mean"
            savefig(p, "right_boundary_tmax_M2_vs_px.pdf")
            savefig(p, "right_boundary_tmax_M2_vs_px.png")
        else
            savefig(p, "right_boundary_tmax_M4_vs_px.pdf")
            savefig(p, "right_boundary_tmax_M4_vs_px.png")
        end
    end
end

function build_8l_vs_4l_table(stats)
    t4 = filter(row -> row.T_max_factor == 4, stats)
    t8 = filter(row -> row.T_max_factor == 8, stats)

    rename!(t4, Dict(
        :B_mean => :B_mean_4L,
        :B_sem => :B_sem_4L,
        :M2_mean => :M2_mean_4L,
        :M4_mean => :M4_mean_4L,
        :n_samples => :n_samples_4L,
    ))

    rename!(t8, Dict(
        :B_mean => :B_mean_8L,
        :B_sem => :B_sem_8L,
        :M2_mean => :M2_mean_8L,
        :M4_mean => :M4_mean_8L,
        :n_samples => :n_samples_8L,
    ))

    keep4 = select(t4, [:L, :P_x, :B_mean_4L, :B_sem_4L, :M2_mean_4L, :M4_mean_4L, :n_samples_4L])
    keep8 = select(t8, [:L, :P_x, :B_mean_8L, :B_sem_8L, :M2_mean_8L, :M4_mean_8L, :n_samples_8L])

    comp = innerjoin(keep8, keep4, on=[:L, :P_x])

    comp.B_diff_8L_minus_4L = comp.B_mean_8L .- comp.B_mean_4L
    comp.B_combined_sem = sqrt.(comp.B_sem_8L .^ 2 .+ comp.B_sem_4L .^ 2)
    comp.B_diff_zscore = comp.B_diff_8L_minus_4L ./ comp.B_combined_sem
    comp.same_within_2sigma = abs.(comp.B_diff_zscore) .<= 2.0

    sort!(comp, [:L, :P_x])
    return comp
end

function print_comparison_summary(comp)
    println("\n=== 8L vs 4L Binder Comparison ===")

    for L in sort(unique(comp.L))
        println("L = $L")
        sub = filter(r -> r.L == L, comp)
        for r in eachrow(sub)
            status = r.same_within_2sigma ? "same (<=2sigma)" : "different (>2sigma)"
            println(
                "  P=$(round(r.P_x, digits=2)): B8=$(round(r.B_mean_8L, digits=4)), " *
                "B4=$(round(r.B_mean_4L, digits=4)), dB=$(round(r.B_diff_8L_minus_4L, digits=4)), " *
                "z=$(round(r.B_diff_zscore, digits=2)) => $status"
            )
        end
    end

    n_same = count(comp.same_within_2sigma)
    n_total = nrow(comp)
    println("\nOverall: $n_same / $n_total points are consistent within 2 sigma.")
end

function main()
    dirs = filter(d -> startswith(d, "right_boundary_tmax_results_") && isdir(d), readdir())
    if isempty(dirs)
        println("ERROR: No results directory found matching 'right_boundary_tmax_results_*'")
        return
    end

    results_dir = last(sort(dirs))
    println("Loading results from: $results_dir")

    df = load_results(results_dir)
    if nrow(df) == 0
        println("ERROR: No data loaded.")
        return
    end

    println("Loaded $(nrow(df)) data points")

    stats = compute_statistics(df)
    println("Computed $(nrow(stats)) grouped statistics rows")

    plot_binder_by_tfactor(stats)
    plot_moment_comparison(stats)

    comp = build_8l_vs_4l_table(stats)
    CSV.write("right_boundary_tmax_8L_vs_4L_comparison.csv", comp)
    print_comparison_summary(comp)

    println("\nOutputs generated:")
    println("  - right_boundary_tmax_L{16,24,32}_renyi2_binder_vs_px.(pdf|png)")
    println("  - right_boundary_tmax_M2_vs_px.(pdf|png)")
    println("  - right_boundary_tmax_M4_vs_px.(pdf|png)")
    println("  - right_boundary_tmax_8L_vs_4L_comparison.csv")
end

main()
