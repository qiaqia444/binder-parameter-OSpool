using JSON
using Statistics
using Plots

files = if !isempty(ARGS)
    ARGS
else
    result_dirs = filter(
        d -> startswith(d, "defect_insertion_results_") && isdir(d),
        readdir(),
    )
    isempty(result_dirs) ? String[] : String[
        file for directory in [last(sort(result_dirs))]
        for L in [16, 24, 32]
        for file in readdir(joinpath(directory, "L$L"), join=true)
        if endswith(file, ".json")
    ]
end
rows = Any[]
for file in files
    append!(rows, JSON.parsefile(file))
end
pairs = filter(row -> get(row, "row_type", "") == "pair", rows)
averages = filter(row -> get(row, "row_type", "") == "trajectory_average", rows)
@assert !isempty(averages) "Pass one or more defect_insertion JSON files."

groups = Dict{Tuple{Int,Float64},Vector{Float64}}()
for row in averages
    key = (Int(row["L"]), Float64(row["control_parameter"]))
    push!(get!(groups, key, Float64[]), Float64(row["deltaF"]))
end
summary = [(L=k[1], control_parameter=k[2], deltaF_mean=mean(filter(isfinite, v)), deltaF_se=length(filter(isfinite, v))>1 ? std(filter(isfinite, v))/sqrt(length(filter(isfinite, v))) : NaN, n_finite=length(filter(isfinite, v)), n_total=length(v), deltaF_over_r=mean(filter(isfinite, v))/(k[1]/2)) for (k,v) in sort(collect(groups))]
for row in summary
    println(row)
end
if !isempty(pairs)
    println("Saved pair rows: ", length(pairs), "; trajectory rows: ", length(averages))
end
for P in sort(unique(row.control_parameter for row in summary))
    selected = filter(row -> row.control_parameter == P, summary)
    L = [row.L for row in selected]
    F = [row.deltaF_mean for row in selected]
    p = plot(L, F; marker=:circle, xlabel="L", ylabel="mean deltaF", label="P=$(P)")
    savefig(p, "defect_insertion_deltaF_vs_L_P$(P).pdf")
    p = plot(log.(L), F; marker=:circle, xlabel="log L", ylabel="mean deltaF", label="P=$(P)")
    savefig(p, "defect_insertion_deltaF_vs_logL_P$(P).pdf")
end
