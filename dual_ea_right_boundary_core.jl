"""
dual_ea_right_boundary_core.jl

Dual Edwards-Anderson (EA) Binder parameter on the right boundary, built on top of
the EXISTING right-boundary doubled-MPS evolution in `renyi2_right_boundary_core.jl`
(`rb_evolve_right_edge_one_trial`). No evolution code is duplicated or modified and the
existing Renyi-2 / EA implementations are untouched.

--------------------------------------------------------------------------------
Reused dynamics (open chain, lambda_x = 0.7, lambda_zz = 0, q_x = q_zz = q)
--------------------------------------------------------------------------------
* Doubled MPS with 2L sites, INTERLEAVED layout: MPS site 2i-1 = bra index of
  physical site i, MPS site 2i = ket index of physical site i.
* Initial state |up..up><up..up| (all Z = +1), trace normalised.
* Open boundary conditions: ZZ dephasing acts on bonds (i, i+1), i = 1..L-1.
* One layer = (1) weak X measurements with exact Born sampling over both Kraus
  branches, (2) X dephasing CHANNEL (superoperator gate, never sampled),
  (3) weak ZZ measurement (absent, lambda_zz = 0), (4) ZZ dephasing CHANNEL.
* T_max = T_max_factor * L layers (default factor 4, same as run_right_boundary_scan.jl).
* Per-trajectory seeds are drawn exactly as in `rb_run_right_edge_point`
  (`Int(rand(MersenneTwister(seed), UInt32))` once per trajectory), so for the same
  `seed` the dual-EA scan analyses the SAME Born trajectories as the Renyi-2 scan.

WHEN THE OBSERVABLE IS MEASURED: on the state returned by the evolution, i.e. after
the last step (4) (ZZ dephasing) of the final layer t = T_max. NOTE this differs from
the left-boundary EA code (`evolve_density_matrix_one_trial`), which returns the state
after step 2 (X dephasing). It matters here: ZZ dephasing on bond (i,i+1) multiplies
<W> by (1-2q) whenever the X-string W contains exactly one of the sites i, i+1.

--------------------------------------------------------------------------------
Observable
--------------------------------------------------------------------------------
For one Born-sampled trajectory density matrix rho_m,

    E_m(O) = Tr(rho_m O) / Tr(rho_m)                (ordinary trace expectation)

Dual endpoints a = 2..L (open chain, N = L-1 internal endpoints), strings

    W_ab = prod_{k=min(a,b)}^{max(a,b)-1} X_k,      W_aa = I   (X_k acts on site k in 2..L-1).

    G2(a,b)     = E_m(W_ab)^2
    G4(a,b,c,d) = E_m(W_uv W_wx)^2,   (u<=v<=w<=x) = sort(a,b,c,d)
    M2_m = sum_{a,b}   G2 / N^2
    M4_m = sum_{a,b,c,d} G4 / N^4              (all ordered tuples, repeats included)
    B_m  = 1 - M4_m / (3 M2_m^2)

W_uv W_wx is the same operator for every pairing of the four endpoints (the X_k commute
and tau_a tau_b = W_ab), so G4 is a symmetric function of the multiset {a,b,c,d}.
Sorted tuples are therefore summed with their exact permutation multiplicity
4!/prod(m!) (m = multiplicities of equal endpoints); repeated endpoints and empty strings
are handled exactly (no special cases). This diagnostic uses ordinary trace expectations
squared: no Tr(rho O rho O)/Tr(rho^2), no Hilbert-Schmidt normalisation, no purity
denominator, and no purity-based rejection of states.

--------------------------------------------------------------------------------
Contraction (exact, no truncation)
--------------------------------------------------------------------------------
Tr(rho O_1 ... O_L) = <<v_{O_1} ... v_{O_L} | rho>> is a (linear!) overlap of the MPS with
a PRODUCT state, so every environment is a bond-dimension-sized VECTOR, not a matrix.
Contracting the bra/ket tensor pair (2i-1, 2i) of each physical site with the vectorised
identity / X gives per-site transfer matrices MI[i], MX[i] (chi_{i-1} x chi_i). All
required expectations are then products
    lid[u] . MX[u..v-1] . MI[v..w-1] . MX[w..x-1] . MI[x..L] . 1.
Left vectors A[u][v], the running identity-propagated vector, and right vectors
R1[w][x] are built incrementally and reused across all endpoint tuples. Cost per
trajectory ~ (N^3/6) chi^2 + (N^4/24) chi flops, memory O(N^2 chi); no dense operator
or density matrix and no exponentially large object is formed. See
`dual_ea_cost_estimate`. The input MPS is only read (never orthogonalised, truncated
or otherwise mutated) and no random numbers are drawn.
"""

using Random
using Statistics
using LinearAlgebra
using ITensors
using ITensorMPS

# Reuse the established right-boundary evolution (guarded so it can also be included
# by a caller that already loaded it; the core defines `const`s).
if !isdefined(Main, :rb_evolve_right_edge_one_trial)
    include(joinpath(@__DIR__, "renyi2_right_boundary_core.jl"))
end

const DUAL_EA_PAULI_X = Float64[0 1; 1 0]
const DUAL_EA_IDENTITY_2 = Matrix{Float64}(I, 2, 2)

"Dual endpoints a = 2..L for an open chain of L sites (N = L-1 of them)."
dual_ea_endpoints(L::Int) = 2:L

# ============================================================
# Transfer matrices of the doubled MPS
# ============================================================
"""
    dual_ea_transfer_matrices(rho::MPS, L) -> (MI, MX)

Per physical site i, `MI[i]` / `MX[i]` are the (chi_{i-1} x chi_i) matrices obtained by
contracting the interleaved (bra, ket) tensors `rho[2i-1] * rho[2i]` with the vectorised
identity / Pauli X. For real symmetric X (and identity) the result does not depend on
which of the two physical indices is called row or column. Boundary matrices have a
trivial (dimension-1) outer link.
"""
function dual_ea_transfer_matrices(rho::MPS, L::Int)
    @assert length(rho) == 2L "Expected a doubled MPS with 2L = $(2L) sites, got $(length(rho))"
    sites = siteinds(rho)

    mats = Vector{Tuple{Matrix,Matrix}}(undef, L)
    for i in 1:L
        s_bra, s_ket = sites[2i - 1], sites[2i]
        pair = rho[2i - 1] * rho[2i]

        # An absent link (chain boundary, or a hand-built product state) has dimension 1.
        left = i == 1 ? nothing : linkind(rho, 2i - 2)
        right = i == L ? nothing : linkind(rho, 2i)
        dl = isnothing(left) ? 1 : dim(left)
        dr = isnothing(right) ? 1 : dim(right)
        outer = Index[l for l in (left, right) if !isnothing(l)]

        to_matrix(O) = reshape(array(pair * ITensor(O, s_bra, s_ket), outer...), dl, dr)

        mats[i] = (to_matrix(DUAL_EA_IDENTITY_2), to_matrix(DUAL_EA_PAULI_X))
    end

    T = mapreduce(m -> promote_type(eltype(m[1]), eltype(m[2])), promote_type, mats)
    MI = [convert(Matrix{T}, m[1]) for m in mats]
    MX = [convert(Matrix{T}, m[2]) for m in mats]
    return MI, MX
end

"""
    dual_ea_pattern_expectation(MI, MX, flips::AbstractVector{Bool}) -> unnormalised Tr(rho prod_k X_k^{flips[k]})

Direct (non-incremental) chain contraction for an arbitrary X-pattern; used for
cross-checks and for single-string evaluation.
"""
function dual_ea_pattern_expectation(MI, MX, flips::AbstractVector{Bool})
    L = length(MI)
    @assert length(flips) == L
    v = ones(eltype(MI[1]), 1)
    for i in 1:L
        M = flips[i] ? MX[i] : MI[i]
        v = transpose(M) * v
    end
    return real(v[1])
end

"X-pattern of W_uv W_wx for endpoints (any order): sites in [u,v) and [w,x) after sorting."
function dual_ea_string_pattern(L::Int, endpoints::NTuple{N,Int}) where {N}
    flips = falses(L)
    sorted = sort(collect(endpoints))
    @assert iseven(length(sorted))
    for p in 1:2:length(sorted)
        for k in sorted[p]:(sorted[p + 1] - 1)
            flips[k] = !flips[k]
        end
    end
    return flips
end

# Permutation multiplicity of a sorted 4-tuple: 4! / prod(run-length!).
@inline function dual_ea_multiplicity4(u::Int, v::Int, w::Int, x::Int)
    if u == x
        return 1
    elseif u == w || v == x
        return 4
    elseif u == v && w == x
        return 6
    elseif u == v || v == w || w == x
        return 12
    else
        return 24
    end
end

# ============================================================
# Exact moments of one density matrix
# ============================================================
"""
    dual_ea_moments(rho::MPS, L) -> NamedTuple

Exact dual-EA moments (M2, M4, B) of the doubled-MPS density matrix `rho`
(endpoints 2..L, see the module docstring). `rho` may be unnormalised; expectations are
divided by Tr(rho) = <bell|rho>. Only requires Tr(rho) to be finite and positive.
Returns (M2, M4, B, trace, n_endpoints, n_sorted_pairs, n_sorted_quads).
"""
function dual_ea_moments(rho::MPS, L::Int)
    @assert L >= 2
    MI, MX = dual_ea_transfer_matrices(rho, L)
    T = eltype(MI[1])
    N = L - 1

    # Right identity vectors: rid[k] = MI[k] ... MI[L] * 1  (column, length chi_{k-1}).
    rid = Vector{Vector{T}}(undef, L + 1)
    rid[L + 1] = ones(T, 1)
    for k in L:-1:1
        rid[k] = MI[k] * rid[k + 1]
    end
    tr = real(rid[1][1])
    if !isfinite(tr) || tr <= 0.0
        error("Invalid density-matrix trace: Tr(rho) = $tr")
    end

    # Left identity vectors: lid[u] = (MI[1] ... MI[u-1])^T  (column, length chi_{u-1}).
    lid = Vector{Vector{T}}(undef, L)
    lid[1] = ones(T, 1)
    for u in 2:L
        lid[u] = transpose(MI[u - 1]) * lid[u - 1]
    end

    # R1[w][:, x-w+1] = MX[w] ... MX[x-1] * rid[x] for x = w..L   (X on sites [w, x)).
    R1 = Vector{Matrix{T}}(undef, L)
    R1[L] = reshape(rid[L], :, 1)
    for w in (L - 1):-1:2
        R1[w] = hcat(rid[w], MX[w] * R1[w + 1])
    end

    sum2 = 0.0
    sum4 = 0.0
    mult2_total = 0
    mult4_total = 0
    n_pairs = 0
    n_quads = 0

    for u in 2:L
        a = lid[u]                       # left vector with X on [u, v), starting at v = u
        for v in u:L
            # Two-point term: E(W_uv) = a . rid[v] / tr, multiplicity 2 for u<v, 1 for u=v.
            e2 = real(transpose(a) * rid[v]) / tr
            m2 = u == v ? 1 : 2
            sum2 += m2 * e2^2
            mult2_total += m2
            n_pairs += 1

            # Four-point terms with first pair (u,v); propagate identity from v to w.
            c = a
            for w in v:L
                vals = transpose(R1[w]) * c     # x = w..L
                for (j, val) in enumerate(vals)
                    x = w + j - 1
                    m4 = dual_ea_multiplicity4(u, v, w, x)
                    e4 = real(val) / tr
                    sum4 += m4 * e4^2
                    mult4_total += m4
                    n_quads += 1
                end
                w < L && (c = transpose(MI[w]) * c)
            end

            v < L && (a = transpose(MX[v]) * a)
        end
    end

    @assert mult2_total == N^2 "2-point multiplicity bookkeeping failed: $mult2_total != $(N^2)"
    @assert mult4_total == N^4 "4-point multiplicity bookkeeping failed: $mult4_total != $(N^4)"

    M2 = sum2 / N^2
    M4 = sum4 / N^4
    B = M2 > 0.0 && isfinite(M2) && isfinite(M4) ? 1.0 - M4 / (3.0 * M2^2) : NaN

    return (
        M2=M2, M4=M4, B=B, trace=tr,
        n_endpoints=N, n_sorted_pairs=n_pairs, n_sorted_quads=n_quads,
    )
end

"""
    dual_ea_cost_estimate(L, chi) -> NamedTuple

Counting-based cost model of one `dual_ea_moments` call (floating point operations,
dominated by the identity-propagation matvecs and the endpoint dot products).
"""
function dual_ea_cost_estimate(L::Int, chi::Int)
    N = L - 1
    n_pairs = N * (N + 1) ÷ 2
    n_quads = binomial(N + 3, 4)
    n_propagations = binomial(N + 2, 3)          # (u<=v<=w) triples
    flops = 2.0 * n_propagations * chi^2 + 2.0 * n_quads * chi + 2.0 * n_pairs * chi^2
    memory_floats = (N * (N + 1) ÷ 2 + 3L) * chi + 2L * chi^2   # R1 columns + vectors + MI/MX
    return (
        n_endpoints=N, n_sorted_pairs=n_pairs, n_sorted_quads=n_quads,
        n_ordered_quads=N^4, approx_flops=flops, approx_memory_bytes=8.0 * memory_floats,
    )
end

# ============================================================
# Statistics across trajectories
# ============================================================
function dual_ea_binder_from_moments(M2::Real, M4::Real)
    return (isfinite(M2) && isfinite(M4) && M2 > 0.0) ? 1.0 - M4 / (3.0 * M2^2) : NaN
end

function dual_ea_se(values::AbstractVector)
    n = length(values)
    return n < 2 ? NaN : std(values; corrected=true) / sqrt(n)
end

"Bootstrap standard error of B_ratio_of_means (resampling whole trajectories)."
function dual_ea_bootstrap_ratio_se(
    M2s::AbstractVector, M4s::AbstractVector; nboot::Int, rng::AbstractRNG,
)
    n = length(M2s)
    (n < 2 || nboot < 2) && return NaN
    samples = Vector{Float64}(undef, nboot)
    for b in 1:nboot
        idx = rand(rng, 1:n, n)
        samples[b] = dual_ea_binder_from_moments(mean(M2s[idx]), mean(M4s[idx]))
    end
    return std(samples; corrected=true)
end

# ============================================================
# One right-edge parameter point (L, q)
# ============================================================
"""
    dual_ea_run_right_edge_point(L, q; lambda_x, lambda_zz, ntrials, T_max, maxdim, cutoff, seed, nboot, verbose)

Run `ntrials` Born-sampled right-edge trajectories (identical dynamics and seed
derivation as `rb_run_right_edge_point`) and evaluate the dual-EA moments on each final
density matrix. The two Binder conventions are kept strictly separate:

  B_mean_of_trials  = mean_m B_m
  B_ratio_of_means  = 1 - mean(M4_m) / (3 mean(M2_m)^2)
"""
function dual_ea_run_right_edge_point(
    L::Int,
    q::Float64;
    lambda_x::Float64,
    lambda_zz::Float64,
    ntrials::Int,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    seed::Int,
    nboot::Int=1000,
    verbose::Bool=false,
)
    @assert L >= 2
    @assert ntrials >= 1
    @assert 0.0 <= q <= 0.5
    @assert isapprox(lambda_zz, 0.0; atol=1e-9) (
        "The right-boundary evolution only implements lambda_zz = 0."
    )

    master_rng = MersenneTwister(seed)   # used ONLY to draw trajectory seeds, as in rb_run_right_edge_point

    trial_seeds = Vector{Int}(undef, ntrials)
    M2s = fill(NaN, ntrials)
    M4s = fill(NaN, ntrials)
    Bs = fill(NaN, ntrials)
    obs_seconds = zeros(ntrials)
    linkdims = zeros(Int, ntrials)
    trace_errors = zeros(ntrials)
    n_completed = 0

    for trial in 1:ntrials
        trial_seeds[trial] = Int(rand(master_rng, UInt32))

        evolved = rb_evolve_right_edge_one_trial(
            L; lambda_x=lambda_x, q=q, T_max=T_max, maxdim=maxdim,
            cutoff=cutoff, seed=trial_seeds[trial],
        )

        t0 = time()
        obs = dual_ea_moments(evolved.rho, L)
        obs_seconds[trial] = time() - t0

        M2s[trial], M4s[trial], Bs[trial] = obs.M2, obs.M4, obs.B
        linkdims[trial] = evolved.max_interphysical_linkdim
        trace_errors[trial] = evolved.max_trace_error
        n_completed += 1

        verbose && println(
            "    trajectory $trial/$ntrials: M2=$(round(obs.M2, digits=6)) " *
            "M4=$(round(obs.M4, digits=6)) B=$(round(obs.B, digits=4)) " *
            "chi=$(linkdims[trial]) obs_time=$(round(obs_seconds[trial], digits=3))s"
        )
    end

    finite = [isfinite(Bs[i]) for i in 1:ntrials]
    n_finite = count(finite)
    boot_rng = MersenneTwister(seed + 918273)   # independent of the dynamics RNG

    M2_mean = mean(M2s)
    M4_mean = mean(M4s)

    return (
        M2s=M2s, M4s=M4s, Bs=Bs, trial_seeds=trial_seeds,
        M2_mean=M2_mean, M4_mean=M4_mean,
        M2_se=dual_ea_se(M2s), M4_se=dual_ea_se(M4s),
        B_mean_of_trials=mean(Bs[finite]),
        B_std_of_trials=n_finite < 2 ? NaN : std(Bs[finite]; corrected=true),
        B_se_of_trials=dual_ea_se(Bs[finite]),
        B_ratio_of_means=dual_ea_binder_from_moments(M2_mean, M4_mean),
        B_ratio_of_means_bootstrap_se=dual_ea_bootstrap_ratio_se(M2s, M4s; nboot=nboot, rng=boot_rng),
        ntrials_requested=ntrials, n_completed=n_completed, n_finite_B=n_finite,
        max_interphysical_linkdim=maximum(linkdims),
        max_trace_error=maximum(trace_errors),
        obs_seconds_mean=mean(obs_seconds), obs_seconds_max=maximum(obs_seconds),
    )
end
