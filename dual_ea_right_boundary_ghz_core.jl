"""
dual_ea_right_boundary_ghz_core.jl

Dual Edwards-Anderson (EA) Binder parameter on the right boundary with a GHZ INITIAL STATE,
|GHZ> = (|00...0> + |11...1>)/sqrt(2), instead of the all-up product state used by
`dual_ea_right_boundary_core.jl`.

Nothing existing is modified. This file only adds the GHZ initial state and a trajectory
driver around code that already exists:
  * dynamics helpers : `rb_*` from `renyi2_right_boundary_core.jl` (weak X layer, X / ZZ
                       dephasing channels, trace normalisation, ...), exactly as used by
                       `rb_evolve_right_edge_one_trial`;
  * observable       : `dual_ea_moments` and the statistics helpers from
                       `dual_ea_right_boundary_core.jl` (state independent: it only reads
                       the doubled MPS).

`dual_ea_ghz_evolve_right_edge_one_trial` below is `rb_evolve_right_edge_one_trial` line for
line (same layer order: weak X -> X dephasing -> ZZ dephasing, same RNG use, same
open-chain conventions) except for the initial state; the GHZ state is built exactly as in
`renyi2_right_boundary_ghz_core.jl` (`rb_ghz_initial_state`). Per-trajectory seeds are drawn
as in `rb_run_right_edge_point` / `renyi2_right_boundary_ghz_core.jl`, so for the same `seed`
this scan analyses the SAME Born trajectories as the Renyi-2 right-boundary GHZ scan.

Observable and measurement time are those of the all-up dual-EA scan (endpoints a = 2..L,
N = L-1, ordinary trace expectations squared, measured after the last ZZ-dephasing step of
layer t = T_max).

Exact limits used as tests (see `dual_ea_right_boundary_ghz_validate.jl`):
  * t = 0: every X-string W_ab with 2 <= a < b <= L acts on at most sites 2..L-1 and never on
    all L sites, so E(W_ab) = <GHZ|W_ab|GHZ> = 0 and only the diagonal a = b survives:
    M2 = 1/N, M4 = (3N^2 - 2N)/N^4, B = 2/(3N).
  * q = 1/2: the final ZZ-dephasing layer removes every non-empty string for EVERY state,
    so the same values hold for every trajectory.
  * q = 0 and long times: the repeated weak X measurements collapse the GHZ superposition
    onto a single X-eigenstate product state, so M2 = M4 = 1 and B = 2/3.
"""

using Random
using Statistics
using LinearAlgebra
using ITensors
using ITensorMPS

# The dual-EA core includes renyi2_right_boundary_core.jl (guarded), so rb_* helpers exist.
if !isdefined(Main, :dual_ea_moments)
    include(joinpath(@__DIR__, "dual_ea_right_boundary_core.jl"))
end

"""
    dual_ea_ghz_initial_state(sites_in)

Unnormalised doubled-MPS density matrix of |GHZ><GHZ|: the sum of the four bra/ket
computational-basis sectors (00, 01, 10, 11), each a product state across bra/ket pairs.
Same construction as `rb_ghz_initial_state` in `renyi2_right_boundary_ghz_core.jl`.
"""
function dual_ea_ghz_initial_state(sites_in)
    sites = collect(sites_in)
    rho_00 = rb_product_density_mps(sites, [1.0, 0.0, 0.0, 0.0])  # |0...0><0...0|
    rho_01 = rb_product_density_mps(sites, [0.0, 1.0, 0.0, 0.0])  # |0...0><1...1|
    rho_10 = rb_product_density_mps(sites, [0.0, 0.0, 1.0, 0.0])  # |1...1><0...0|
    rho_11 = rb_product_density_mps(sites, [0.0, 0.0, 0.0, 1.0])  # |1...1><1...1|
    return (rho_00 + rho_01 + rho_10 + rho_11) / 2
end

"""
    dual_ea_ghz_evolve_right_edge_one_trial(L; lambda_x, q, T_max, maxdim, cutoff, seed)

One Born-sampled right-edge (lambda_zz = 0) trajectory starting from the GHZ state; identical
to `rb_evolve_right_edge_one_trial` except for the initial state.
"""
function dual_ea_ghz_evolve_right_edge_one_trial(
    L::Int;
    lambda_x::Float64,
    q::Float64,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    seed::Int,
    initial_state=dual_ea_ghz_initial_state,   # swappable so tests can run the same code path from all-up
)
    @assert L >= 2
    @assert 0.0 <= q <= 0.5
    @assert T_max >= 0

    rng = MersenneTwister(seed)

    sites = siteinds("Qubit", 2L)
    trace_bra = rb_bell_state(sites)
    rho = rb_trace_normalize(initial_state(sites), trace_bra)

    x_dephasing_gates = rb_build_x_dephasing_gates(sites, L, q)
    zz_dephasing_gates = rb_build_zz_dephasing_gates(sites, L, q)

    max_trace_error = abs(rb_doubled_trace(rho, trace_bra) - 1.0)
    max_cross_site_bond = rb_max_interphysical_linkdim(rho, L)

    for _ in 1:T_max
        # 1. Weak X measurements.
        rho, _ = rb_apply_weak_x_layer(
            rho, sites, L, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )

        # 2. X dephasing, q_x = q.
        rho = rb_apply_channel_layer(
            rho, x_dephasing_gates, trace_bra; maxdim=maxdim, cutoff=cutoff,
        )

        # 3. Weak ZZ measurements are absent because lambda_zz = 0.

        # 4. ZZ dephasing remains present, q_zz = q.
        rho = rb_apply_channel_layer(
            rho, zz_dephasing_gates, trace_bra; maxdim=maxdim, cutoff=cutoff,
        )

        max_trace_error = max(
            max_trace_error, abs(rb_doubled_trace(rho, trace_bra) - 1.0),
        )
        max_cross_site_bond = max(
            max_cross_site_bond, rb_max_interphysical_linkdim(rho, L),
        )
    end

    return (
        rho=rho,
        sites=sites,
        trace_bra=trace_bra,
        max_trace_error=max_trace_error,
        max_interphysical_linkdim=max_cross_site_bond,
    )
end

"""
    dual_ea_ghz_run_right_edge_point(L, q; lambda_x, lambda_zz, ntrials, T_max, maxdim, cutoff, seed, nboot, verbose)

Run `ntrials` Born-sampled GHZ-initial-state trajectories at (L, q) and evaluate the dual-EA
moments on each final density matrix. Same as `dual_ea_run_right_edge_point` apart from the
initial state. The two Binder conventions are kept strictly separate:

  B_mean_of_trials  = mean_m B_m
  B_ratio_of_means  = 1 - mean(M4_m) / (3 mean(M2_m)^2)
"""
function dual_ea_ghz_run_right_edge_point(
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

    for trial in 1:ntrials
        trial_seeds[trial] = Int(rand(master_rng, UInt32))

        evolved = dual_ea_ghz_evolve_right_edge_one_trial(
            L; lambda_x=lambda_x, q=q, T_max=T_max, maxdim=maxdim,
            cutoff=cutoff, seed=trial_seeds[trial],
        )

        t0 = time()
        obs = dual_ea_moments(evolved.rho, L)
        obs_seconds[trial] = time() - t0

        M2s[trial], M4s[trial], Bs[trial] = obs.M2, obs.M4, obs.B
        linkdims[trial] = evolved.max_interphysical_linkdim
        trace_errors[trial] = evolved.max_trace_error

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
        B_mean_of_trials=n_finite > 0 ? mean(Bs[finite]) : NaN,
        B_std_of_trials=n_finite < 2 ? NaN : std(Bs[finite]; corrected=true),
        B_se_of_trials=dual_ea_se(Bs[finite]),
        B_ratio_of_means=dual_ea_binder_from_moments(M2_mean, M4_mean),
        B_ratio_of_means_bootstrap_se=dual_ea_bootstrap_ratio_se(M2s, M4s; nboot=nboot, rng=boot_rng),
        ntrials_requested=ntrials, n_completed=ntrials, n_finite_B=n_finite,
        max_interphysical_linkdim=maximum(linkdims),
        max_trace_error=maximum(trace_errors),
        obs_seconds_mean=mean(obs_seconds), obs_seconds_max=maximum(obs_seconds),
    )
end
