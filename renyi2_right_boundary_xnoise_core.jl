"""
renyi2_right_boundary_xnoise_core.jl

Right-boundary Renyi-2 Binder measured AFTER THE LAST X-DEPHASING step (before the last
ZZ-dephasing layer), i.e. at the same point of the circuit where the left-boundary EA code
(`evolve_density_matrix_one_trial`, "ρ_after_X_noise") measures.

This file is a thin layer over the production `renyi2_right_boundary_core.jl`, which is
included unchanged and whose helpers are ALL reused (initial state, Born-sampled weak X
layer, X/ZZ channel gates, trace normalisation, replica-overlap observable
`rb_renyi2_binder_one_trajectory`). Nothing there is modified or copied.

Circuit (identical to the production core; lambda_x = 0.7, lambda_zz = 0, q_x = q_zz = q,
open chain, initial state |up..up><up..up|), per step t = 1..T_max:
  1. weak X measurement (exact Born sampling)
  2. X dephasing channel
  3. weak ZZ measurement   -- absent (lambda_zz = 0)
  4. ZZ dephasing channel
Production observable time:  after step 4 of t = T_max   (`rb_evolve_right_edge_one_trial`).
THIS FILE'S observable time: after step 2 of t = T_max   (steps 3-4 of the final step are skipped).

Only the final ZZ-dephasing layer is skipped: the first T_max-1 steps and the final step's
measurement + X dephasing are the same operations in the same order with the same RNG
(`MersenneTwister(seed)` is consumed only by the Born sampling), so for the same seed the state
returned here equals the production trajectory just before its last ZZ layer, and applying
that ZZ layer to it reproduces the production state bit-for-bit (tested in
`renyi2_right_boundary_xnoise_validate.jl`). Per-point and per-trajectory seeds are drawn
exactly as in `rb_run_right_edge_point`, so with the same `--seed` the X-noise scan and the
production (after-ZZ) scan analyse the SAME Born trajectories, observed at two times.

Neither observation time is "the corrected one"; the comparison probes the temporal boundary
of the finite-time circuit. Observable definition, per-trajectory Binder averaging,
trace normalisation and validity handling are the production ones.
"""

using Random
using Statistics
using LinearAlgebra
using ITensors
using ITensorMPS

# Production core, included unchanged (guarded: it defines `const`s).
if !isdefined(Main, :rb_evolve_right_edge_one_trial)
    include(joinpath(@__DIR__, "renyi2_right_boundary_core.jl"))
end

"""
    rb_xn_evolve_right_edge_one_trial(L; lambda_x, q, T_max, maxdim, cutoff, seed)

Born-sampled right-edge trajectory for `T_max` steps; the last step stops after the X-dephasing
channel. Returns the state at that point (`rho`) plus the same diagnostics as the production
evolution (trace error, maximal inter-physical bond dimension, measurement record).
"""
function rb_xn_evolve_right_edge_one_trial(
    L::Int;
    lambda_x::Float64,
    q::Float64,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    seed::Int,
)
    @assert L >= 2
    @assert 0.0 <= q <= 0.5
    @assert T_max >= 1

    rng = MersenneTwister(seed)

    sites = siteinds("Qubit", 2L)
    trace_bra = rb_bell_state(sites)
    rho = rb_trace_normalize(rb_initial_repetition_code_state(sites), trace_bra)

    x_dephasing_gates = rb_build_x_dephasing_gates(sites, L, q)
    zz_dephasing_gates = rb_build_zz_dephasing_gates(sites, L, q)

    max_trace_error = abs(rb_doubled_trace(rho, trace_bra) - 1.0)
    max_cross_site_bond = rb_max_interphysical_linkdim(rho, L)
    outcomes = zeros(Int, T_max, L)

    for t in 1:T_max
        # 1. Weak X measurements.
        rho, outcomes[t, :] = rb_apply_weak_x_layer(
            rho, sites, L, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )

        # 2. X dephasing, q_x = q.
        rho = rb_apply_channel_layer(
            rho, x_dephasing_gates, trace_bra; maxdim=maxdim, cutoff=cutoff,
        )

        max_trace_error = max(
            max_trace_error, abs(rb_doubled_trace(rho, trace_bra) - 1.0),
        )
        max_cross_site_bond = max(
            max_cross_site_bond, rb_max_interphysical_linkdim(rho, L),
        )

        # Observation point: after the X dephasing of the LAST step.
        t == T_max && break

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
        outcomes=outcomes,
        seed=seed,
        max_trace_error=max_trace_error,
        max_interphysical_linkdim=max_cross_site_bond,
    )
end

"""
    rb_xn_run_right_edge_point(L, q; lambda_x, lambda_zz, ntrials, T_max, maxdim, cutoff,
                               obs_maxdim, obs_cutoff, seed, nboot)

Same as `rb_run_right_edge_point` (same seed derivation, same observable, same summary
fields), except that the observable is evaluated after the last X-dephasing step. Additionally
returns the per-trajectory M2, M4, B2, purity and the trial seeds. Trajectories whose
observable is non-finite are NOT hidden: they are counted in `n_invalid` and stored as NaN.
"""
function rb_xn_run_right_edge_point(
    L::Int,
    q::Float64;
    lambda_x::Float64,
    lambda_zz::Float64,
    ntrials::Int,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    obs_maxdim::Int,
    obs_cutoff::Float64,
    seed::Int,
    nboot::Int,
)
    @assert ntrials >= 1
    @assert 0.0 <= q <= 0.5
    @assert isapprox(lambda_zz, 0.0; atol=1e-9) (
        "This module only implements the right edge (lambda_zz = 0); " *
        "weak ZZ measurements are not implemented."
    )

    master_rng = MersenneTwister(seed)

    trial_seeds = Vector{Int}(undef, ntrials)
    M2s = fill(NaN, ntrials)
    M4s = fill(NaN, ntrials)
    B2s = fill(NaN, ntrials)
    purities = fill(NaN, ntrials)
    cross_site_dims = Vector{Int}(undef, ntrials)
    trace_errors = Vector{Float64}(undef, ntrials)

    for trial in 1:ntrials
        trial_seeds[trial] = Int(rand(master_rng, UInt32))

        evolved = rb_xn_evolve_right_edge_one_trial(
            L; lambda_x=lambda_x, q=q, T_max=T_max, maxdim=maxdim,
            cutoff=cutoff, seed=trial_seeds[trial],
        )

        cross_site_dims[trial] = evolved.max_interphysical_linkdim
        trace_errors[trial] = evolved.max_trace_error

        observable = rb_renyi2_binder_one_trajectory(
            evolved.rho, evolved.sites, L; maxdim=obs_maxdim, cutoff=obs_cutoff,
        )
        M2s[trial] = observable.M2
        M4s[trial] = observable.M4
        B2s[trial] = observable.B2
        purities[trial] = observable.purity
    end

    valid = [
        isfinite(M2s[i]) && isfinite(M4s[i]) && isfinite(B2s[i]) &&
        isfinite(purities[i]) && M2s[i] > 0 && purities[i] > 0
        for i in eachindex(B2s)
    ]

    n_valid = count(valid)
    n_invalid = ntrials - n_valid

    common = (
        ntrials=ntrials, n_valid=n_valid, n_invalid=n_invalid,
        max_interphysical_linkdim=maximum(cross_site_dims),
        max_trace_error=maximum(trace_errors),
        trial_seeds=trial_seeds, M2_per_trajectory=M2s, M4_per_trajectory=M4s,
        B2_per_trajectory=B2s, purity_per_trajectory=purities,
    )

    if n_valid < 2
        return merge((
            B=NaN, B_mean_of_trials=NaN, B_std_of_trials=NaN,
            M2_bar=NaN, M4_bar=NaN, purity_bar=NaN,
            B_bootstrap_se=NaN, B_ci_low=NaN, B_ci_high=NaN,
            B2_ratio_of_mean_moments=NaN,
        ), common)
    end

    valid_M2 = M2s[valid]
    valid_M4 = M4s[valid]
    valid_B2 = B2s[valid]
    valid_purity = purities[valid]

    # Primary estimate: trajectory-averaged Binder (production convention).
    B2_mean = mean(valid_B2)

    boot = rb_bootstrap_mean(
        valid_B2; nboot=nboot, rng=MersenneTwister(seed + 918273),
    )

    # Diagnostic only; not the primary average.
    M2_mean = mean(valid_M2)
    M4_mean = mean(valid_M4)

    return merge((
        B=B2_mean,
        B_mean_of_trials=B2_mean,
        B_std_of_trials=std(valid_B2; corrected=true),
        M2_bar=M2_mean,
        M4_bar=M4_mean,
        purity_bar=mean(valid_purity),
        B_bootstrap_se=boot.standard_error,
        B_ci_low=boot.ci_low,
        B_ci_high=boot.ci_high,
        B2_ratio_of_mean_moments=rb_binder_from_moments(M2_mean, M4_mean),
    ), common)
end
