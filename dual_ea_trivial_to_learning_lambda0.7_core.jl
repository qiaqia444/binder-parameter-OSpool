"""
dual_ea_trivial_to_learning_lambda0.7_core.jl

Dual Edwards-Anderson (EA) Binder parameter for the "trivial-to-learning" slice at
lambda = 0.7 (delta = 0.7: lambda_x = delta*lambda = 0.49, lambda_zz = delta*(1-lambda) = 0.21),
q_x = q_zz = q (open chain, BOTH weak-measurement layers active).

Nothing is duplicated or modified:
  * dynamics   : `lt_evolve_one_trial` from `renyi2_learning_to_trivial_lambda0.7_core.jl`
                 (the lambda-generic doubled-MPS Born-sampling evolution; layer order
                 weak X -> X dephasing -> weak ZZ -> ZZ dephasing, initial state
                 |up..up><up..up|, interleaved bra/ket layout, T_max layers);
  * observable : `dual_ea_moments` from `dual_ea_right_boundary_core.jl`
                 (exact dual-EA moments M2, M4 and B = 1 - M4/(3 M2^2) of one trajectory
                 density matrix, ordinary trace expectations E(O) = Tr(rho O)/Tr(rho)
                 squared, dual endpoints a = 2..L, N = L-1; see that file's docstring).

The observable is evaluated on the state returned by the evolution, i.e. after the last
ZZ-dephasing step of layer t = T_max (same measurement time as the right-boundary
dual-EA scan).

Per-trajectory seeds are drawn exactly as in `lt_run_fixed_lambda_point`
(`Int(rand(MersenneTwister(seed), UInt32))` once per trajectory), so for the same `seed`
this scan analyses the SAME Born trajectories as the Renyi-2 learning-to-trivial scan at
the same (lambda_x, lambda_zz).

Exact limit used as a test: at q = 1/2 the final ZZ-dephasing layer removes every
non-empty X-string on endpoints 2:L, so for every trajectory M2 = 1/N,
M4 = (3N^2 - 2N)/N^4, B = 2/(3N). Unlike the right boundary, q = 0 is NOT an X-eigenstate
product state here (the ZZ measurements compete with the X measurements), so B(q=0) is
not assumed to be 2/3.
"""

using Random
using Statistics
using LinearAlgebra
using ITensors
using ITensorMPS

if !isdefined(Main, :lt_evolve_one_trial)
    include(joinpath(@__DIR__, "renyi2_learning_to_trivial_lambda0.7_core.jl"))
end
if !isdefined(Main, :dual_ea_moments)
    include(joinpath(@__DIR__, "dual_ea_right_boundary_core.jl"))
end

"""
    dual_ea_tl_run_point(L, q; lambda_x, lambda_zz, ntrials, T_max, maxdim, cutoff, seed, nboot, verbose)

Run `ntrials` Born-sampled trajectories at (L, q) and evaluate the dual-EA moments on each
final density matrix. The two Binder conventions are kept strictly separate:

  B_mean_of_trials  = mean_m B_m
  B_ratio_of_means  = 1 - mean(M4_m) / (3 mean(M2_m)^2)
"""
function dual_ea_tl_run_point(
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

    master_rng = MersenneTwister(seed)   # used ONLY to draw trajectory seeds, as in lt_run_fixed_lambda_point

    trial_seeds = Vector{Int}(undef, ntrials)
    M2s = fill(NaN, ntrials)
    M4s = fill(NaN, ntrials)
    Bs = fill(NaN, ntrials)
    obs_seconds = zeros(ntrials)
    linkdims = zeros(Int, ntrials)
    trace_errors = zeros(ntrials)

    for trial in 1:ntrials
        trial_seeds[trial] = Int(rand(master_rng, UInt32))

        evolved = lt_evolve_one_trial(
            L; lambda_x=lambda_x, lambda_zz=lambda_zz, q=q, T_max=T_max,
            maxdim=maxdim, cutoff=cutoff, seed=trial_seeds[trial],
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
