"""
renyi2_learning_to_trivial_lambda0.9_core.jl

Shared, proposal-aligned physics core for the fixed-lambda=0.9 interior slice
("learning-to-trivial" transition scan). This module is a direct translation
of the validated notebook

    lambda_0p7_renyi2_proposal_aligned_v3.ipynb

into a plain Julia include file, following the same pattern as
`renyi2_learning_to_trivial_lambda0.7_core.jl` and `renyi2_right_boundary_core.jl`,
so that every entry point for this slice (cluster jobs via
`run_learning_to_trivial_lambda0.9_scan.jl`) uses exactly the same, single
implementation.

Exact phase-diagram parameterization:

    lambda_x = delta * lambda,  lambda_zz = delta * (1 - lambda),  q_x = q_zz = q

with delta = 0.7. The fixed interior slice studied here is lambda = 0.9, i.e.

    lambda_x = 0.63,  lambda_zz = 0.07,  q_x = q_zz = q in [0, 1/2].

Layer order per full time step (BOTH weak-measurement layers are active here
since lambda_zz != 0):
  1. weak X measurements, exact conditional Born sampling;
  2. X dephasing channel, strength q_x = q;
  3. weak ZZ measurements, exact conditional Born sampling;
  4. ZZ dephasing channel, strength q_zz = q.

Observable: the proposal Renyi-2 Binder is computed PER Born-sampled
trajectory from the replica-overlap operator

    Q = sum_i Z_i^bra Z_i^ket

(NOT a bra-only magnetization). The reported estimate is the trajectory
average

    B2_bar = mean_m B2(m).
"""

using Random
using Statistics
using LinearAlgebra
using ITensors
using ITensorMPS

# ============================================================
# Pauli matrices / constants
# ============================================================
const LT_sigma_x = Float64[0 1; 1 0]
const LT_sigma_z = Float64[1 0; 0 -1]
const LT_identity_2 = Matrix{Float64}(I, 2, 2)

# ============================================================
# Doubled-MPS utilities (bra/ket pairs: site 2i-1 = bra, 2i = ket)
# ============================================================
function lt_pair_mps(local_vec::AbstractVector, s_bra::Index, s_ket::Index)
    @assert length(local_vec) == 4
    return MPS(collect(local_vec), [s_bra, s_ket])
end

function lt_product_density_mps(sites_in, local_vec::AbstractVector)
    sites = collect(sites_in)
    @assert iseven(length(sites))

    tensors = ITensor[]
    for n in 1:2:length(sites)
        pair = lt_pair_mps(local_vec, sites[n], sites[n + 1])
        push!(tensors, pair[1], pair[2])
    end

    return MPS(tensors)
end

lt_bell_state(sites_in) =
    lt_product_density_mps(collect(sites_in), [1.0, 0.0, 0.0, 1.0])

lt_initial_repetition_code_state(sites_in) =
    lt_product_density_mps(collect(sites_in), [1.0, 0.0, 0.0, 0.0])

lt_doubled_trace(rho::MPS, trace_bra::MPS) = real(inner(trace_bra, rho))
lt_hilbert_schmidt_norm(rho::MPS) = real(inner(rho, rho))

function lt_trace_normalize(rho::MPS, trace_bra::MPS; atol::Float64=1e-13)
    tr_rho = lt_doubled_trace(rho, trace_bra)

    if !isfinite(tr_rho) || tr_rho <= atol
        error("Invalid density-matrix trace: Tr(rho) = $tr_rho")
    end

    return rho / tr_rho
end

function lt_physical_bonds(L::Int)
    @assert L >= 2
    return [(i, i + 1) for i in 1:(L - 1)]
end

function lt_bond_dimension_or_one(rho::MPS, bond::Int)
    link = linkind(rho, bond)
    return isnothing(link) ? 1 : dim(link)
end

function lt_max_interphysical_linkdim(rho::MPS, L::Int)
    L <= 1 && return 1
    return maximum(lt_bond_dimension_or_one(rho, 2i) for i in 1:(L - 1))
end

# ============================================================
# Weak X / ZZ measurements: exact conditional Born sampling
# (analytic single-branch probabilities, as in the notebook)
# ============================================================
function lt_bounded_pauli_expectation(
    rho::MPS,
    trace_bra::MPS,
    operator_gates;
    context::AbstractString,
    physical_tol::Float64=1e-7,
)
    input_trace = lt_doubled_trace(rho, trace_bra)

    if !isfinite(input_trace) || input_trace <= 0
        error("$context: invalid input trace $input_trace.")
    end

    # operator_gates are one-site Pauli gates, so this contraction does not
    # require bond-dimension truncation.
    rho_with_operator = apply(operator_gates, rho)

    expectation = real(inner(trace_bra, rho_with_operator)) / input_trace

    if !isfinite(expectation)
        error("$context: non-finite Pauli expectation $expectation.")
    end

    if abs(expectation) > 1 + physical_tol
        error(
            "$context: unphysical Pauli expectation $expectation. " *
            "Increase maxdim and/or decrease cutoff."
        )
    end

    return clamp(expectation, -1.0, 1.0)
end

function lt_binary_weak_measurement_probabilities(
    expectation::Real,
    lambda::Float64;
    context::AbstractString,
    probability_tol::Float64=1e-12,
)
    @assert 0.0 <= lambda <= 1.0

    coefficient = lambda / (1 + lambda^2)

    probabilities = Float64[
        0.5 + coefficient * expectation,
        0.5 - coefficient * expectation,
    ]

    if any(p -> !isfinite(p) || p < -probability_tol, probabilities)
        error("$context: invalid analytic Born probabilities $(probabilities).")
    end

    probabilities = max.(probabilities, 0.0)
    probabilities ./= sum(probabilities)

    return probabilities
end

function lt_x_measurement_probabilities(
    rho::MPS, sites, i::Int, lambda_x::Float64, trace_bra::MPS,
)
    bra_i = 2i - 1

    expectation = lt_bounded_pauli_expectation(
        rho, trace_bra, [op(LT_sigma_x, sites[bra_i])];
        context="X measurement at site $i",
    )

    probabilities = lt_binary_weak_measurement_probabilities(
        expectation, lambda_x; context="X measurement at site $i",
    )

    return expectation, probabilities
end

function lt_apply_selected_x_measurement(
    rho::MPS, sites, i::Int, lambda_x::Float64, outcome::Int;
    maxdim::Int, cutoff::Float64,
)
    @assert outcome in (0, 1)

    bra_i = 2i - 1
    ket_i = 2i

    K = (LT_identity_2 + (-1)^outcome * lambda_x * LT_sigma_x) / sqrt(2 * (1 + lambda_x^2))

    K_bra = op(K, sites[bra_i])
    K_ket = op(conj(K), sites[ket_i])

    return apply([K_bra, K_ket], rho; cutoff=cutoff, maxdim=maxdim)
end

function lt_apply_weak_x_measurement_site(
    rho::MPS, sites, i::Int, lambda_x::Float64, trace_bra::MPS, rng::AbstractRNG;
    maxdim::Int, cutoff::Float64,
)
    _, probabilities = lt_x_measurement_probabilities(rho, sites, i, lambda_x, trace_bra)

    selected_index = rand(rng) < probabilities[1] ? 1 : 2
    outcome = selected_index - 1

    rho_selected = lt_apply_selected_x_measurement(
        rho, sites, i, lambda_x, outcome; maxdim=maxdim, cutoff=cutoff,
    )

    return lt_trace_normalize(rho_selected, trace_bra), outcome, probabilities
end

function lt_apply_weak_x_layer(
    rho::MPS, sites, L::Int, lambda_x::Float64, trace_bra::MPS, rng::AbstractRNG;
    maxdim::Int, cutoff::Float64,
)
    outcomes = Vector{Int}(undef, L)

    for i in 1:L
        rho, outcomes[i], _ = lt_apply_weak_x_measurement_site(
            rho, sites, i, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )
    end

    return rho, outcomes
end

function lt_zz_measurement_probabilities(
    rho::MPS, sites, i::Int, j::Int, lambda_zz::Float64, trace_bra::MPS,
)
    @assert j == i + 1

    bra_i = 2i - 1
    bra_j = 2j - 1

    expectation = lt_bounded_pauli_expectation(
        rho, trace_bra,
        [op(LT_sigma_z, sites[bra_i]), op(LT_sigma_z, sites[bra_j])];
        context="ZZ measurement on bond ($i,$j)",
    )

    probabilities = lt_binary_weak_measurement_probabilities(
        expectation, lambda_zz; context="ZZ measurement on bond ($i,$j)",
    )

    return expectation, probabilities
end

function lt_apply_selected_zz_measurement(
    rho::MPS, sites, i::Int, j::Int, lambda_zz::Float64, outcome::Int;
    maxdim::Int, cutoff::Float64,
)
    @assert j == i + 1
    @assert outcome in (0, 1)

    bra_i, ket_i = 2i - 1, 2i
    bra_j, ket_j = 2j - 1, 2j

    identity_bra = op(LT_identity_2, sites[bra_i]) * op(LT_identity_2, sites[bra_j])
    zz_bra = op(LT_sigma_z, sites[bra_i]) * op(LT_sigma_z, sites[bra_j])

    identity_ket = op(conj(LT_identity_2), sites[ket_i]) * op(conj(LT_identity_2), sites[ket_j])
    zz_ket = op(conj(LT_sigma_z), sites[ket_i]) * op(conj(LT_sigma_z), sites[ket_j])

    normalization = sqrt(2 * (1 + lambda_zz^2))
    sign = (-1)^outcome

    K_bra = (identity_bra + sign * lambda_zz * zz_bra) / normalization
    K_ket = (identity_ket + sign * lambda_zz * zz_ket) / normalization

    return apply([K_bra, K_ket], rho; cutoff=cutoff, maxdim=maxdim)
end

function lt_apply_weak_zz_measurement_bond(
    rho::MPS, sites, i::Int, j::Int, lambda_zz::Float64, trace_bra::MPS, rng::AbstractRNG;
    maxdim::Int, cutoff::Float64,
)
    _, probabilities = lt_zz_measurement_probabilities(rho, sites, i, j, lambda_zz, trace_bra)

    selected_index = rand(rng) < probabilities[1] ? 1 : 2
    outcome = selected_index - 1

    rho_selected = lt_apply_selected_zz_measurement(
        rho, sites, i, j, lambda_zz, outcome; maxdim=maxdim, cutoff=cutoff,
    )

    return lt_trace_normalize(rho_selected, trace_bra), outcome, probabilities
end

function lt_apply_weak_zz_layer(
    rho::MPS, sites, L::Int, lambda_zz::Float64, trace_bra::MPS, rng::AbstractRNG;
    maxdim::Int, cutoff::Float64,
)
    bonds = lt_physical_bonds(L)
    outcomes = Vector{Int}(undef, length(bonds))

    for (n, (i, j)) in enumerate(bonds)
        rho, outcomes[n], _ = lt_apply_weak_zz_measurement_bond(
            rho, sites, i, j, lambda_zz, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )
    end

    return rho, outcomes
end

# ============================================================
# X / ZZ dephasing superoperator gates
# ============================================================
function lt_build_x_dephasing_gates(sites, L::Int, q_x::Float64)
    @assert 0.0 <= q_x <= 0.5

    gates = Vector{ITensor}(undef, L)

    for i in 1:L
        bra_i = 2i - 1
        ket_i = 2i

        identity_gate = op(LT_identity_2, sites[bra_i]) * op(conj(LT_identity_2), sites[ket_i])
        flip_gate = op(LT_sigma_x, sites[bra_i]) * op(conj(LT_sigma_x), sites[ket_i])

        gates[i] = (1 - q_x) * identity_gate + q_x * flip_gate
    end

    return gates
end

function lt_build_zz_dephasing_gates(sites, L::Int, q_zz::Float64)
    @assert 0.0 <= q_zz <= 0.5

    bonds = lt_physical_bonds(L)
    gates = Vector{ITensor}(undef, length(bonds))

    for (n, (i, j)) in enumerate(bonds)
        bra_i, ket_i = 2i - 1, 2i
        bra_j, ket_j = 2j - 1, 2j

        identity_gate =
            op(LT_identity_2, sites[bra_i]) * op(LT_identity_2, sites[bra_j]) *
            op(conj(LT_identity_2), sites[ket_i]) * op(conj(LT_identity_2), sites[ket_j])

        zz_gate =
            op(LT_sigma_z, sites[bra_i]) * op(LT_sigma_z, sites[bra_j]) *
            op(conj(LT_sigma_z), sites[ket_i]) * op(conj(LT_sigma_z), sites[ket_j])

        gates[n] = (1 - q_zz) * identity_gate + q_zz * zz_gate
    end

    return gates
end

function lt_apply_channel_layer(
    rho::MPS, gates::Vector{ITensor}, trace_bra::MPS; maxdim::Int, cutoff::Float64,
)
    for gate in gates
        rho = apply([gate], rho; cutoff=cutoff, maxdim=maxdim)
    end

    return lt_trace_normalize(rho, trace_bra)
end

# ============================================================
# One Born-sampled trajectory at the fixed lambda_x/lambda_zz slice
# ============================================================
"""
    lt_evolve_one_trial(L; lambda_x, lambda_zz, q, T_max, maxdim, cutoff, seed)

Evolve one Born-sampled trajectory with both weak-measurement layers active:
weak X -> X dephasing -> weak ZZ -> ZZ dephasing, repeated `T_max` times.
"""
function lt_evolve_one_trial(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    q::Float64,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    seed::Int,
)
    @assert L >= 2
    @assert 0.0 <= q <= 0.5
    @assert T_max >= 0

    rng = MersenneTwister(seed)

    sites = siteinds("Qubit", 2L)
    trace_bra = lt_bell_state(sites)
    rho = lt_trace_normalize(lt_initial_repetition_code_state(sites), trace_bra)

    x_dephasing_gates = lt_build_x_dephasing_gates(sites, L, q)
    zz_dephasing_gates = lt_build_zz_dephasing_gates(sites, L, q)

    max_trace_error = abs(lt_doubled_trace(rho, trace_bra) - 1.0)
    max_cross_site_bond = lt_max_interphysical_linkdim(rho, L)

    for _ in 1:T_max
        # 1. Weak X measurements.
        rho, _ = lt_apply_weak_x_layer(
            rho, sites, L, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )

        # 2. X dephasing, q_x = q.
        rho = lt_apply_channel_layer(
            rho, x_dephasing_gates, trace_bra; maxdim=maxdim, cutoff=cutoff,
        )

        # 3. Weak ZZ measurements.
        rho, _ = lt_apply_weak_zz_layer(
            rho, sites, L, lambda_zz, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )

        # 4. ZZ dephasing, q_zz = q.
        rho = lt_apply_channel_layer(
            rho, zz_dephasing_gates, trace_bra; maxdim=maxdim, cutoff=cutoff,
        )

        max_trace_error = max(
            max_trace_error, abs(lt_doubled_trace(rho, trace_bra) - 1.0),
        )
        max_cross_site_bond = max(
            max_cross_site_bond, lt_max_interphysical_linkdim(rho, L),
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

# ============================================================
# Proposal Renyi-2 moments and Binder (replica-overlap operator)
# ============================================================
function lt_build_replica_overlap_mpo(sites, L::Int)
    os = OpSum()
    for i in 1:L
        os += 1.0, "Z", 2i - 1, "Z", 2i
    end
    return MPO(os, sites)
end

function lt_binder_from_moments(M2::Real, M4::Real; tol::Float64=1e-13)
    if !isfinite(M2) || !isfinite(M4) || M2 <= tol
        return NaN
    end
    return 1.0 - M4 / (3.0 * M2^2)
end

"""
    lt_renyi2_binder_one_trajectory(rho, sites, L; maxdim, cutoff)

Compute the replica-overlap Renyi-2 moments/Binder for a single Born-sampled
trajectory's doubled-MPS state `rho`.
"""
function lt_renyi2_binder_one_trajectory(
    rho::MPS, sites, L::Int; maxdim::Int, cutoff::Float64, norm_tol::Float64=1e-14,
)
    Q = lt_build_replica_overlap_mpo(sites, L)

    hs_norm = lt_hilbert_schmidt_norm(rho)
    if !isfinite(hs_norm) || hs_norm <= norm_tol
        error("Invalid doubled-state norm: <rho|rho> = $hs_norm")
    end

    psi1 = apply(Q, rho; cutoff=cutoff, maxdim=maxdim)
    psi2 = apply(Q, psi1; cutoff=cutoff, maxdim=maxdim)

    numerator_2 = real(inner(psi1, psi1))
    numerator_4 = real(inner(psi2, psi2))

    M2 = numerator_2 / (L^2 * hs_norm)
    M4 = numerator_4 / (L^4 * hs_norm)
    B2 = lt_binder_from_moments(M2, M4)

    return (M2=M2, M4=M4, B2=B2, purity=hs_norm)
end

# ============================================================
# Bootstrap of the trajectory-averaged proposal observable
# ============================================================
function lt_bootstrap_mean(
    values::AbstractVector; nboot::Int=1000, rng::AbstractRNG=MersenneTwister(1234),
)
    finite_values = collect(filter(isfinite, values))
    n = length(finite_values)

    if n < 2
        return (standard_error=NaN, ci_low=NaN, ci_high=NaN)
    end

    samples = Vector{Float64}(undef, nboot)
    for b in 1:nboot
        indices = rand(rng, 1:n, n)
        samples[b] = mean(finite_values[indices])
    end

    return (
        standard_error=std(samples; corrected=true),
        ci_low=quantile(samples, 0.025),
        ci_high=quantile(samples, 0.975),
    )
end

# ============================================================
# One parameter point (L, q): ntrials Born-sampled trajectories
# ============================================================
"""
    lt_run_fixed_lambda_point(L, q; lambda_x, lambda_zz, ntrials, T_max,
                               maxdim, cutoff, obs_maxdim, obs_cutoff, seed, nboot)

Run `ntrials` independent Born-sampled trajectories at fixed (L, q) with both
weak-measurement layers active, and return the proposal-aligned Renyi-2
Binder summary (trajectory-averaged B2, diagnostics, bootstrap uncertainty).
"""
function lt_run_fixed_lambda_point(
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

    master_rng = MersenneTwister(seed)

    M2s = Vector{Float64}(undef, ntrials)
    M4s = Vector{Float64}(undef, ntrials)
    B2s = Vector{Float64}(undef, ntrials)
    purities = Vector{Float64}(undef, ntrials)
    cross_site_dims = Vector{Int}(undef, ntrials)
    trace_errors = Vector{Float64}(undef, ntrials)

    for trial in 1:ntrials
        trial_seed = Int(rand(master_rng, UInt32))

        evolved = lt_evolve_one_trial(
            L; lambda_x=lambda_x, lambda_zz=lambda_zz, q=q, T_max=T_max,
            maxdim=maxdim, cutoff=cutoff, seed=trial_seed,
        )

        observable = lt_renyi2_binder_one_trajectory(
            evolved.rho, evolved.sites, L; maxdim=obs_maxdim, cutoff=obs_cutoff,
        )

        M2s[trial] = observable.M2
        M4s[trial] = observable.M4
        B2s[trial] = observable.B2
        purities[trial] = observable.purity
        cross_site_dims[trial] = evolved.max_interphysical_linkdim
        trace_errors[trial] = evolved.max_trace_error
    end

    valid = [
        isfinite(M2s[i]) && isfinite(M4s[i]) && isfinite(B2s[i]) &&
        isfinite(purities[i]) && M2s[i] > 0 && purities[i] > 0
        for i in eachindex(B2s)
    ]

    n_valid = count(valid)
    n_invalid = ntrials - n_valid

    if n_valid < 2
        return (
            B=NaN, B_mean_of_trials=NaN, B_std_of_trials=NaN,
            M2_bar=NaN, M4_bar=NaN, purity_bar=NaN,
            B_bootstrap_se=NaN, B_ci_low=NaN, B_ci_high=NaN,
            B2_ratio_of_mean_moments=NaN,
            ntrials=ntrials, n_valid=n_valid, n_invalid=n_invalid,
            max_interphysical_linkdim=maximum(cross_site_dims),
            max_trace_error=maximum(trace_errors),
        )
    end

    valid_M2 = M2s[valid]
    valid_M4 = M4s[valid]
    valid_B2 = B2s[valid]
    valid_purity = purities[valid]

    # Proposal observable: the trajectory-averaged Binder.
    B2_mean = mean(valid_B2)

    boot = lt_bootstrap_mean(valid_B2; nboot=nboot, rng=MersenneTwister(seed + 918273))

    # Diagnostic only; not the proposal average.
    M2_mean = mean(valid_M2)
    M4_mean = mean(valid_M4)
    B2_ratio_of_mean_moments = lt_binder_from_moments(M2_mean, M4_mean)

    return (
        B=B2_mean,
        B_mean_of_trials=B2_mean,
        B_std_of_trials=std(valid_B2; corrected=true),
        M2_bar=M2_mean,
        M4_bar=M4_mean,
        purity_bar=mean(valid_purity),
        B_bootstrap_se=boot.standard_error,
        B_ci_low=boot.ci_low,
        B_ci_high=boot.ci_high,
        B2_ratio_of_mean_moments=B2_ratio_of_mean_moments,
        ntrials=ntrials,
        n_valid=n_valid,
        n_invalid=n_invalid,
        max_interphysical_linkdim=maximum(cross_site_dims),
        max_trace_error=maximum(trace_errors),
    )
end
