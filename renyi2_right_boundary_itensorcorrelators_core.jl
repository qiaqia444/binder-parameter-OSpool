"""
renyi2_right_boundary_core.jl

Shared, proposal-aligned physics core for the right-boundary (measurement-induced
transition) scan, using ITensorCorrelators.jl for the Rényi-2 moments. This module is a direct translation of the validated notebook

    right_boundary_renyi2_proposal_aligned_v4.ipynb

into a plain Julia include file so that every right-boundary entry point
(cluster jobs via `run_right_boundary_scan.jl`, local runs via
`run_right_boundary_simple.jl`) uses exactly the same, single implementation.

Exact phase-diagram parameterization (see notebook for derivation):

    lambda_x = delta * lambda,  lambda_zz = delta * (1 - lambda),  q_x = q_zz = q

with delta = 0.7. The right edge is lambda = 1, i.e.

    lambda_x = 0.7,  lambda_zz = 0,  q_x = q_zz = q in [0, 1/2].

Layer order per full time step:
  1. weak X measurements, exact Born sampling over both outcome branches;
  2. X dephasing channel, strength q_x = q;
  3. weak ZZ measurements -- absent here because lambda_zz = 0;
  4. ZZ dephasing channel, strength q_zz = q (still present even though the
     weak ZZ measurement layer is absent).

Observable: the proposal Renyi-2 Binder is computed PER Born-sampled
trajectory from the replica-overlap operator

    Q = sum_i Z_i^bra Z_i^ket

(NOT a bra-only magnetization -- that would be a different, non-proposal
observable). The reported estimate is the trajectory average

    B2_bar = mean_m B2(m).

This file intentionally implements only the lambda_zz = 0 (right-edge) slice:
there is no weak ZZ measurement code path, matching the notebook exactly.
"""

# Install ITensorCorrelators.jl once in the active Julia environment:
#
#   using Pkg
#   Pkg.add(url="https://github.com/ITensor/ITensorCorrelators.jl.git")
#
using Random
using Statistics
using LinearAlgebra
using ITensors
using ITensorMPS
using ITensorCorrelators

# ============================================================
# Pauli matrices / constants
# ============================================================
const RB_sigma_x = Float64[0 1; 1 0]
const RB_sigma_z = Float64[1 0; 0 -1]
const RB_identity_2 = Matrix{Float64}(I, 2, 2)

# ============================================================
# Doubled-MPS utilities (bra/ket pairs: site 2i-1 = bra, 2i = ket)
# ============================================================
function rb_pair_mps(local_vec::AbstractVector, s_bra::Index, s_ket::Index)
    @assert length(local_vec) == 4
    return MPS(collect(local_vec), [s_bra, s_ket])
end

function rb_product_density_mps(sites_in, local_vec::AbstractVector)
    sites = collect(sites_in)
    @assert iseven(length(sites))

    tensors = ITensor[]
    for n in 1:2:length(sites)
        pair = rb_pair_mps(local_vec, sites[n], sites[n + 1])
        push!(tensors, pair[1], pair[2])
    end

    return MPS(tensors)
end

rb_bell_state(sites_in) =
    rb_product_density_mps(collect(sites_in), [1.0, 0.0, 0.0, 1.0])

rb_initial_repetition_code_state(sites_in) =
    rb_product_density_mps(collect(sites_in), [1.0, 0.0, 0.0, 0.0])

rb_doubled_trace(rho::MPS, trace_bra::MPS) = real(inner(trace_bra, rho))
rb_hilbert_schmidt_norm(rho::MPS) = real(inner(rho, rho))

function rb_trace_normalize(rho::MPS, trace_bra::MPS; atol::Float64=1e-13)
    tr_rho = rb_doubled_trace(rho, trace_bra)

    if !isfinite(tr_rho) || tr_rho <= atol
        error("Invalid density-matrix trace: Tr(rho) = $tr_rho")
    end

    return rho / tr_rho
end

function rb_physical_bonds(L::Int)
    @assert L >= 2
    return [(i, i + 1) for i in 1:(L - 1)]
end

function rb_bond_dimension_or_one(rho::MPS, bond::Int)
    link = linkind(rho, bond)
    return isnothing(link) ? 1 : dim(link)
end

function rb_max_interphysical_linkdim(rho::MPS, L::Int)
    L <= 1 && return 1
    return maximum(rb_bond_dimension_or_one(rho, 2i) for i in 1:(L - 1))
end

# ============================================================
# Weak X measurement: exact Born sampling on both branches
# ============================================================
function rb_x_measurement_candidates(
    rho::MPS,
    sites,
    i::Int,
    lambda_x::Float64,
    trace_bra::MPS;
    maxdim::Int,
    cutoff::Float64,
    negative_weight_tol::Float64=1e-10,
    completeness_warn_tol::Float64=5e-8,
    completeness_fail_tol::Float64=1e-6,
)
    @assert 0.0 <= lambda_x <= 1.0

    bra_i = 2i - 1
    ket_i = 2i

    states = Vector{MPS}(undef, 2)
    weights = zeros(Float64, 2)

    for outcome in 0:1
        K = (
            RB_identity_2 + (-1)^outcome * lambda_x * RB_sigma_x
        ) / sqrt(2 * (1 + lambda_x^2))

        K_bra = op(K, sites[bra_i])
        K_ket = op(conj(K), sites[ket_i])

        rho_m = apply([K_bra, K_ket], rho; cutoff=cutoff, maxdim=maxdim)
        weight = rb_doubled_trace(rho_m, trace_bra)

        if !isfinite(weight) || weight < -negative_weight_tol
            error(
                "Invalid Born weight at site $i, outcome $outcome: $weight. " *
                "Increase maxdim or reduce cutoff."
            )
        end

        states[outcome + 1] = rho_m
        weights[outcome + 1] = max(weight, 0.0)
    end

    total_weight = sum(weights)
    input_trace = rb_doubled_trace(rho, trace_bra)

    if total_weight <= negative_weight_tol
        error("Both X-measurement outcomes have zero numerical weight.")
    end

    completeness_error = abs(total_weight - input_trace)
    relative_error = completeness_error / max(abs(input_trace), eps(Float64))

    if relative_error > completeness_fail_tol
        error(
            "Kraus completeness drift is too large: " *
            "sum(weights)=$total_weight, Tr(rho)=$input_trace, " *
            "relative_error=$relative_error. Increase maxdim and/or " *
            "decrease cutoff."
        )
    elseif relative_error > completeness_warn_tol
        @warn(
            "Small Kraus completeness drift from MPS truncation",
            site=i,
            relative_error=relative_error,
            maxdim=maxdim,
            cutoff=cutoff,
        )
    end

    return states, weights ./ total_weight
end

function rb_apply_weak_x_measurement_site(
    rho::MPS,
    sites,
    i::Int,
    lambda_x::Float64,
    trace_bra::MPS,
    rng::AbstractRNG;
    maxdim::Int,
    cutoff::Float64,
)
    states, probabilities = rb_x_measurement_candidates(
        rho, sites, i, lambda_x, trace_bra; maxdim=maxdim, cutoff=cutoff,
    )

    selected = rand(rng) < probabilities[1] ? 1 : 2
    rho_selected = states[selected]

    return rb_trace_normalize(rho_selected, trace_bra), selected - 1, probabilities
end

function rb_apply_weak_x_layer(
    rho::MPS,
    sites,
    L::Int,
    lambda_x::Float64,
    trace_bra::MPS,
    rng::AbstractRNG;
    maxdim::Int,
    cutoff::Float64,
)
    outcomes = Vector{Int}(undef, L)

    for i in 1:L
        rho, outcomes[i], _ = rb_apply_weak_x_measurement_site(
            rho, sites, i, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )
    end

    return rho, outcomes
end

# ============================================================
# X / ZZ dephasing superoperator gates
# ============================================================
function rb_build_x_dephasing_gates(sites, L::Int, q_x::Float64)
    @assert 0.0 <= q_x <= 0.5

    gates = Vector{ITensor}(undef, L)

    for i in 1:L
        bra_i = 2i - 1
        ket_i = 2i

        identity_gate =
            op(RB_identity_2, sites[bra_i]) * op(conj(RB_identity_2), sites[ket_i])

        flip_gate =
            op(RB_sigma_x, sites[bra_i]) * op(conj(RB_sigma_x), sites[ket_i])

        gates[i] = (1 - q_x) * identity_gate + q_x * flip_gate
    end

    return gates
end

function rb_build_zz_dephasing_gates(sites, L::Int, q_zz::Float64)
    @assert 0.0 <= q_zz <= 0.5

    bonds = rb_physical_bonds(L)
    gates = Vector{ITensor}(undef, length(bonds))

    for (n, (i, j)) in enumerate(bonds)
        bra_i, ket_i = 2i - 1, 2i
        bra_j, ket_j = 2j - 1, 2j

        identity_gate =
            op(RB_identity_2, sites[bra_i]) * op(RB_identity_2, sites[bra_j]) *
            op(conj(RB_identity_2), sites[ket_i]) * op(conj(RB_identity_2), sites[ket_j])

        zz_gate =
            op(RB_sigma_z, sites[bra_i]) * op(RB_sigma_z, sites[bra_j]) *
            op(conj(RB_sigma_z), sites[ket_i]) * op(conj(RB_sigma_z), sites[ket_j])

        gates[n] = (1 - q_zz) * identity_gate + q_zz * zz_gate
    end

    return gates
end

function rb_apply_channel_layer(
    rho::MPS, gates::Vector{ITensor}, trace_bra::MPS; maxdim::Int, cutoff::Float64,
)
    for gate in gates
        rho = apply([gate], rho; cutoff=cutoff, maxdim=maxdim)
    end

    return rb_trace_normalize(rho, trace_bra)
end

# ============================================================
# One Born-sampled right-edge trajectory
# ============================================================
"""
    rb_evolve_right_edge_one_trial(L; lambda_x, q, T_max, maxdim, cutoff, seed)

Evolve one Born-sampled trajectory on the right edge (lambda_zz = 0):
weak X measurement -> X dephasing -> ZZ dephasing, repeated `T_max` times.
"""
function rb_evolve_right_edge_one_trial(
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
    @assert T_max >= 0

    rng = MersenneTwister(seed)

    sites = siteinds("Qubit", 2L)
    trace_bra = rb_bell_state(sites)
    rho = rb_trace_normalize(rb_initial_repetition_code_state(sites), trace_bra)

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

# ============================================================
# Proposal Renyi-2 moments and Binder (replica-overlap operator)
# ============================================================
function rb_build_replica_overlap_mpo(sites, L::Int)
    os = OpSum()
    for i in 1:L
        os += 1.0, "Z", 2i - 1, "Z", 2i
    end
    return MPO(os, sites)
end

function rb_binder_from_moments(M2::Real, M4::Real; tol::Float64=1e-13)
    if !isfinite(M2) || !isfinite(M4) || M2 <= tol
        return NaN
    end
    return 1.0 - M4 / (3.0 * M2^2)
end

"""
    rb_renyi2_binder_one_trajectory(rho, sites, L; maxdim, cutoff)

Compute the replica-overlap Renyi-2 moments/Binder for a single Born-sampled
trajectory's doubled-MPS state `rho`.
"""
function rb_renyi2_binder_one_trajectory(
    rho::MPS,
    sites,
    L::Int;
    maxdim::Int,
    cutoff::Float64,
    norm_tol::Float64=1e-14,
)
    Q = rb_build_replica_overlap_mpo(sites, L)

    hs_norm = rb_hilbert_schmidt_norm(rho)
    if !isfinite(hs_norm) || hs_norm <= norm_tol
        error("Invalid doubled-state norm: <rho|rho> = $hs_norm")
    end

    psi1 = apply(Q, rho; cutoff=cutoff, maxdim=maxdim)
    psi2 = apply(Q, psi1; cutoff=cutoff, maxdim=maxdim)

    numerator_2 = real(inner(psi1, psi1))
    numerator_4 = real(inner(psi2, psi2))

    M2 = numerator_2 / (L^2 * hs_norm)
    M4 = numerator_4 / (L^4 * hs_norm)
    B2 = rb_binder_from_moments(M2, M4)

    return (M2=M2, M4=M4, B2=B2, purity=hs_norm)
end


# ============================================================
# ITensorCorrelators.jl implementation of the same Rényi-2 Binder
# ============================================================

"""
    rb_sum_correlators_batched(psi, operators, site_tuples; batch_size)

Evaluate and sum a list of n-point correlators in batches. The official
ITensorCorrelators.jl interface is

    correlator(psi, ("Op1", ..., "OpN"), [(i1, ..., iN), ...])

and returns a dictionary keyed by site tuples.
"""
function rb_sum_correlators_batched(
    psi::MPS,
    operators::Tuple,
    site_tuples::AbstractVector;
    batch_size::Int=2048,
)
    isempty(site_tuples) && return 0.0
    @assert batch_size >= 1

    total = 0.0

    for first_index in 1:batch_size:length(site_tuples)
        last_index = min(first_index + batch_size - 1, length(site_tuples))
        batch = site_tuples[first_index:last_index]

        values_dict = correlator(psi, operators, batch)

        for value in values(values_dict)
            total += real(value)
        end
    end

    return total
end


"""
    rb_renyi2_binder_one_trajectory_correlator(
        rho, sites, L;
        pair_batch_size=4096,
        quartet_batch_size=1024,
        norm_tol=1e-14,
    )

Compute exactly the same proposal Rényi-2 moments and Binder as
`rb_renyi2_binder_one_trajectory`, but use ITensorCorrelators.jl instead of
constructing and applying the global Q MPO.

The doubled-MPS ordering is

    physical site i -> bra site 2i-1, ket site 2i,

and the local replica-overlap variable is

    r_i = Z_(2i-1) Z_(2i).

For commuting r_i with r_i^2 = 1,

    Q = sum_i r_i,
    <Q^2> = L + 2 sum_{i<j} <r_i r_j>,

and

    <Q^4> = 3L^2 - 2L
             + 4(3L-4) sum_{i<j} <r_i r_j>
             + 24 sum_{i<j<k<l} <r_i r_j r_k r_l>.
"""
function rb_renyi2_binder_one_trajectory_correlator(
    rho::MPS,
    sites,
    L::Int;
    pair_batch_size::Int=4096,
    quartet_batch_size::Int=1024,
    norm_tol::Float64=1e-14,
)
    @assert L >= 1
    @assert length(sites) == 2L

    hs_norm = rb_hilbert_schmidt_norm(rho)
    if !isfinite(hs_norm) || hs_norm <= norm_tol
        error("Invalid doubled-state norm: <rho|rho> = $hs_norm")
    end

    # Normalize so correlator() directly returns
    # <rho|O|rho>/<rho|rho>.
    psi = deepcopy(rho / sqrt(hs_norm))

    # --------------------------------------------------------
    # Pair sum: sum_{i<j} <r_i r_j>
    # Each r_i contributes Z on bra and Z on ket, so this is
    # a four-point correlator on the 2L-site doubled MPS.
    # --------------------------------------------------------
    pair_sites = NTuple{4,Int}[]

    if L >= 2
        sizehint!(pair_sites, L * (L - 1) ÷ 2)

        for i in 1:(L - 1)
            for j in (i + 1):L
                push!(pair_sites, (2i - 1, 2i, 2j - 1, 2j))
            end
        end
    end

    pair_sum = rb_sum_correlators_batched(
        psi,
        ("Z", "Z", "Z", "Z"),
        pair_sites;
        batch_size=pair_batch_size,
    )

    Q2_expectation = Float64(L) + 2.0 * pair_sum

    # --------------------------------------------------------
    # Quartet sum: sum_{i<j<k<l} <r_i r_j r_k r_l>
    # This is an eight-point correlator on the doubled MPS.
    # --------------------------------------------------------
    quartet_sites = NTuple{8,Int}[]

    if L >= 4
        nquartets = L * (L - 1) * (L - 2) * (L - 3) ÷ 24
        sizehint!(quartet_sites, nquartets)

        for i in 1:(L - 3)
            for j in (i + 1):(L - 2)
                for k in (j + 1):(L - 1)
                    for l in (k + 1):L
                        push!(
                            quartet_sites,
                            (
                                2i - 1, 2i,
                                2j - 1, 2j,
                                2k - 1, 2k,
                                2l - 1, 2l,
                            ),
                        )
                    end
                end
            end
        end
    end

    quartet_sum = rb_sum_correlators_batched(
        psi,
        ("Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z"),
        quartet_sites;
        batch_size=quartet_batch_size,
    )

    Q4_expectation =
        3.0 * L^2 -
        2.0 * L +
        4.0 * (3.0 * L - 4.0) * pair_sum +
        24.0 * quartet_sum

    M2 = Q2_expectation / L^2
    M4 = Q4_expectation / L^4
    B2 = rb_binder_from_moments(M2, M4)

    return (
        M2=M2,
        M4=M4,
        B2=B2,
        purity=hs_norm,
        Q2_expectation=Q2_expectation,
        Q4_expectation=Q4_expectation,
        pair_sum=pair_sum,
        quartet_sum=quartet_sum,
        n_pair_correlators=length(pair_sites),
        n_quartet_correlators=length(quartet_sites),
    )
end


"""
    rb_compare_mpo_and_correlator(rho, sites, L; ...)

Independent validation helper. It compares the original global-MPO method with
ITensorCorrelators.jl on the same final doubled MPS.
"""
function rb_compare_mpo_and_correlator(
    rho::MPS,
    sites,
    L::Int;
    mpo_maxdim::Int=2048,
    mpo_cutoff::Float64=1e-14,
    pair_batch_size::Int=4096,
    quartet_batch_size::Int=1024,
)
    mpo_result = rb_renyi2_binder_one_trajectory(
        rho,
        sites,
        L;
        maxdim=mpo_maxdim,
        cutoff=mpo_cutoff,
    )

    correlator_result = rb_renyi2_binder_one_trajectory_correlator(
        rho,
        sites,
        L;
        pair_batch_size=pair_batch_size,
        quartet_batch_size=quartet_batch_size,
    )

    return (
        mpo=mpo_result,
        correlator=correlator_result,
        delta_M2=correlator_result.M2 - mpo_result.M2,
        delta_M4=correlator_result.M4 - mpo_result.M4,
        delta_B2=correlator_result.B2 - mpo_result.B2,
    )
end

# ============================================================
# Bootstrap of the trajectory-averaged proposal observable
# ============================================================
function rb_bootstrap_mean(
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
# One right-edge parameter point (L, q): ntrials Born-sampled trajectories
# ============================================================
"""
    rb_run_right_edge_point(L, q; lambda_x, lambda_zz, ntrials, T_max, maxdim,
                             cutoff, obs_maxdim, obs_cutoff, seed, nboot)

Run `ntrials` independent Born-sampled trajectories at fixed (L, q) on the
right edge (lambda_zz = 0) and return the proposal-aligned Renyi-2 Binder
summary (trajectory-averaged B2, diagnostics, and bootstrap uncertainty).
"""
function rb_run_right_edge_point(
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

    M2s = Vector{Float64}(undef, ntrials)
    M4s = Vector{Float64}(undef, ntrials)
    B2s = Vector{Float64}(undef, ntrials)
    purities = Vector{Float64}(undef, ntrials)
    cross_site_dims = Vector{Int}(undef, ntrials)
    trace_errors = Vector{Float64}(undef, ntrials)

    for trial in 1:ntrials
        trial_seed = Int(rand(master_rng, UInt32))

        evolved = rb_evolve_right_edge_one_trial(
            L; lambda_x=lambda_x, q=q, T_max=T_max, maxdim=maxdim,
            cutoff=cutoff, seed=trial_seed,
        )

        observable = rb_renyi2_binder_one_trajectory_correlator(
            evolved.rho,
            evolved.sites,
            L;
            pair_batch_size=4096,
            quartet_batch_size=1024,
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

    boot = rb_bootstrap_mean(
        valid_B2; nboot=nboot, rng=MersenneTwister(seed + 918273),
    )

    # Diagnostic only; not the proposal average.
    M2_mean = mean(valid_M2)
    M4_mean = mean(valid_M4)
    B2_ratio_of_mean_moments = rb_binder_from_moments(M2_mean, M4_mean)

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
