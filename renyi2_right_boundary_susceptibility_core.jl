"""
renyi2_right_boundary_susceptibility_core.jl

Shared, proposal-aligned physics core for the right-boundary (measurement-induced
transition) scan, using ITensorCorrelators.jl for the Rényi-2 overlap susceptibility. This module is a direct translation of the validated notebook

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

Observable: the Rényi-2 overlap susceptibility is computed PER Born-sampled
trajectory from the replica-overlap operator

    Q = sum_i Z_i^bra Z_i^ket

(NOT a bra-only magnetization). For each trajectory,

    chi2(m) = <Q^2>_m / L = L * M2(m),

and the reported estimate is the trajectory average

    chi2_bar = mean_m chi2(m).

This is the replica-overlap / Edwards-Anderson susceptibility. A connected
version, [<Q^2> - <Q>^2]/L, is also returned as a diagnostic.

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
# Rényi-2 overlap susceptibility (replica-overlap operator)
# ============================================================
function rb_build_replica_overlap_mpo(sites, L::Int)
    os = OpSum()
    for i in 1:L
        os += 1.0, "Z", 2i - 1, "Z", 2i
    end
    return MPO(os, sites)
end

"""
    rb_renyi2_susceptibility_one_trajectory(rho, sites, L; maxdim, cutoff)

Compute the Rényi-2 replica-overlap susceptibility using the global MPO

    Q = sum_i Z_(2i-1) Z_(2i).

The primary observable is

    chi2 = <Q^2> / L = L * M2,

where all expectation values are normalized by <rho|rho>. The connected
susceptibility is returned as a diagnostic:

    chi2_connected = (<Q^2> - <Q>^2) / L.
"""
function rb_renyi2_susceptibility_one_trajectory(
    rho::MPS,
    sites,
    L::Int;
    maxdim::Int,
    cutoff::Float64,
    norm_tol::Float64=1e-14,
)
    @assert L >= 1

    Q = rb_build_replica_overlap_mpo(sites, L)

    hs_norm = rb_hilbert_schmidt_norm(rho)
    if !isfinite(hs_norm) || hs_norm <= norm_tol
        error("Invalid doubled-state norm: <rho|rho> = $hs_norm")
    end

    Qrho = apply(Q, rho; cutoff=cutoff, maxdim=maxdim)

    Q_expectation = real(inner(rho, Qrho)) / hs_norm
    Q2_expectation = real(inner(Qrho, Qrho)) / hs_norm

    M2 = Q2_expectation / L^2
    chi2 = Q2_expectation / L
    chi2_connected = (Q2_expectation - Q_expectation^2) / L
    overlap_density = Q_expectation / L

    return (
        chi=chi2,
        chi2=chi2,
        chi_connected=chi2_connected,
        M2=M2,
        overlap_density=overlap_density,
        Q_expectation=Q_expectation,
        Q2_expectation=Q2_expectation,
        purity=hs_norm,
    )
end

# ============================================================
# ITensorCorrelators.jl implementation
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
    batch_size::Int=4096,
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
    rb_renyi2_susceptibility_one_trajectory_correlator(
        rho, sites, L; single_batch_size=4096, pair_batch_size=4096
    )

Compute the same Rényi-2 overlap susceptibility with ITensorCorrelators.jl.
For

    r_i = Z_(2i-1) Z_(2i),
    Q   = sum_i r_i,

and r_i^2 = 1,

    <Q>   = sum_i <r_i>,
    <Q^2> = L + 2 sum_{i<j} <r_i r_j>.

Thus only two-point and four-point correlators on the doubled MPS are needed;
no eight-point quartet correlators are required.
"""
function rb_renyi2_susceptibility_one_trajectory_correlator(
    rho::MPS,
    sites,
    L::Int;
    single_batch_size::Int=4096,
    pair_batch_size::Int=4096,
    norm_tol::Float64=1e-14,
)
    @assert L >= 1
    @assert length(sites) == 2L

    hs_norm = rb_hilbert_schmidt_norm(rho)
    if !isfinite(hs_norm) || hs_norm <= norm_tol
        error("Invalid doubled-state norm: <rho|rho> = $hs_norm")
    end

    # Normalize so correlator() directly returns <rho|O|rho>/<rho|rho>.
    psi = deepcopy(rho / sqrt(hs_norm))

    # <Q> = sum_i <r_i>; each r_i is a two-point correlator on the doubled MPS.
    single_sites = NTuple{2,Int}[]
    sizehint!(single_sites, L)
    for i in 1:L
        push!(single_sites, (2i - 1, 2i))
    end

    Q_expectation = rb_sum_correlators_batched(
        psi,
        ("Z", "Z"),
        single_sites;
        batch_size=single_batch_size,
    )

    # sum_{i<j} <r_i r_j>; each term is a four-point correlator.
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

    M2 = Q2_expectation / L^2
    chi2 = Q2_expectation / L
    chi2_connected = (Q2_expectation - Q_expectation^2) / L
    overlap_density = Q_expectation / L

    return (
        chi=chi2,
        chi2=chi2,
        chi_connected=chi2_connected,
        M2=M2,
        overlap_density=overlap_density,
        Q_expectation=Q_expectation,
        Q2_expectation=Q2_expectation,
        pair_sum=pair_sum,
        purity=hs_norm,
        n_single_correlators=length(single_sites),
        n_pair_correlators=length(pair_sites),
    )
end

"""
    rb_compare_susceptibility_mpo_and_correlator(rho, sites, L; ...)

Validate the correlator implementation against the independent global-MPO
calculation on the same final doubled MPS.
"""
function rb_compare_susceptibility_mpo_and_correlator(
    rho::MPS,
    sites,
    L::Int;
    mpo_maxdim::Int=2048,
    mpo_cutoff::Float64=1e-14,
    single_batch_size::Int=4096,
    pair_batch_size::Int=4096,
)
    mpo_result = rb_renyi2_susceptibility_one_trajectory(
        rho,
        sites,
        L;
        maxdim=mpo_maxdim,
        cutoff=mpo_cutoff,
    )

    correlator_result = rb_renyi2_susceptibility_one_trajectory_correlator(
        rho,
        sites,
        L;
        single_batch_size=single_batch_size,
        pair_batch_size=pair_batch_size,
    )

    return (
        mpo=mpo_result,
        correlator=correlator_result,
        delta_chi=correlator_result.chi - mpo_result.chi,
        delta_chi_connected=(
            correlator_result.chi_connected - mpo_result.chi_connected
        ),
        delta_M2=correlator_result.M2 - mpo_result.M2,
        delta_Q=correlator_result.Q_expectation - mpo_result.Q_expectation,
        delta_Q2=correlator_result.Q2_expectation - mpo_result.Q2_expectation,
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
    rb_run_right_edge_susceptibility_point(
        L, q; lambda_x, lambda_zz, ntrials, T_max, maxdim, cutoff,
        obs_maxdim, obs_cutoff, seed, nboot
    )

Run `ntrials` independent Born-sampled trajectories at fixed `(L, q)` on the
right edge (`lambda_zz = 0`) and return the trajectory-averaged Rényi-2
overlap susceptibility

    chi2_bar = mean_m [ <Q^2>_m / L ].

`obs_maxdim` and `obs_cutoff` are retained in the interface for compatibility
with the old Binder runner, although the default correlator measurement does
not need them. They can still be used by the MPO validation helper.
"""
function rb_run_right_edge_susceptibility_point(
    L::Int,
    q::Float64;
    lambda_x::Float64,
    lambda_zz::Float64,
    ntrials::Int,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    obs_maxdim::Int=2048,
    obs_cutoff::Float64=1e-14,
    seed::Int,
    nboot::Int=1000,
    single_batch_size::Int=4096,
    pair_batch_size::Int=4096,
)
    @assert ntrials >= 1
    @assert 0.0 <= q <= 0.5
    @assert isapprox(lambda_zz, 0.0; atol=1e-9) (
        "This module only implements the right edge (lambda_zz = 0); " *
        "weak ZZ measurements are not implemented."
    )

    # Keep these arguments visible and type-checked for compatibility.
    @assert obs_maxdim >= 1
    @assert obs_cutoff >= 0.0

    master_rng = MersenneTwister(seed)

    chis = Vector{Float64}(undef, ntrials)
    chi_connecteds = Vector{Float64}(undef, ntrials)
    M2s = Vector{Float64}(undef, ntrials)
    overlap_densities = Vector{Float64}(undef, ntrials)
    purities = Vector{Float64}(undef, ntrials)
    cross_site_dims = Vector{Int}(undef, ntrials)
    trace_errors = Vector{Float64}(undef, ntrials)

    for trial in 1:ntrials
        trial_seed = Int(rand(master_rng, UInt32))

        evolved = rb_evolve_right_edge_one_trial(
            L;
            lambda_x=lambda_x,
            q=q,
            T_max=T_max,
            maxdim=maxdim,
            cutoff=cutoff,
            seed=trial_seed,
        )

        observable = rb_renyi2_susceptibility_one_trajectory_correlator(
            evolved.rho,
            evolved.sites,
            L;
            single_batch_size=single_batch_size,
            pair_batch_size=pair_batch_size,
        )

        chis[trial] = observable.chi
        chi_connecteds[trial] = observable.chi_connected
        M2s[trial] = observable.M2
        overlap_densities[trial] = observable.overlap_density
        purities[trial] = observable.purity
        cross_site_dims[trial] = evolved.max_interphysical_linkdim
        trace_errors[trial] = evolved.max_trace_error
    end

    valid = [
        isfinite(chis[i]) &&
        isfinite(chi_connecteds[i]) &&
        isfinite(M2s[i]) &&
        isfinite(overlap_densities[i]) &&
        isfinite(purities[i]) &&
        chis[i] >= -1e-12 &&
        purities[i] > 0
        for i in eachindex(chis)
    ]

    n_valid = count(valid)
    n_invalid = ntrials - n_valid

    if n_valid < 2
        return (
            chi=NaN,
            susceptibility=NaN,
            chi_mean_of_trials=NaN,
            chi_std_of_trials=NaN,
            chi_connected_bar=NaN,
            M2_bar=NaN,
            overlap_density_bar=NaN,
            purity_bar=NaN,
            chi_bootstrap_se=NaN,
            chi_ci_low=NaN,
            chi_ci_high=NaN,
            ntrials=ntrials,
            n_valid=n_valid,
            n_invalid=n_invalid,
            max_interphysical_linkdim=maximum(cross_site_dims),
            max_trace_error=maximum(trace_errors),
        )
    end

    valid_chi = chis[valid]
    valid_chi_connected = chi_connecteds[valid]
    valid_M2 = M2s[valid]
    valid_overlap_density = overlap_densities[valid]
    valid_purity = purities[valid]

    chi_mean = mean(valid_chi)
    boot = rb_bootstrap_mean(
        valid_chi;
        nboot=nboot,
        rng=MersenneTwister(seed + 918273),
    )

    return (
        chi=chi_mean,
        susceptibility=chi_mean,
        chi_mean_of_trials=chi_mean,
        chi_std_of_trials=std(valid_chi; corrected=true),
        chi_connected_bar=mean(valid_chi_connected),
        M2_bar=mean(valid_M2),
        overlap_density_bar=mean(valid_overlap_density),
        purity_bar=mean(valid_purity),
        chi_bootstrap_se=boot.standard_error,
        chi_ci_low=boot.ci_low,
        chi_ci_high=boot.ci_high,
        ntrials=ntrials,
        n_valid=n_valid,
        n_invalid=n_invalid,
        max_interphysical_linkdim=maximum(cross_site_dims),
        max_trace_error=maximum(trace_errors),
    )
end

# Backward-compatible entry-point name for existing scan scripts.
# The returned NamedTuple now contains susceptibility fields instead of Binder fields.
rb_run_right_edge_point(args...; kwargs...) =
    rb_run_right_edge_susceptibility_point(args...; kwargs...)
