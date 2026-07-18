
"""
lambda_scan_susceptibilities_itensorcorrelators.jl

Scan the paper-like path

    lambda_x  = delta * lambda
    lambda_zz = delta * (1 - lambda)
    q_x = q_zz = q

so increasing `lambda` strengthens the weak X measurement and weakens the
weak ZZ measurement.

The code evolves Born-sampled doubled-MPS density-matrix trajectories and
returns both

    kappa_EA = (1/L) sum_ij <Z_i Z_j>^2

and

    kappa_2  = (1/L) sum_ij
               <rho|r_i r_j|rho>/<rho|rho>,
    r_i = Z_i^bra Z_i^ket.

The Rényi-2 susceptibility uses ITensorCorrelators.jl.
"""

using Random
using Statistics
using LinearAlgebra
using ITensors
using ITensorMPS
using ITensorCorrelators

const LS_X  = Float64[0 1; 1 0]
const LS_Z  = Float64[1 0; 0 -1]
const LS_I2 = Matrix{Float64}(I, 2, 2)

# ---------------------------------------------------------------------------
# Doubled-MPS utilities
# Site ordering: (1_bra, 1_ket, 2_bra, 2_ket, ...)
# ---------------------------------------------------------------------------

function ls_pair_mps(local_vec::AbstractVector, s_bra::Index, s_ket::Index)
    @assert length(local_vec) == 4
    return MPS(collect(local_vec), [s_bra, s_ket])
end

function ls_product_density_mps(sites_in, local_vec::AbstractVector)
    sites = collect(sites_in)
    @assert iseven(length(sites))

    tensors = ITensor[]
    for n in 1:2:length(sites)
        pair = ls_pair_mps(local_vec, sites[n], sites[n + 1])
        push!(tensors, pair[1], pair[2])
    end
    return MPS(tensors)
end

# |I>> on every physical site, used as the trace bra.
ls_trace_state(sites) =
    ls_product_density_mps(sites, [1.0, 0.0, 0.0, 1.0])

# Initial physical state |0...0><0...0|.
ls_initial_state(sites) =
    ls_product_density_mps(sites, [1.0, 0.0, 0.0, 0.0])

ls_trace(rho::MPS, trace_bra::MPS) = real(inner(trace_bra, rho))
ls_hs_norm(rho::MPS) = real(inner(rho, rho))

function ls_trace_normalize(
    rho::MPS,
    trace_bra::MPS;
    atol::Float64=1e-13,
)
    trrho = ls_trace(rho, trace_bra)
    if !isfinite(trrho) || trrho <= atol
        error("Invalid density-matrix trace: Tr(rho) = $trrho")
    end
    return rho / trrho
end

ls_physical_bonds(L::Int) = [(i, i + 1) for i in 1:(L - 1)]

function ls_bond_dimension_or_one(rho::MPS, bond::Int)
    link = linkind(rho, bond)
    return isnothing(link) ? 1 : dim(link)
end

function ls_max_interphysical_linkdim(rho::MPS, L::Int)
    L <= 1 && return 1
    return maximum(ls_bond_dimension_or_one(rho, 2i) for i in 1:(L - 1))
end

# ---------------------------------------------------------------------------
# Generic Born sampling helper
# ---------------------------------------------------------------------------

function ls_sample_candidates(
    states::Vector{MPS},
    weights::Vector{Float64},
    trace_bra::MPS,
    rng::AbstractRNG;
    negative_weight_tol::Float64=1e-10,
    completeness_fail_tol::Float64=1e-6,
)
    @assert length(states) == length(weights) == 2

    for n in eachindex(weights)
        if !isfinite(weights[n]) || weights[n] < -negative_weight_tol
            error("Invalid Born weight: $(weights[n])")
        end
        weights[n] = max(weights[n], 0.0)
    end

    total_weight = sum(weights)
    total_weight > negative_weight_tol ||
        error("Both measurement outcomes have zero numerical weight.")

    probabilities = weights ./ total_weight

    input_trace = sum(weights)
    if abs(sum(probabilities) - 1.0) > completeness_fail_tol
        error("Born probabilities do not sum to one.")
    end

    selected = rand(rng) < probabilities[1] ? 1 : 2
    rho_selected = ls_trace_normalize(states[selected], trace_bra)

    return rho_selected, selected - 1, probabilities, input_trace
end

# ---------------------------------------------------------------------------
# Weak X measurement
#
# K_m^X = [I + (-1)^m lambda_x X] / sqrt(2(1 + lambda_x^2))
# ---------------------------------------------------------------------------

function ls_x_measurement_candidates(
    rho::MPS,
    sites,
    i::Int,
    lambda_x::Float64,
    trace_bra::MPS;
    maxdim::Int,
    cutoff::Float64,
)
    @assert 0.0 <= lambda_x <= 1.0

    bra_i, ket_i = 2i - 1, 2i
    denom = sqrt(2.0 * (1.0 + lambda_x^2))

    states = Vector{MPS}(undef, 2)
    weights = zeros(Float64, 2)

    for outcome in 0:1
        s = (-1.0)^outcome

        K_bra =
            (op(LS_I2, sites[bra_i]) +
             s * lambda_x * op(LS_X, sites[bra_i])) / denom

        K_ket =
            (op(conj(LS_I2), sites[ket_i]) +
             s * lambda_x * op(conj(LS_X), sites[ket_i])) / denom

        supergate = K_bra * K_ket
        rho_m = apply([supergate], rho; cutoff=cutoff, maxdim=maxdim)

        states[outcome + 1] = rho_m
        weights[outcome + 1] = ls_trace(rho_m, trace_bra)
    end

    return states, weights
end

function ls_apply_weak_x_layer(
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
        states, weights = ls_x_measurement_candidates(
            rho, sites, i, lambda_x, trace_bra;
            maxdim=maxdim, cutoff=cutoff,
        )

        rho, outcomes[i], _, _ = ls_sample_candidates(
            states, weights, trace_bra, rng,
        )
    end

    return rho, outcomes
end

# ---------------------------------------------------------------------------
# Weak ZZ measurement
#
# K_m^ZZ = [I + (-1)^m lambda_zz Z_i Z_j]
#          / sqrt(2(1 + lambda_zz^2))
# ---------------------------------------------------------------------------

function ls_zz_measurement_candidates(
    rho::MPS,
    sites,
    i::Int,
    j::Int,
    lambda_zz::Float64,
    trace_bra::MPS;
    maxdim::Int,
    cutoff::Float64,
)
    @assert j == i + 1
    @assert 0.0 <= lambda_zz <= 1.0

    bra_i, ket_i = 2i - 1, 2i
    bra_j, ket_j = 2j - 1, 2j
    denom = sqrt(2.0 * (1.0 + lambda_zz^2))

    I_bra =
        op(LS_I2, sites[bra_i]) *
        op(LS_I2, sites[bra_j])

    ZZ_bra =
        op(LS_Z, sites[bra_i]) *
        op(LS_Z, sites[bra_j])

    I_ket =
        op(conj(LS_I2), sites[ket_i]) *
        op(conj(LS_I2), sites[ket_j])

    ZZ_ket =
        op(conj(LS_Z), sites[ket_i]) *
        op(conj(LS_Z), sites[ket_j])

    states = Vector{MPS}(undef, 2)
    weights = zeros(Float64, 2)

    for outcome in 0:1
        s = (-1.0)^outcome

        K_bra = (I_bra + s * lambda_zz * ZZ_bra) / denom
        K_ket = (I_ket + s * lambda_zz * ZZ_ket) / denom

        # Four doubled sites are contiguous for a nearest-neighbor
        # physical bond: (i_bra, i_ket, j_bra, j_ket).
        supergate = K_bra * K_ket
        rho_m = apply([supergate], rho; cutoff=cutoff, maxdim=maxdim)

        states[outcome + 1] = rho_m
        weights[outcome + 1] = ls_trace(rho_m, trace_bra)
    end

    return states, weights
end

function ls_apply_weak_zz_layer(
    rho::MPS,
    sites,
    L::Int,
    lambda_zz::Float64,
    trace_bra::MPS,
    rng::AbstractRNG;
    maxdim::Int,
    cutoff::Float64,
)
    bonds = ls_physical_bonds(L)
    outcomes = Vector{Int}(undef, length(bonds))

    for (n, (i, j)) in enumerate(bonds)
        states, weights = ls_zz_measurement_candidates(
            rho, sites, i, j, lambda_zz, trace_bra;
            maxdim=maxdim, cutoff=cutoff,
        )

        rho, outcomes[n], _, _ = ls_sample_candidates(
            states, weights, trace_bra, rng,
        )
    end

    return rho, outcomes
end

# ---------------------------------------------------------------------------
# Dephasing channels
#
# X channel:  rho -> (1-q_x) rho + q_x X rho X
# ZZ channel: rho -> (1-q_zz) rho + q_zz ZZ rho ZZ
# ---------------------------------------------------------------------------

function ls_build_x_dephasing_gates(sites, L::Int, q_x::Float64)
    @assert 0.0 <= q_x <= 0.5

    gates = Vector{ITensor}(undef, L)
    for i in 1:L
        bra_i, ket_i = 2i - 1, 2i

        identity_gate =
            op(LS_I2, sites[bra_i]) *
            op(conj(LS_I2), sites[ket_i])

        flip_gate =
            op(LS_X, sites[bra_i]) *
            op(conj(LS_X), sites[ket_i])

        gates[i] = (1.0 - q_x) * identity_gate + q_x * flip_gate
    end
    return gates
end

function ls_build_zz_dephasing_gates(sites, L::Int, q_zz::Float64)
    @assert 0.0 <= q_zz <= 0.5

    bonds = ls_physical_bonds(L)
    gates = Vector{ITensor}(undef, length(bonds))

    for (n, (i, j)) in enumerate(bonds)
        bra_i, ket_i = 2i - 1, 2i
        bra_j, ket_j = 2j - 1, 2j

        identity_gate =
            op(LS_I2, sites[bra_i]) *
            op(LS_I2, sites[bra_j]) *
            op(conj(LS_I2), sites[ket_i]) *
            op(conj(LS_I2), sites[ket_j])

        zz_gate =
            op(LS_Z, sites[bra_i]) *
            op(LS_Z, sites[bra_j]) *
            op(conj(LS_Z), sites[ket_i]) *
            op(conj(LS_Z), sites[ket_j])

        gates[n] = (1.0 - q_zz) * identity_gate + q_zz * zz_gate
    end
    return gates
end

function ls_apply_channel_layer(
    rho::MPS,
    gates::Vector{ITensor},
    trace_bra::MPS;
    maxdim::Int,
    cutoff::Float64,
)
    for gate in gates
        rho = apply([gate], rho; cutoff=cutoff, maxdim=maxdim)
    end
    return ls_trace_normalize(rho, trace_bra)
end

# ---------------------------------------------------------------------------
# One trajectory at a general point on the lambda path
# ---------------------------------------------------------------------------

function ls_evolve_one_trial(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    q_x::Float64,
    q_zz::Float64,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    seed::Int,
)
    @assert L >= 2
    @assert 0.0 <= lambda_x <= 1.0
    @assert 0.0 <= lambda_zz <= 1.0
    @assert 0.0 <= q_x <= 0.5
    @assert 0.0 <= q_zz <= 0.5
    @assert T_max >= 0

    rng = MersenneTwister(seed)

    sites = siteinds("Qubit", 2L)
    trace_bra = ls_trace_state(sites)
    rho = ls_trace_normalize(ls_initial_state(sites), trace_bra)

    x_dephasing_gates = ls_build_x_dephasing_gates(sites, L, q_x)
    zz_dephasing_gates = ls_build_zz_dephasing_gates(sites, L, q_zz)

    max_trace_error = abs(ls_trace(rho, trace_bra) - 1.0)
    max_linkdim = ls_max_interphysical_linkdim(rho, L)

    for _ in 1:T_max
        # 1. Weak X measurements.
        rho, _ = ls_apply_weak_x_layer(
            rho, sites, L, lambda_x, trace_bra, rng;
            maxdim=maxdim, cutoff=cutoff,
        )

        # 2. X dephasing.
        rho = ls_apply_channel_layer(
            rho, x_dephasing_gates, trace_bra;
            maxdim=maxdim, cutoff=cutoff,
        )

        # 3. Weak ZZ measurements.
        rho, _ = ls_apply_weak_zz_layer(
            rho, sites, L, lambda_zz, trace_bra, rng;
            maxdim=maxdim, cutoff=cutoff,
        )

        # 4. ZZ dephasing.
        rho = ls_apply_channel_layer(
            rho, zz_dephasing_gates, trace_bra;
            maxdim=maxdim, cutoff=cutoff,
        )

        max_trace_error = max(
            max_trace_error,
            abs(ls_trace(rho, trace_bra) - 1.0),
        )
        max_linkdim = max(
            max_linkdim,
            ls_max_interphysical_linkdim(rho, L),
        )
    end

    return (
        rho=rho,
        sites=sites,
        trace_bra=trace_bra,
        max_trace_error=max_trace_error,
        max_interphysical_linkdim=max_linkdim,
    )
end

# ---------------------------------------------------------------------------
# ITensorCorrelators helper
# ---------------------------------------------------------------------------

function ls_sum_correlators_batched(
    psi::MPS,
    operators::Tuple,
    site_tuples::AbstractVector;
    batch_size::Int=4096,
)
    isempty(site_tuples) && return 0.0
    @assert batch_size >= 1

    total = 0.0
    for first_index in 1:batch_size:length(site_tuples)
        last_index = min(
            first_index + batch_size - 1,
            length(site_tuples),
        )
        batch = site_tuples[first_index:last_index]
        values_dict = correlator(psi, operators, batch)

        for value in values(values_dict)
            total += real(value)
        end
    end
    return total
end

# ---------------------------------------------------------------------------
# Panel (b): Rényi-2 susceptibility
#
# r_i = Z_(2i-1) Z_(2i)
# kappa_2 = [L + 2 sum_{i<j}<r_i r_j>] / L
# ---------------------------------------------------------------------------

function ls_kappa2_one_trajectory(
    rho::MPS,
    sites,
    L::Int;
    pair_batch_size::Int=4096,
    norm_tol::Float64=1e-14,
)
    hs_norm = ls_hs_norm(rho)
    if !isfinite(hs_norm) || hs_norm <= norm_tol
        error("Invalid Hilbert-Schmidt norm: <rho|rho> = $hs_norm")
    end

    psi = deepcopy(rho / sqrt(hs_norm))

    pair_sites = NTuple{4,Int}[]
    sizehint!(pair_sites, L * (L - 1) ÷ 2)

    for i in 1:(L - 1)
        for j in (i + 1):L
            push!(pair_sites, (2i - 1, 2i, 2j - 1, 2j))
        end
    end

    pair_sum = ls_sum_correlators_batched(
        psi,
        ("Z", "Z", "Z", "Z"),
        pair_sites;
        batch_size=pair_batch_size,
    )

    Q2 = Float64(L) + 2.0 * pair_sum
    kappa2 = Q2 / L
    M2 = Q2 / L^2

    return (
        kappa2=kappa2,
        M2=M2,
        Q2_expectation=Q2,
        pair_sum=pair_sum,
        purity=hs_norm,
    )
end

# ---------------------------------------------------------------------------
# Panel (a): Edwards-Anderson susceptibility
#
# kappa_EA = [L + 2 sum_{i<j}<Z_i Z_j>^2] / L
#
# Here <...> is the ordinary trace expectation within one trajectory.
# ---------------------------------------------------------------------------

function ls_physical_zz_expectation(
    rho::MPS,
    trace_bra::MPS,
    sites,
    i::Int,
    j::Int,
)
    bra_i, bra_j = 2i - 1, 2j - 1

    # One-site Z gates do not increase MPS bond dimension.
    rho_zz = apply(
        [
            op(LS_Z, sites[bra_i]),
            op(LS_Z, sites[bra_j]),
        ],
        rho;
        cutoff=0.0,
    )

    return real(inner(trace_bra, rho_zz))
end

function ls_kappaEA_one_trajectory(
    rho::MPS,
    trace_bra::MPS,
    sites,
    L::Int,
)
    pair_square_sum = 0.0

    for i in 1:(L - 1)
        for j in (i + 1):L
            cij = ls_physical_zz_expectation(
                rho, trace_bra, sites, i, j,
            )
            pair_square_sum += cij^2
        end
    end

    kappaEA = (Float64(L) + 2.0 * pair_square_sum) / L

    return (
        kappaEA=kappaEA,
        pair_square_sum=pair_square_sum,
    )
end

# ---------------------------------------------------------------------------
# Average over Born-sampled trajectories at one (L, lambda) point
# ---------------------------------------------------------------------------

function ls_run_lambda_point(
    L::Int,
    lambda::Real;
    delta::Float64=0.7,
    q::Float64=0.1,
    ntrials::Int=20,
    T_max::Int=2L,
    maxdim::Int=256,
    cutoff::Float64=1e-10,
    pair_batch_size::Int=4096,
    seed::Int=1234,
    compute_EA::Bool=true,
)
    lambda = Float64(lambda)

    @assert 0.0 <= lambda <= 1.0
    @assert 0.0 <= delta <= 1.0
    @assert 0.0 <= q <= 0.5
    @assert ntrials >= 1

    lambda_x = delta * lambda
    lambda_zz = delta * (1.0 - lambda)

    master_rng = MersenneTwister(seed)

    kappa2_trials = Vector{Float64}(undef, ntrials)
    kappaEA_trials = fill(NaN, ntrials)
    purity_trials = Vector{Float64}(undef, ntrials)
    linkdims = Vector{Int}(undef, ntrials)
    trace_errors = Vector{Float64}(undef, ntrials)

    for trial in 1:ntrials
        trial_seed = Int(rand(master_rng, UInt32))

        evolved = ls_evolve_one_trial(
            L;
            lambda_x=lambda_x,
            lambda_zz=lambda_zz,
            q_x=q,
            q_zz=q,
            T_max=T_max,
            maxdim=maxdim,
            cutoff=cutoff,
            seed=trial_seed,
        )

        obs2 = ls_kappa2_one_trajectory(
            evolved.rho,
            evolved.sites,
            L;
            pair_batch_size=pair_batch_size,
        )

        kappa2_trials[trial] = obs2.kappa2
        purity_trials[trial] = obs2.purity
        linkdims[trial] = evolved.max_interphysical_linkdim
        trace_errors[trial] = evolved.max_trace_error

        if compute_EA
            obsEA = ls_kappaEA_one_trajectory(
                evolved.rho,
                evolved.trace_bra,
                evolved.sites,
                L,
            )
            kappaEA_trials[trial] = obsEA.kappaEA
        end
    end

    kappa2_mean = mean(kappa2_trials)
    kappa2_se = ntrials > 1 ?
        std(kappa2_trials; corrected=true) / sqrt(ntrials) : NaN

    kappaEA_mean = compute_EA ? mean(kappaEA_trials) : NaN
    kappaEA_se = compute_EA && ntrials > 1 ?
        std(kappaEA_trials; corrected=true) / sqrt(ntrials) : NaN

    return (
        L=L,
        lambda=lambda,
        lambda_x=lambda_x,
        lambda_zz=lambda_zz,
        delta=delta,
        q=q,
        kappa2=kappa2_mean,
        kappa2_se=kappa2_se,
        kappaEA=kappaEA_mean,
        kappaEA_se=kappaEA_se,
        purity=mean(purity_trials),
        ntrials=ntrials,
        T_max=T_max,
        max_interphysical_linkdim=maximum(linkdims),
        max_trace_error=maximum(trace_errors),
        kappa2_trials=kappa2_trials,
        kappaEA_trials=kappaEA_trials,
    )
end

# ---------------------------------------------------------------------------
# Full scan
# ---------------------------------------------------------------------------

function ls_run_scan(
    L_values,
    lambda_values;
    delta::Float64=0.7,
    q::Float64=0.1,
    ntrials::Int=20,
    T_factor::Int=2,
    maxdim::Int=256,
    cutoff::Float64=1e-10,
    pair_batch_size::Int=4096,
    seed::Int=1234,
    compute_EA::Bool=true,
    verbose::Bool=true,
)
    results = NamedTuple[]

    for L in L_values
        for lambda in lambda_values
            point_seed =
                seed +
                100000 * Int(L) +
                round(Int, 10000 * Float64(lambda))

            verbose && println(
                "L=$L, lambda=$(round(Float64(lambda), digits=4)), " *
                "lambda_x=$(round(delta * Float64(lambda), digits=4)), " *
                "lambda_zz=$(round(delta * (1-Float64(lambda)), digits=4))",
            )

            result = ls_run_lambda_point(
                Int(L),
                Float64(lambda);
                delta=delta,
                q=q,
                ntrials=ntrials,
                T_max=T_factor * Int(L),
                maxdim=maxdim,
                cutoff=cutoff,
                pair_batch_size=pair_batch_size,
                seed=point_seed,
                compute_EA=compute_EA,
            )
            push!(results, result)
        end
    end

    return results
end
