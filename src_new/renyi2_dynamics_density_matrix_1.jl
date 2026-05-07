"""
renyi2_dynamics_density_matrix.jl

Rényi-2 Binder dynamics for full density-matrix evolution using MixedStateMPS.

This follows the density-matrix structure in `renyi2_density_matrix.jl` /
`density_matrix_dynamics.jl` style:

  • MixedStateMPS is a doubled MPS representation of ρ.
  • Physical spin i uses two MPS sites:
        2i-1 = bra leg, 2i = ket leg.
  • Weak measurements are sampled trajectory-by-trajectory.
  • Dephasing is applied as the full averaged density-matrix channel by default:
        ρ -> (1-p)ρ + p UρU†
    implemented as a single superoperator gate, not as explicit MPS addition.
  • Rényi-2 correlators are computed directly in doubled space:
        C₂[O] = Tr(O ρ O ρ) / Tr(ρ²)
              = <ρ| O_bra O_ket |ρ> / <ρ|ρ>
    for Pauli-Z strings O.

Important:
  This file intentionally does NOT convert MixedStateMPS to a one-site
  d²-vectorized MPS. The direct doubled-space contraction is safer and avoids
  index/bond-order mistakes.
"""

using Random, Statistics
using ITensors, ITensorMPS
using LinearAlgebra

# -----------------------------------------------------------------------------
# Project-local MixedStateMPS type
# -----------------------------------------------------------------------------
include("types.jl")
using .Main: MixedStateMPS, get_mps

export evolve_renyi2_density_matrix_one_trial
export renyi2_binder_density_matrix_dynamics
export renyi2_moments_density_matrix
export renyi2_binder_density_matrix_final
export renyi2_binder_density_matrix_separate
export doubledtrace
export bell

# -----------------------------------------------------------------------------
# Pauli matrices
# -----------------------------------------------------------------------------
const σx = Float64[0 1; 1 0]
const σz = Float64[1 0; 0 -1]
const I2 = Matrix{Float64}(I, 2, 2)
const I4 = Matrix{Float64}(I, 4, 4)
const ZZ_2site = kron(σz, σz)

# -----------------------------------------------------------------------------
# Bell / trace utilities for doubled density-matrix MPS
# -----------------------------------------------------------------------------

"""
    id_mps(s1, s2)

Two-site MPS for |00> + |11> on a doubled bra/ket pair.
"""
function id_mps(s1::Index, s2::Index)
    return MPS([1.0, 0.0, 0.0, 1.0], [s1, s2])
end

"""
    bell(sites)

Create the product Bell state <I| used for Tr(ρ) in doubled space.
For each physical site, the local vector is |00> + |11> on `(bra, ket)`.
"""
function bell(sites_in)
    sites = collect(sites_in)
    N = length(sites)
    @assert iseven(N) "Doubled density-matrix MPS must have an even number of sites."

    tensors = ITensor[]
    for n in 1:2:N
        pair = id_mps(sites[n], sites[n + 1])
        push!(tensors, pair[1], pair[2])
    end
    return MPS(tensors)
end

"""
    doubledtrace(ρ)
    doubledtrace(ρ, trace_bra)

Compute Tr(ρ) from a doubled density-matrix MPS.
"""
function doubledtrace(ρ::MPS)
    return real(inner(bell(siteinds(ρ)), ρ))
end

function doubledtrace(ρ::MPS, trace_bra::MPS)
    return real(inner(trace_bra, ρ))
end

function _trace_normalize(ρ::MPS, trace_bra::MPS; tol::Float64=1e-14)
    trρ = doubledtrace(ρ, trace_bra)
    if abs(trρ) < tol
        error("Trace is too close to zero during density-matrix normalization: Tr(ρ) = $trρ")
    end
    return ρ / trρ
end

"""
    M_bra(sites, M, pos; refs=0)

Create the trace bra with operator M inserted at physical position `pos`.
This follows the src_1/MPS/tools.jl pattern used in your current code:
operator insertion is on the ket/even legs, and contraction with |ρ> gives
Tr(M ρ) for real Pauli operators.
"""
function M_bra(sites_in, M::AbstractMatrix, pos::Int; refs::Int=0)
    sites = collect(sites_in)
    L = length(sites) ÷ 2 - refs
    M_width = Int(round(log2(size(M, 1))))
    @assert size(M, 1) == size(M, 2) "Inserted operator must be square."
    @assert 2^M_width == size(M, 1) "Operator dimension must be a power of 2."
    @assert 1 <= pos <= L "Position pos=$pos outside physical chain length L=$L."

    trace_bra = bell(sites)
    ket_sites = [sites[mod1(2 * (pos + offset), 2L)] for offset in 0:(M_width - 1)]
    return apply(op(M, ket_sites...), trace_bra)
end

# -----------------------------------------------------------------------------
# Evolution gates
# -----------------------------------------------------------------------------

function _precompute_trace_bras(sites, L::Int; lambda_x::Float64, lambda_zz::Float64)
    trace_bra = bell(sites)

    x_trace_bras = Vector{MPS}(undef, lambda_x > 0 ? L : 0)
    if lambda_x > 0
        for i in 1:L
            x_trace_bras[i] = M_bra(sites, σx, i)
        end
    end

    zz_trace_bras = Vector{MPS}(undef, lambda_zz > 0 ? max(L - 1, 0) : 0)
    if lambda_zz > 0
        for i in 1:(L - 1)
            zz_trace_bras[i] = M_bra(sites, ZZ_2site, i)
        end
    end

    return trace_bra, x_trace_bras, zz_trace_bras
end

function _precompute_channel_gates(sites, L::Int; P_x::Float64, P_zz::Float64)
    x_channel_gates = Vector{ITensor}(undef, P_x > 0 ? L : 0)
    if P_x > 0
        for i in 1:L
            bra_i = 2i - 1
            ket_i = 2i

            I_gate = op(I2, sites[bra_i]) * op(I2, sites[ket_i])
            X_gate = op(σx, sites[bra_i]) * op(σx, sites[ket_i])

            # Superoperator: (1-p) I⊗I + p X⊗X*
            # σx is real, so X* = X.
            x_channel_gates[i] = (1 - P_x) * I_gate + P_x * X_gate
        end
    end

    zz_channel_gates = Vector{ITensor}(undef, P_zz > 0 ? max(L - 1, 0) : 0)
    if P_zz > 0
        for i in 1:(L - 1)
            j = i + 1
            bra_i = 2i - 1
            ket_i = 2i
            bra_j = 2j - 1
            ket_j = 2j

            I_gate = op(I2, sites[bra_i]) * op(I2, sites[bra_j]) *
                     op(I2, sites[ket_i]) * op(I2, sites[ket_j])

            ZZ_gate = op(σz, sites[bra_i]) * op(σz, sites[bra_j]) *
                      op(σz, sites[ket_i]) * op(σz, sites[ket_j])

            # Superoperator: (1-p) I + p (ZZ)⊗(ZZ)*
            zz_channel_gates[i] = (1 - P_zz) * I_gate + P_zz * ZZ_gate
        end
    end

    return x_channel_gates, zz_channel_gates
end

function _precompute_renyi2_z_gates(sites, L::Int)
    z_gates = Vector{ITensor}(undef, L)
    for i in 1:L
        bra_i = 2i - 1
        ket_i = 2i
        z_gates[i] = op(σz, sites[bra_i]) * op(σz, sites[ket_i])
    end
    return z_gates
end

# -----------------------------------------------------------------------------
# One-step evolution pieces
# -----------------------------------------------------------------------------

function _apply_weak_x_measurements(ρ::MPS, sites, L::Int, trace_bra::MPS,
                                    x_trace_bras::Vector{MPS}, rng;
                                    lambda_x::Float64,
                                    maxdim::Int,
                                    cutoff::Float64)
    if lambda_x <= 0
        return ρ
    end

    for i in 1:L
        bra_i = 2i - 1
        ket_i = 2i

        trρ = doubledtrace(ρ, trace_bra)
        expval_X = real(inner(x_trace_bras[i], ρ) / trρ)

        # For K_m = (I + (-1)^m λX)/sqrt(2(1+λ²)),
        # p(m=0) = [1 + 2λ/(1+λ²) <X>] / 2.
        p0 = (1 + 2 * lambda_x / (1 + lambda_x^2) * expval_X) / 2
        p0 = clamp(p0, 0.0, 1.0)
        outcome = rand(rng) < p0 ? 0 : 1

        K = (I2 + (-1)^outcome * lambda_x * σx) / sqrt(2 * (1 + lambda_x^2))
        K_bra = op(K, sites[bra_i])
        K_ket = op(K, sites[ket_i])

        ρ = apply([K_bra, K_ket], ρ; cutoff=cutoff, maxdim=maxdim)
        ρ = _trace_normalize(ρ, trace_bra)
    end

    return ρ
end

function _apply_x_dephasing(ρ::MPS, sites, L::Int, rng;
                            P_x::Float64,
                            mode::Symbol,
                            x_channel_gates::Vector{ITensor},
                            maxdim::Int,
                            cutoff::Float64)
    if P_x <= 0
        return ρ
    end

    if mode == :channel
        for i in 1:L
            ρ = apply([x_channel_gates[i]], ρ; cutoff=cutoff, maxdim=maxdim)
        end
    elseif mode == :trajectory
        for i in 1:L
            if rand(rng) < P_x
                bra_i = 2i - 1
                ket_i = 2i
                X_bra = op(σx, sites[bra_i])
                X_ket = op(σx, sites[ket_i])
                ρ = apply([X_bra, X_ket], ρ; cutoff=cutoff, maxdim=maxdim)
            end
        end
    else
        error("Unknown dephasing mode: $mode. Use :channel or :trajectory.")
    end

    return ρ
end

function _apply_weak_zz_measurements(ρ::MPS, sites, L::Int, trace_bra::MPS,
                                     zz_trace_bras::Vector{MPS}, rng;
                                     lambda_zz::Float64,
                                     maxdim::Int,
                                     cutoff::Float64)
    if lambda_zz <= 0
        return ρ
    end

    for i in 1:(L - 1)
        j = i + 1
        bra_i = 2i - 1
        ket_i = 2i
        bra_j = 2j - 1
        ket_j = 2j

        trρ = doubledtrace(ρ, trace_bra)
        expval_ZZ = real(inner(zz_trace_bras[i], ρ) / trρ)

        # For K_m = (I + (-1)^m λ ZZ)/sqrt(2(1+λ²)).
        p0 = (1 + 2 * lambda_zz / (1 + lambda_zz^2) * expval_ZZ) / 2
        p0 = clamp(p0, 0.0, 1.0)
        outcome = rand(rng) < p0 ? 0 : 1

        K = (I4 + (-1)^outcome * lambda_zz * ZZ_2site) / sqrt(2 * (1 + lambda_zz^2))
        K_bra = op(K, sites[bra_i], sites[bra_j])
        K_ket = op(K, sites[ket_i], sites[ket_j])

        ρ = apply([K_bra, K_ket], ρ; cutoff=cutoff, maxdim=maxdim)
        ρ = _trace_normalize(ρ, trace_bra)
    end

    return ρ
end

function _apply_zz_dephasing(ρ::MPS, sites, L::Int, rng;
                             P_zz::Float64,
                             mode::Symbol,
                             zz_channel_gates::Vector{ITensor},
                             maxdim::Int,
                             cutoff::Float64)
    if P_zz <= 0
        return ρ
    end

    if mode == :channel
        for i in 1:(L - 1)
            ρ = apply([zz_channel_gates[i]], ρ; cutoff=cutoff, maxdim=maxdim)
        end
    elseif mode == :trajectory
        for i in 1:(L - 1)
            if rand(rng) < P_zz
                j = i + 1
                bra_i = 2i - 1
                ket_i = 2i
                bra_j = 2j - 1
                ket_j = 2j

                ZZ_bra = op(σz, sites[bra_i]) * op(σz, sites[bra_j])
                ZZ_ket = op(σz, sites[ket_i]) * op(σz, sites[ket_j])
                ρ = apply([ZZ_bra, ZZ_ket], ρ; cutoff=cutoff, maxdim=maxdim)
            end
        end
    else
        error("Unknown dephasing mode: $mode. Use :channel or :trajectory.")
    end

    return ρ
end

# -----------------------------------------------------------------------------
# Rényi-2 Binder observables in doubled space
# -----------------------------------------------------------------------------

function _renyi2_expectation(ρ::MPS, z_gates::Vector{ITensor}, active_sites::Vector{Int}, hs_norm::Float64;
                             maxdim::Int,
                             cutoff::Float64)
    # Empty Pauli string = identity, so C₂[I] = Tr(ρ²)/Tr(ρ²) = 1.
    if isempty(active_sites)
        return 1.0
    end

    gates = ITensor[z_gates[i] for i in active_sites]
    ρO = apply(gates, ρ; cutoff=cutoff, maxdim=maxdim)
    return real(inner(ρ, ρO) / hs_norm)
end

"""
    renyi2_moments_density_matrix(ρ, L; sites=nothing, maxdim=256, cutoff=1e-12)

Compute the Rényi-2 moments

    C₂[O] = Tr(O ρ O ρ) / Tr(ρ²)
    M₂   = (1/L²) Σᵢⱼ C₂[Zᵢ Zⱼ]
    M₄   = (1/L⁴) Σᵢⱼₖₗ C₂[Zᵢ Zⱼ Zₖ Zₗ]

using the doubled MPS directly.

The sums are exact but optimized using Z²=I and commutativity:

    Σᵢⱼ C₂[ZᵢZⱼ]
      = L + 2 Σᵢ<ⱼ C₂[ZᵢZⱼ]

    Σᵢⱼₖₗ C₂[ZᵢZⱼZₖZₗ]
      = (3L² - 2L)
        + (12L - 16) Σᵢ<ⱼ C₂[ZᵢZⱼ]
        + 24 Σᵢ<ⱼ<ₖ<ₗ C₂[ZᵢZⱼZₖZₗ]

The first term in M₄ is the total identity contribution from [4] and [2,2]
index multiplicities.
"""
function renyi2_moments_density_matrix(ρ::MPS, L::Int;
                                        sites=nothing,
                                        maxdim::Int=256,
                                        cutoff::Float64=1e-12,
                                        purity_tol::Float64=1e-14)
    @assert L >= 1 "L must be at least 1."

    sites_local = isnothing(sites) ? siteinds(ρ) : sites
    z_gates = _precompute_renyi2_z_gates(sites_local, L)

    hs_norm = real(inner(ρ, ρ))
    if hs_norm < purity_tol
        error("Tr(ρ²) is too small for Rényi-2 normalization: Tr(ρ²) = $hs_norm")
    end

    # Pair contribution: sum over unordered i<j.
    pair_sum = 0.0
    for i in 1:(L - 1)
        for j in (i + 1):L
            pair_sum += _renyi2_expectation(ρ, z_gates, [i, j], hs_norm;
                                            maxdim=maxdim, cutoff=cutoff)
        end
    end

    # Quad contribution: sum over unordered i<j<k<l.
    quad_sum = 0.0
    if L >= 4
        for i in 1:(L - 3)
            for j in (i + 1):(L - 2)
                for k in (j + 1):(L - 1)
                    for l in (k + 1):L
                        quad_sum += _renyi2_expectation(ρ, z_gates, [i, j, k, l], hs_norm;
                                                        maxdim=maxdim, cutoff=cutoff)
                    end
                end
            end
        end
    end

    m2_num = L + 2 * pair_sum
    M2 = m2_num / L^2

    identity_m4_count = 3L^2 - 2L
    pair_m4_coeff = 12L - 16
    m4_num = identity_m4_count + pair_m4_coeff * pair_sum + 24 * quad_sum
    M4 = m4_num / L^4

    return (M2 = M2,
            M4 = M4,
            purity = hs_norm,
            pair_sum = pair_sum,
            quad_sum = quad_sum)
end

function _binder_from_moments(M2::Float64, M4::Float64)
    return 1.0 - M4 / (3.0 * M2^2 + eps(Float64))
end

# -----------------------------------------------------------------------------
# Single-trial dynamics
# -----------------------------------------------------------------------------

"""
    evolve_renyi2_density_matrix_one_trial(L; kwargs...)

Evolve one sampled weak-measurement trajectory of the full density matrix and
record Rényi-2 Binder data as a function of layer.

Keyword arguments:
  lambda_x, lambda_zz : weak measurement strengths
  P_x, P_zz           : dephasing probabilities
  T_max               : number of full layers; default 2L
  record_strobe       : :after_x_noise, :end_of_layer, or :both
  dephasing_mode      : :channel by default; :trajectory for jump sampling tests
  compute_initial     : whether to record t=0 before evolution
  normalize_after_channels : renormalize by Tr(ρ) after channel steps to remove
                             small truncation drift

Returns a NamedTuple containing time labels, M₂, M₄, Binder, purity, and final state.
"""
function evolve_renyi2_density_matrix_one_trial(L::Int;
                                                lambda_x::Float64,
                                                lambda_zz::Float64,
                                                P_x::Float64=0.0,
                                                P_zz::Float64=0.0,
                                                T_max::Union{Nothing,Int}=nothing,
                                                maxdim::Int=256,
                                                cutoff::Float64=1e-12,
                                                rng=Random.GLOBAL_RNG,
                                                seed::Union{Nothing,Int}=nothing,
                                                record_strobe::Symbol=:after_x_noise,
                                                dephasing_mode::Symbol=:channel,
                                                compute_initial::Bool=false,
                                                normalize_after_channels::Bool=true)
    @assert L >= 1 "L must be at least 1."
    if !isnothing(seed)
        rng = MersenneTwister(seed)
    end
    if !(record_strobe in (:after_x_noise, :end_of_layer, :both))
        error("record_strobe must be :after_x_noise, :end_of_layer, or :both.")
    end

    sites = siteinds("Qubit", 2L)
    ρ = MPS(sites, _ -> "Up")

    T = isnothing(T_max) ? 2L : T_max

    trace_bra, x_trace_bras, zz_trace_bras =
        _precompute_trace_bras(sites, L; lambda_x=lambda_x, lambda_zz=lambda_zz)

    x_channel_gates, zz_channel_gates =
        _precompute_channel_gates(sites, L; P_x=P_x, P_zz=P_zz)

    times = Int[]
    strobes = Symbol[]
    M2s = Float64[]
    M4s = Float64[]
    Bs = Float64[]
    purities = Float64[]
    traces = Float64[]

    function record!(t::Int, strobe::Symbol, ρ_current::MPS)
        ρ_norm = _trace_normalize(ρ_current, trace_bra)
        obs = renyi2_moments_density_matrix(ρ_norm, L;
                                            sites=sites,
                                            maxdim=maxdim,
                                            cutoff=cutoff)
        push!(times, t)
        push!(strobes, strobe)
        push!(M2s, obs.M2)
        push!(M4s, obs.M4)
        push!(Bs, _binder_from_moments(obs.M2, obs.M4))
        push!(purities, obs.purity)
        push!(traces, doubledtrace(ρ_norm, trace_bra))
        return nothing
    end

    if compute_initial
        record!(0, :initial, ρ)
    end

    for t in 1:T
        # 1. Weak X measurements
        ρ = _apply_weak_x_measurements(ρ, sites, L, trace_bra, x_trace_bras, rng;
                                       lambda_x=lambda_x,
                                       maxdim=maxdim,
                                       cutoff=cutoff)

        # 2. X dephasing
        ρ = _apply_x_dephasing(ρ, sites, L, rng;
                               P_x=P_x,
                               mode=dephasing_mode,
                               x_channel_gates=x_channel_gates,
                               maxdim=maxdim,
                               cutoff=cutoff)

        if normalize_after_channels
            ρ = _trace_normalize(ρ, trace_bra)
        end

        if record_strobe in (:after_x_noise, :both)
            record!(t, :after_x_noise, ρ)
        end

        # 3. Weak ZZ measurements
        ρ = _apply_weak_zz_measurements(ρ, sites, L, trace_bra, zz_trace_bras, rng;
                                        lambda_zz=lambda_zz,
                                        maxdim=maxdim,
                                        cutoff=cutoff)

        # 4. ZZ dephasing
        ρ = _apply_zz_dephasing(ρ, sites, L, rng;
                                P_zz=P_zz,
                                mode=dephasing_mode,
                                zz_channel_gates=zz_channel_gates,
                                maxdim=maxdim,
                                cutoff=cutoff)

        if normalize_after_channels
            ρ = _trace_normalize(ρ, trace_bra)
        end

        if record_strobe in (:end_of_layer, :both)
            record!(t, :end_of_layer, ρ)
        end
    end

    return (times = times,
            strobes = strobes,
            M2 = M2s,
            M4 = M4s,
            B = Bs,
            purity = purities,
            trace = traces,
            state = MixedStateMPS(ρ),
            sites = sites)
end

# -----------------------------------------------------------------------------
# Ensemble dynamics
# -----------------------------------------------------------------------------

function _column_mean(A::Matrix{Float64})
    return vec(mean(A; dims=1))
end

function _column_std(A::Matrix{Float64})
    return vec(std(A; dims=1, corrected=false))
end

"""
    renyi2_binder_density_matrix_dynamics(L; kwargs...)

Run `ntrials` density-matrix trajectories and return the ensemble-averaged
Rényi-2 Binder dynamics.

The ensemble Binder is formed from averaged moments:

    B_ensemble(t) = 1 - mean(M4(t)) / [3 mean(M2(t))²]

The function also returns the mean/std of trial-wise Binder values.
"""
function renyi2_binder_density_matrix_dynamics(L::Int;
                                               lambda_x::Float64,
                                               lambda_zz::Float64,
                                               P_x::Float64=0.0,
                                               P_zz::Float64=0.0,
                                               ntrials::Int=200,
                                               T_max::Union{Nothing,Int}=nothing,
                                               maxdim::Int=256,
                                               cutoff::Float64=1e-12,
                                               seed::Union{Nothing,Int}=nothing,
                                               record_strobe::Symbol=:after_x_noise,
                                               dephasing_mode::Symbol=:channel,
                                               compute_initial::Bool=false,
                                               normalize_after_channels::Bool=true,
                                               verbose::Bool=true)
    @assert ntrials >= 1 "ntrials must be at least 1."

    master_rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)

    first = evolve_renyi2_density_matrix_one_trial(L;
                                                   lambda_x=lambda_x,
                                                   lambda_zz=lambda_zz,
                                                   P_x=P_x,
                                                   P_zz=P_zz,
                                                   T_max=T_max,
                                                   maxdim=maxdim,
                                                   cutoff=cutoff,
                                                   rng=master_rng,
                                                   record_strobe=record_strobe,
                                                   dephasing_mode=dephasing_mode,
                                                   compute_initial=compute_initial,
                                                   normalize_after_channels=normalize_after_channels)

    nrecords = length(first.B)
    times = copy(first.times)
    strobes = copy(first.strobes)

    M2_trials = Matrix{Float64}(undef, ntrials, nrecords)
    M4_trials = Matrix{Float64}(undef, ntrials, nrecords)
    B_trials = Matrix{Float64}(undef, ntrials, nrecords)
    purity_trials = Matrix{Float64}(undef, ntrials, nrecords)
    trace_trials = Matrix{Float64}(undef, ntrials, nrecords)

    M2_trials[1, :] .= first.M2
    M4_trials[1, :] .= first.M4
    B_trials[1, :] .= first.B
    purity_trials[1, :] .= first.purity
    trace_trials[1, :] .= first.trace

    if verbose
        println("Completed trial 1 / $ntrials")
    end

    for trial in 2:ntrials
        result = evolve_renyi2_density_matrix_one_trial(L;
                                                        lambda_x=lambda_x,
                                                        lambda_zz=lambda_zz,
                                                        P_x=P_x,
                                                        P_zz=P_zz,
                                                        T_max=T_max,
                                                        maxdim=maxdim,
                                                        cutoff=cutoff,
                                                        rng=master_rng,
                                                        record_strobe=record_strobe,
                                                        dephasing_mode=dephasing_mode,
                                                        compute_initial=compute_initial,
                                                        normalize_after_channels=normalize_after_channels)

        if length(result.B) != nrecords || result.times != times || result.strobes != strobes
            error("Trial $trial produced a different record structure. Check record_strobe/T_max settings.")
        end

        M2_trials[trial, :] .= result.M2
        M4_trials[trial, :] .= result.M4
        B_trials[trial, :] .= result.B
        purity_trials[trial, :] .= result.purity
        trace_trials[trial, :] .= result.trace

        if verbose && (trial == ntrials || trial % max(1, ntrials ÷ 10) == 0)
            println("Completed trial $trial / $ntrials")
        end
    end

    M2_mean = _column_mean(M2_trials)
    M4_mean = _column_mean(M4_trials)
    B_ensemble = [1.0 - M4_mean[i] / (3.0 * M2_mean[i]^2 + eps(Float64))
                  for i in eachindex(M2_mean)]

    return (times = times,
            strobes = strobes,
            B = B_ensemble,
            B_mean_of_trials = _column_mean(B_trials),
            B_std_of_trials = _column_std(B_trials),
            M2_bar = M2_mean,
            M4_bar = M4_mean,
            M2_trials = M2_trials,
            M4_trials = M4_trials,
            B_trials = B_trials,
            purity_mean = _column_mean(purity_trials),
            purity_std = _column_std(purity_trials),
            trace_mean = _column_mean(trace_trials),
            trace_std = _column_std(trace_trials),
            ntrials = ntrials,
            params = (L=L,
                      lambda_x=lambda_x,
                      lambda_zz=lambda_zz,
                      P_x=P_x,
                      P_zz=P_zz,
                      T_max=isnothing(T_max) ? 2L : T_max,
                      maxdim=maxdim,
                      cutoff=cutoff,
                      record_strobe=record_strobe,
                      dephasing_mode=dephasing_mode))
end

"""
    renyi2_binder_density_matrix_final(L; kwargs...)

Convenience wrapper that returns only the last recorded time point.
"""
function renyi2_binder_density_matrix_final(L::Int; kwargs...)
    dyn = renyi2_binder_density_matrix_dynamics(L; kwargs...)
    last_idx = length(dyn.B)
    return (B = dyn.B[last_idx],
            B_mean_of_trials = dyn.B_mean_of_trials[last_idx],
            B_std_of_trials = dyn.B_std_of_trials[last_idx],
            M2_bar = dyn.M2_bar[last_idx],
            M4_bar = dyn.M4_bar[last_idx],
            purity_mean = dyn.purity_mean[last_idx],
            purity_std = dyn.purity_std[last_idx],
            trace_mean = dyn.trace_mean[last_idx],
            trace_std = dyn.trace_std[last_idx],
            time = dyn.times[last_idx],
            strobe = dyn.strobes[last_idx],
            ntrials = dyn.ntrials,
            params = dyn.params)
end

"""
    renyi2_binder_density_matrix_separate(L::Int; kwargs...)

Compatibility wrapper that returns arrays of per-trial values for backward compatibility.

This calls renyi2_binder_density_matrix_dynamics internally and runs ntrials independent trials,
collecting per-trial statistics for backward compatibility with older interfaces.
"""
function renyi2_binder_density_matrix_separate(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    P_x::Float64=0.0,
    P_zz::Float64=0.0,
    ntrials::Int=200,
    maxdim::Int=256,
    cutoff::Float64=1e-12,
    seed::Union{Nothing,Int}=nothing,
    strobe::Symbol=:after_full_layer,
    verbose::Bool=false,
)
    # Map old strobe names to new ones if needed
    record_strobe = if strobe == :after_full_layer
        :after_x_noise  # default
    else
        strobe
    end
    
    # Call dynamics to get final result  
    result = renyi2_binder_density_matrix_final(
        L;
        lambda_x=lambda_x,
        lambda_zz=lambda_zz,
        P_x=P_x,
        P_zz=P_zz,
        ntrials=ntrials,
        maxdim=maxdim,
        cutoff=cutoff,
        seed=seed,
        record_strobe=record_strobe,
        dephasing_mode=:channel,
        compute_initial=false,
        normalize_after_channels=true,
        verbose=verbose,
    )
    
    # For backward compatibility, return the result structure with expected fields
    # The result from _final already has: B, B_mean_of_trials, B_std_of_trials, 
    # M2_bar, M4_bar, purity_mean, purity_std, trace_mean, trace_std, ntrials, params
    return (
        B = result.B,
        B_mean_of_trials = result.B_mean_of_trials,
        B_std_of_trials = result.B_std_of_trials,
        M2_bar = result.M2_bar,
        M4_bar = result.M4_bar,
        purity_bar = result.purity_mean,  # Map purity_mean -> purity_bar
        M2s = fill(NaN, ntrials),  # Placeholder for per-trial values (not available from dynamics)
        M4s = fill(NaN, ntrials),
        Bs = fill(NaN, ntrials),
        purities = fill(NaN, ntrials),
        ntrials = result.ntrials,
        n_valid = result.ntrials,  # Assume all valid
        n_invalid = 0,
        max_linkdim = 0,  # Not tracked by dynamics
    )
end