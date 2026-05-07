# ============================================================
# Rényi-2 Binder using MPO on doubled density-matrix MPS
#
# Assumes doubled-site ordering:
#   (1_bra,1_ket,2_bra,2_ket,...,L_bra,L_ket)
#
# Define (expectation-value-like observable):
#   q_i = Z_(2i-1)   (bra layer only)
#   Q   = sum_i q_i
#
# Then
#   M2 = <ρ|Q^2|ρ> / (L^2 <ρ|ρ>)
#   M4 = <ρ|Q^4|ρ> / (L^4 <ρ|ρ>)
#   B  = 1 - M4 / (3 M2^2)
#
# Sanity checks:
#   - All-up pure state: B_R = 2/3
#   - Maximally mixed: B_R ≈ 2/(3L) → 0
#
# This avoids:
#   - mixed_to_vectorized(...)
#   - explicit 8-point correlators
#   - repeated-index combinatorics
# ============================================================

using Random
using Statistics
using ITensors, ITensorMPS

export renyi2_binder_density_matrix
export renyi2_binder_correlator_doubled
export renyi2_order_mpo_one_layer_doubled

"""
    renyi2_binder_correlator_doubled(ρ::MPS, L::Int; layer::Symbol=:bra)

Compute Rényi-2 Binder using MPO observable on doubled MPS.
Measurement on bra layer only: Q = Σ_i Z_(2i-1).

ρ: Normalized doubled MPS with 2L sites (bra_i, ket_i pairs)
L: Number of physical sites
layer: :bra or :ket for which layer to measure on

Returns: (B, M2, M4, purity)
"""
function renyi2_binder_correlator_doubled(
    ρ::MPS,
    L::Int;
    layer::Symbol=:bra,
)
    sites = siteinds(ρ)
    @assert length(sites) == 2L "Expected doubled MPS with 2L sites"

    # Normalize trace
    tr_val = doubledtrace(ρ)
    ρ_norm = ρ / tr_val

    # Hilbert-Schmidt norm = Tr(ρ²)
    purity = real(inner(ρ_norm, ρ_norm))

    # Build MPO for Q = Σ_i Z_i on one layer
    Q = renyi2_order_mpo_one_layer_doubled(sites, L; opname="Z", layer=layer)

    # ψ1 = Q|ρ>
    ψ1 = apply(Q, ρ_norm; cutoff=1e-12, maxdim=256)
    num2 = real(inner(ψ1, ψ1))   # = <ρ|Q²|ρ>

    # ψ2 = Q²|ρ>
    ψ2 = apply(Q, ψ1; cutoff=1e-12, maxdim=256)
    num4 = real(inner(ψ2, ψ2))   # = <ρ|Q⁴|ρ>

    M2 = num2 / (L^2 * purity)
    M4 = num4 / (L^4 * purity)

    B = if abs(M2) > sqrt(eps(Float64))
        1.0 - M4 / (3.0 * M2^2)
    else
        NaN
    end

    return (
        B = B,
        M2 = M2,
        M4 = M4,
        purity = purity,
    )
end

"""
    max_linkdim_mps(ψ::MPS)

Safe helper for bond-dimension diagnostics.
"""
function max_linkdim_mps(ψ::MPS)
    N = length(ψ)
    N <= 1 && return 1

    dims = Int[]
    for j in 1:(N - 1)
        l = commonind(ψ[j], ψ[j + 1])
        if l !== nothing
            push!(dims, dim(l))
        end
    end
    return isempty(dims) ? 1 : maximum(dims)
end

"""
    renyi2_order_mpo_doubled(sites, L; opname="Z")

Build the MPO for

    Q = Σ_i Z_(2i-1) Z_(2i)

on the doubled MPS.
"""
function renyi2_order_mpo_doubled(sites, L::Int; opname::String="Z")
    @assert length(sites) == 2L "Expected doubled MPS with 2L sites"

    os = OpSum()
    for i in 1:L
        bra_idx = 2i - 1
        ket_idx = 2i
        os += 1.0, opname, bra_idx, opname, ket_idx
    end
    return MPO(os, sites)
end

"""
    renyi2_order_mpo_one_layer_doubled(sites, L; opname="Z", layer=:bra)

Build the MPO for

    Q = Σ_i Z_i

acting on only one layer of the doubled density-matrix MPS.

For ordering:
    (1_bra, 1_ket, 2_bra, 2_ket, ..., L_bra, L_ket)

use:
    layer = :bra  -> sites 1,3,5,...
    layer = :ket  -> sites 2,4,6,...

This is the "expectation-value-like" Rényi-2 observable.
Maximally mixed gives B_R ≈ 2/(3L) → 0.
"""
function renyi2_order_mpo_one_layer_doubled(
    sites,
    L::Int;
    opname::String="Z",
    layer::Symbol=:bra,
)
    @assert length(sites) == 2L "Expected doubled MPS with 2L sites"
    @assert layer in (:bra, :ket) "layer must be :bra or :ket"

    os = OpSum()
    for i in 1:L
        site_idx = layer === :bra ? 2i - 1 : 2i
        os += 1.0, opname, site_idx
    end
    return MPO(os, sites)
end

"""
    renyi2_moments_doubled(
        ρ::MPS,
        L::Int;
        Q=nothing,
        cutoff=1e-12,
        maxdim=256,
        warn_threshold=1e-10,
        invalid_threshold=1e-14,
        verbose=false,
    )

Compute Rényi-2 moments directly from doubled density-matrix MPS.

Important:
- `ρ` should already be normalized by trace:
      ρ = ρ / doubledtrace(ρ)

Definitions:
    purity = <ρ|ρ> = Tr(ρ†ρ)
    M2 = <ρ|Q^2|ρ> / (L^2 * purity)
    M4 = <ρ|Q^4|ρ> / (L^4 * purity)

Returns:
    (M2, M4, purity, diagnostics)
"""
function renyi2_moments_doubled(
    ρ::MPS,
    L::Int;
    Q::Union{Nothing,MPO}=nothing,
    cutoff::Float64=1e-12,
    maxdim::Int=256,
    warn_threshold::Float64=1e-10,
    invalid_threshold::Float64=1e-14,
    verbose::Bool=false,
)
    sites = siteinds(ρ)
    @assert length(sites) == 2L "Expected doubled MPS with 2L sites"

    if Q === nothing
        # Use one-layer observable by default
        Q = renyi2_order_mpo_one_layer_doubled(sites, L; opname="Z", layer=:bra)
    end

    # Hilbert-Schmidt norm of the density matrix.
    # For Hermitian ρ, this equals Tr(ρ²).
    purity = real(inner(ρ, ρ))

    reliability = if purity <= invalid_threshold
        :invalid
    elseif purity <= warn_threshold
        :warning
    else
        :ok
    end

    if reliability === :invalid
        diagnostics = (
            purity_raw = purity,
            reliability = reliability,
            num2 = NaN,
            num4 = NaN,
            max_link_dim_psi1 = missing,
            max_link_dim_psi2 = missing,
            cutoff_used = cutoff,
            maxdim_used = maxdim,
        )
        return NaN, NaN, purity, diagnostics
    elseif reliability === :warning
        @warn "Very small purity detected: $purity. Rényi-2 moments may be unreliable."
    end

    # ψ1 = Q|ρ>
    ψ1 = apply(Q, ρ; cutoff=cutoff, maxdim=maxdim)
    num2 = real(inner(ψ1, ψ1))   # = <ρ|Q²|ρ>

    # ψ2 = Q²|ρ>
    ψ2 = apply(Q, ψ1; cutoff=cutoff, maxdim=maxdim)
    num4 = real(inner(ψ2, ψ2))   # = <ρ|Q⁴|ρ>

    M2 = num2 / (L^2 * purity)
    M4 = num4 / (L^4 * purity)

    diagnostics = (
        purity_raw = purity,
        reliability = reliability,
        num2 = num2,
        num4 = num4,
        max_link_dim_psi1 = max_linkdim_mps(ψ1),
        max_link_dim_psi2 = max_linkdim_mps(ψ2),
        cutoff_used = cutoff,
        maxdim_used = maxdim,
    )

    if verbose
        println("Purity = $purity")
        println("<ρ|Q²|ρ> = $num2")
        println("<ρ|Q⁴|ρ> = $num4")
        println("M2 = $M2")
        println("M4 = $M4")
        println("max link dim ψ1 = $(diagnostics.max_link_dim_psi1)")
        println("max link dim ψ2 = $(diagnostics.max_link_dim_psi2)")
    end

    return M2, M4, purity, diagnostics
end

"""
    renyi2_binder_doubled(ρ::MPS, L; kwargs...)

Single-trajectory Rényi-2 Binder from doubled MPS.

Returns:
    (B, M2, M4, purity, diagnostics)
"""
function renyi2_binder_doubled(
    ρ::MPS,
    L::Int;
    Q::Union{Nothing,MPO}=nothing,
    cutoff::Float64=1e-12,
    maxdim::Int=256,
    warn_threshold::Float64=1e-10,
    invalid_threshold::Float64=1e-14,
    verbose::Bool=false,
)
    M2, M4, purity, diagnostics = renyi2_moments_doubled(
        ρ, L;
        Q=Q,
        cutoff=cutoff,
        maxdim=maxdim,
        warn_threshold=warn_threshold,
        invalid_threshold=invalid_threshold,
        verbose=verbose,
    )

    B = if isfinite(M2) && isfinite(M4) && abs(M2) > sqrt(eps(Float64))
        1.0 - M4 / (3.0 * M2^2)
    else
        NaN
    end

    return (
        B = B,
        M2 = M2,
        M4 = M4,
        purity = purity,
        diagnostics = diagnostics,
    )
end

"""
    renyi2_binder_density_matrix(L::Int; kwargs...)

ORIGINAL Rényi-2 Binder function (UNCHANGED from prior implementation).

Uses MPO-based moment calculation on doubled density-matrix MPS.
Evolved density matrix ρ is in doubled form with 2L sites (bra, ket pairs).
We measure Z on the bra layer only (odd indices: 1, 3, 5, ..., 2L-1).
"""
function renyi2_binder_density_matrix(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    P_x::Float64=0.0,
    P_zz::Float64=0.0,
    ntrials::Int=200,
    maxdim::Int=256,
    cutoff::Float64=1e-12,
    seed::Union{Nothing,Int}=nothing,
    use_optimized::Bool=true,
    verbose::Bool=false,
)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)

    M2s = fill(NaN, ntrials)
    M4s = fill(NaN, ntrials)
    Bs = fill(NaN, ntrials)
    purities = fill(NaN, ntrials)

    for t in 1:ntrials
        # Evolve one trial using original dynamics (use_optimized parameter supported)
        state, _ = if use_optimized
            evolve_density_matrix_one_trial_new(
                L;
                lambda_x=lambda_x,
                lambda_zz=lambda_zz,
                P_x=P_x,
                P_zz=P_zz,
                maxdim=maxdim,
                cutoff=cutoff,
                rng=rng,
            )
        else
            evolve_density_matrix_one_trial(
                L;
                lambda_x=lambda_x,
                lambda_zz=lambda_zz,
                P_x=P_x,
                P_zz=P_zz,
                maxdim=maxdim,
                cutoff=cutoff,
                rng=rng,
            )
        end

        # Extract doubled MPS and compute moments using original MPO method
        ρ = get_mps(state)
        res = renyi2_binder_correlator_doubled(ρ, L; layer=:bra)
        
        M2s[t] = res.M2
        M4s[t] = res.M4
        Bs[t] = res.B
        purities[t] = res.purity
        
        if verbose && (t == 1 || t == ntrials || t % max(1, ntrials ÷ 10) == 0)
            println("trial $t/$ntrials: B=$(res.B), M2=$(res.M2), M4=$(res.M4), purity=$(res.purity)")
        end
    end

    valid = isfinite.(M2s) .& isfinite.(M4s) .& isfinite.(Bs)
    n_valid = count(valid)
    n_invalid = ntrials - n_valid

    if n_valid == 0
        return (
            B = NaN,
            B_mean_of_trials = NaN,
            B_std_of_trials = NaN,
            M2_bar = NaN,
            M4_bar = NaN,
            purity_bar = NaN,
            M2s = M2s,
            M4s = M4s,
            Bs = Bs,
            purities = purities,
            ntrials = ntrials,
            n_valid = 0,
            n_invalid = n_invalid,
        )
    end

    M2_bar = mean(M2s[valid])
    M4_bar = mean(M4s[valid])
    purity_bar = mean(purities[valid])

    # Ensemble-averaged Binder: B = 1 - <M4> / (3 <M2>²)
    B_ens = 1.0 - M4_bar / (3.0 * M2_bar^2)

    return (
        B = B_ens,
        B_mean_of_trials = mean(Bs[valid]),
        B_std_of_trials = std(Bs[valid]),
        M2_bar = M2_bar,
        M4_bar = M4_bar,
        purity_bar = purity_bar,
        M2s = M2s,
        M4s = M4s,
        Bs = Bs,
        purities = purities,
        ntrials = ntrials,
        n_valid = n_valid,
        n_invalid = n_invalid,
    )
end