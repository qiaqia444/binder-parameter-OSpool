"""
Quantum channels for density matrix evolution.

This module implements dephasing and decoherence channels that operate on
density matrices (DiagonalStateMPS and MixedStateMPS).

For pure state quantum trajectories (PureStateMPS), dephasing is implemented
stochastically. For density matrices, ALL Kraus operators are applied with
proper weights.
"""

using ITensors, ITensorMPS
using Random
using LinearAlgebra

include("types.jl")
using .Main: PureStateMPS, DiagonalStateMPS, MixedStateMPS, StateMPS
using .Main: get_mps, normalize_state

export apply_x_dephasing_channel, apply_zz_dephasing_channel
export apply_dephasing_layer

# ========================================
# Pure State Trajectory (Stochastic)
# ========================================

"""
    apply_x_dephasing_channel(state::PureStateMPS, site::Int, P_x::Float64; kwargs...)

Apply X dephasing channel to a pure state using quantum trajectory method.
Stochastically samples which Kraus operator to apply.

For pure states: ℰ(|ψ⟩) = (1-P_x)|ψ⟩ or X|ψ⟩ (sampled)
"""
function apply_x_dephasing_channel(state::PureStateMPS, site::Int, P_x::Float64;
                                  maxdim::Int=256, cutoff::Float64=1e-12, 
                                  rng=Random.GLOBAL_RNG)
    P_x <= 0 && return state
    
    ψ = get_mps(state)
    sites = siteinds(ψ)
    
    if rand(rng) < P_x
        # Apply Kraus operator K₁ = √P_x · X
        X_gate = 2 * op("Sx", sites[site])  # Pauli X = 2*Sx
        ψ_new = apply(X_gate, ψ; maxdim=maxdim, cutoff=cutoff)
        return normalize_state(PureStateMPS(ψ_new))
    else
        # Apply Kraus operator K₀ = √(1-P_x) · I (identity)
        return state
    end
end

"""
    apply_zz_dephasing_channel(state::PureStateMPS, i::Int, j::Int, P_zz::Float64; kwargs...)

Apply ZZ dephasing channel to a pure state using quantum trajectory method.
"""
function apply_zz_dephasing_channel(state::PureStateMPS, i::Int, j::Int, P_zz::Float64;
                                   maxdim::Int=256, cutoff::Float64=1e-12,
                                   rng=Random.GLOBAL_RNG)
    P_zz <= 0 && return state
    abs(i-j) != 1 && return state
    
    ψ = get_mps(state)
    sites = siteinds(ψ)
    
    if rand(rng) < P_zz
        # Apply Kraus operator K₁ = √P_zz · (Z⊗Z)
        Z_i = 2 * op("Sz", sites[i])
        Z_j = 2 * op("Sz", sites[j])
        ZZ = Z_i * Z_j
        ψ_new = apply(ZZ, ψ; maxdim=maxdim, cutoff=cutoff)
        return normalize_state(PureStateMPS(ψ_new))
    else
        return state
    end
end

# ========================================
# Density Matrix Evolution (Deterministic)
# ========================================

"""
    apply_x_dephasing_channel(state::DiagonalStateMPS, site::Int, P_x::Float64; kwargs...)

Apply X dephasing channel to a diagonal density matrix.

For diagonal states where X²=I:
    ρ → (1-P_x)ρ + P_x·X·ρ = [(1-P_x)I + P_x·X]·ρ

This applies BOTH Kraus operators with proper weights (not sampling).
"""
function apply_x_dephasing_channel(state::DiagonalStateMPS, site::Int, P_x::Float64;
                                  maxdim::Int=256, cutoff::Float64=1e-12,
                                  rng=nothing)
    P_x <= 0 && return state
    
    ψ = get_mps(state)
    sites = siteinds(ψ)
    
    # For diagonal states: ρ → (1-p)ρ + p·X·ρ
    # Since X² = I, this simplifies to: [(1-p)I + p·X]·ρ
    X = 2 * op("Sx", sites[site])
    Id = op("Id", sites[site])
    
    gate = (1-P_x)*Id + P_x*X
    ψ_new = apply(gate, ψ; maxdim=maxdim, cutoff=cutoff)
    
    return normalize_state(DiagonalStateMPS(ψ_new))
end

"""
    apply_zz_dephasing_channel(state::DiagonalStateMPS, i::Int, j::Int, P_zz::Float64; kwargs...)

Apply ZZ dephasing channel to a diagonal density matrix.

For diagonal states where (Z⊗Z)²=I:
    ρ → (1-P_zz)ρ + P_zz·(Z⊗Z)·ρ = [(1-P_zz)I + P_zz·(Z⊗Z)]·ρ
"""
function apply_zz_dephasing_channel(state::DiagonalStateMPS, i::Int, j::Int, P_zz::Float64;
                                   maxdim::Int=256, cutoff::Float64=1e-12,
                                   rng=nothing)
    P_zz <= 0 && return state
    abs(i-j) != 1 && return state
    
    ψ = get_mps(state)
    sites = siteinds(ψ)
    
    # For diagonal states: ρ → (1-p)ρ + p·(Z⊗Z)·ρ
    Z_i = 2 * op("Sz", sites[i])
    Z_j = 2 * op("Sz", sites[j])
    ZZ = Z_i * Z_j
    II = op("Id", sites[i]) * op("Id", sites[j])
    
    gate = (1-P_zz)*II + P_zz*ZZ
    ψ_new = apply(gate, ψ; maxdim=maxdim, cutoff=cutoff)
    
    return normalize_state(DiagonalStateMPS(ψ_new))
end

"""
    apply_x_dephasing_channel(state::MixedStateMPS, site::Int, P_x::Float64; kwargs...)

Apply X dephasing channel to a mixed state (doubled MPS representation).

For general mixed states:
    ρ → (1-P_x)ρ + P_x·X·ρ·X†
    
In doubled representation, apply to both bra and ket parts.
"""
function apply_x_dephasing_channel(state::MixedStateMPS, site::Int, P_x::Float64;
                                  maxdim::Int=256, cutoff::Float64=1e-12,
                                  rng=nothing)
    P_x <= 0 && return state
    
    ρ = get_mps(state)
    sites = siteinds(ρ)
    L = length(sites) ÷ 2
    
    # Apply X to both bra and ket indices
    X_bra = 2 * op("Sx", sites[2*site-1])
    X_ket = 2 * op("Sx", sites[2*site])
    I_bra = op("Id", sites[2*site-1])
    I_ket = op("Id", sites[2*site])
    
    gate = (1-P_x)*(I_bra*I_ket) + P_x*(X_bra*X_ket)
    ρ_new = apply(gate, ρ; maxdim=maxdim, cutoff=cutoff)
    
    return normalize_state(MixedStateMPS(ρ_new))
end

"""
    apply_zz_dephasing_channel(state::MixedStateMPS, i::Int, j::Int, P_zz::Float64; kwargs...)

Apply ZZ dephasing channel to a mixed state.
"""
function apply_zz_dephasing_channel(state::MixedStateMPS, i::Int, j::Int, P_zz::Float64;
                                   maxdim::Int=256, cutoff::Float64=1e-12,
                                   rng=nothing)
    P_zz <= 0 && return state
    abs(i-j) != 1 && return state
    
    ρ = get_mps(state)
    sites = siteinds(ρ)
    L = length(sites) ÷ 2
    
    # Apply ZZ to both bra and ket indices
    Z_i_bra = 2 * op("Sz", sites[2*i-1])
    Z_j_bra = 2 * op("Sz", sites[2*j-1])
    Z_i_ket = 2 * op("Sz", sites[2*i])
    Z_j_ket = 2 * op("Sz", sites[2*j])
    
    ZZ_bra = Z_i_bra * Z_j_bra
    ZZ_ket = Z_i_ket * Z_j_ket
    II_bra = op("Id", sites[2*i-1]) * op("Id", sites[2*j-1])
    II_ket = op("Id", sites[2*i]) * op("Id", sites[2*j])
    
    gate = (1-P_zz)*(II_bra*II_ket) + P_zz*(ZZ_bra*ZZ_ket)
    ρ_new = apply(gate, ρ; maxdim=maxdim, cutoff=cutoff)
    
    return normalize_state(MixedStateMPS(ρ_new))
end

# ========================================
# Layer-based Application (Like src_1)
# ========================================

"""
    apply_dephasing_layer(state::StateMPS, operator::String, p::Float64, 
                          positions::AbstractVector; kwargs...)

Apply dephasing channel to multiple positions.

# Arguments
- `state`: Quantum state (Pure, Diagonal, or Mixed)
- `operator`: "X" or "ZZ"
- `p`: Dephasing probability
- `positions`: Sites to apply dephasing (for ZZ, these are bond indices)

# Dispatch
- PureStateMPS: Stochastic (quantum trajectory)
- DiagonalStateMPS/MixedStateMPS: Deterministic (full density matrix)
"""
function apply_dephasing_layer(state::StateMPS, operator::String, p::Float64,
                               positions::AbstractVector; kwargs...)
    if operator == "X"
        for pos in positions
            state = apply_x_dephasing_channel(state, pos, p; kwargs...)
        end
    elseif operator == "ZZ"
        for i in positions
            state = apply_zz_dephasing_channel(state, i, i+1, p; kwargs...)
        end
    else
        error("Unknown operator: $operator. Use 'X' or 'ZZ'")
    end
    return state
end

