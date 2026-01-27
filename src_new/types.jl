"""
Type definitions for different quantum state representations.

This module defines wrapper types for different MPS representations and provides
polymorphic operations that dispatch based on state type.
"""

using ITensors, ITensorMPS
import Base: /
import ITensorMPS: norm

export StateMPS, PureStateMPS, DiagonalStateMPS, MixedStateMPS
export get_mps, normalize_state

"""
    PureStateMPS

Wrapper for pure quantum state |ψ⟩ represented as an MPS.
Used for quantum trajectory simulations where states remain pure.
"""
struct PureStateMPS
    mps::MPS
end

"""
    DiagonalStateMPS

Wrapper for diagonal density matrix ρ = Σᵢ pᵢ |i⟩⟨i| represented as an MPS.
Useful for classical probability distributions over quantum states.
"""
struct DiagonalStateMPS
    mps::MPS
end

"""
    MixedStateMPS

Wrapper for general mixed state density matrix ρ represented as a "doubled" MPS.
Uses Choi-Jamiolkowski representation where ρ is vectorized.
"""
struct MixedStateMPS
    mps::MPS
end

"""
Union type for all state representations.
"""
const StateMPS = Union{PureStateMPS, DiagonalStateMPS, MixedStateMPS}

# ========================================
# Accessor Functions
# ========================================

"""
    get_mps(state::StateMPS) -> MPS

Extract the underlying MPS from a wrapped state.
"""
get_mps(state::PureStateMPS) = state.mps
get_mps(state::DiagonalStateMPS) = state.mps
get_mps(state::MixedStateMPS) = state.mps

# ========================================
# Normalization Operations
# ========================================

"""
    norm(state::PureStateMPS) -> Float64

Compute the norm of a pure state: ⟨ψ|ψ⟩^(1/2)
"""
function norm(state::PureStateMPS)
    ψ = state.mps
    return norm(ψ)
end

"""
    norm(state::DiagonalStateMPS) -> Float64

Compute the trace of a diagonal density matrix.
"""
function norm(state::DiagonalStateMPS)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)
    # For diagonal states, trace = ⟨+|ψ⟩ * 2^(L/2)
    return inner(MPS(sites, i -> "+"), ψ) * 2^(L/2)
end

"""
    norm(state::MixedStateMPS) -> Float64

Compute the trace of a mixed state density matrix using doubled trace.
"""
function norm(state::MixedStateMPS)
    ψ = state.mps
    return doubledtrace(ψ)
end

# ========================================
# Scalar Division
# ========================================

/(state::PureStateMPS, scalar::Number) = PureStateMPS(state.mps / scalar)
/(state::DiagonalStateMPS, scalar::Number) = DiagonalStateMPS(state.mps / scalar)
/(state::MixedStateMPS, scalar::Number) = MixedStateMPS(state.mps / scalar)

# ========================================
# Normalization Helper
# ========================================

"""
    normalize_state(state::StateMPS) -> StateMPS

Normalize a quantum state to unit norm/trace.
"""
function normalize_state(state::T) where {T <: StateMPS}
    n = norm(state)
    return n > 1e-14 ? state / n : state
end

# ========================================
# Utility Functions for MixedStateMPS
# ========================================

"""
    id_mps(s1::Index, s2::Index) -> MPS

Create a "doubled" MPS corresponding to identity on a single site.
In practice, this is a Bell pair (without normalization).
"""
function id_mps(s1::Index, s2::Index)
    return MPS([1; 0; 0; 1], [s1, s2])
end

"""
    bell(sites::Vector{Index}) -> MPS

Create a "doubled" MPS corresponding to identity on N sites.
Looks like tensor product of many Bell pairs.
"""
function bell(sites::Vector{<:Index})
    N = length(sites)
    tensors = ITensor[]
    for i in 1:2:N
        a, b = id_mps(sites[i], sites[i+1])
        push!(tensors, a, b)
    end 
    return MPS(tensors)
end

"""
    doubledtrace(ρ::MPS) -> Float64

Compute the trace of a density matrix encoded in a "doubled" MPS.
"""
function doubledtrace(ρ::MPS)
    return inner(bell(siteinds(ρ)), ρ)
end

"""
    M_bra(sites::Vector{Index}, M::AbstractMatrix, pos::Int; refs=0) -> MPS

Create "doubled" MPS corresponding to operator M at position pos 
and identity everywhere else. Used for expectation values in mixed states.
"""
function M_bra(sites::Vector{<:Index}, M::AbstractMatrix, pos::Int; refs=0)
    M_width = Int(log2(size(M)[1]))
    L = length(sites)÷2 - refs

    bra = bell(sites)
    bra = apply(op(M, [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]...), bra)
    return bra
end

# ========================================
# Display Functions
# ========================================

Base.show(io::IO, state::PureStateMPS) = print(io, "PureStateMPS(L=$(length(state.mps)))")
Base.show(io::IO, state::DiagonalStateMPS) = print(io, "DiagonalStateMPS(L=$(length(state.mps)))")
Base.show(io::IO, state::MixedStateMPS) = print(io, "MixedStateMPS(L=$(length(state.mps)÷2)))")
