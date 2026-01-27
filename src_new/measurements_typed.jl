"""
Type-safe measurement operations with multiple dispatch.

This module implements measurement protocols that dispatch based on state type
(PureStateMPS, DiagonalStateMPS, MixedStateMPS) following the architecture from src_1.
"""

using Random
using ITensors, ITensorMPS

include("types.jl")
using .Main: StateMPS, PureStateMPS, DiagonalStateMPS, MixedStateMPS
using .Main: get_mps, normalize_state, M_bra, doubledtrace

export expval, measure, measure_with_outcome

# ========================================
# Expectation Values
# ========================================

"""
    expval(state::PureStateMPS, M::AbstractMatrix, pos::Int; kwargs...)

Compute expectation value ⟨ψ|M|ψ⟩ for a pure state.
"""
function expval(state::PureStateMPS, M::AbstractMatrix, pos::Int; 
                cutoff=1E-8, maxdim=200, refs=0)
    ψ = get_mps(state)
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Mψ = apply(op(M, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...), ψ; 
               cutoff=cutoff, maxdim=maxdim)
    val = inner(ψ, Mψ) / norm(ψ)^2
    return val
end

"""
    expval(state::DiagonalStateMPS, M::AbstractMatrix, pos::Int; kwargs...)

Compute expectation value Tr(ρM) for a diagonal density matrix.
"""
function expval(state::DiagonalStateMPS, M::AbstractMatrix, pos::Int; 
                cutoff=1E-8, maxdim=200, refs=0)
    ψ = get_mps(state)
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Mψ = apply(op(M, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...), ψ; 
               cutoff=cutoff, maxdim=maxdim)
    val = inner(MPS(sites, i -> "+"), Mψ) / inner(MPS(sites, i -> "+"), ψ)
    return val
end

"""
    expval(state::MixedStateMPS, M::AbstractMatrix, pos::Int; kwargs...)

Compute expectation value Tr(ρM) for a mixed state density matrix.
"""
function expval(state::MixedStateMPS, M::AbstractMatrix, pos::Int; 
                cutoff=1E-8, maxdim=200, refs=0)
    ρ = get_mps(state)
    sites = siteinds(ρ)
    
    bra = M_bra(sites, M, pos; refs=refs)
    val = inner(bra, ρ) / doubledtrace(ρ)
    return val
end

# ========================================
# Measurements with Fixed Outcome
# ========================================

"""
    measure(state::PureStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; kwargs...)

Apply weak measurement with fixed outcome m to a pure state.

# Arguments
- `state`: Pure quantum state
- `M`: Measurement operator
- `λ`: Measurement strength (0 ≤ λ ≤ 1)
- `pos`: Position to measure
- `m`: Measurement outcome (false=0, true=1)

# Returns
- Normalized state after measurement
"""
function measure(state::PureStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; 
                 cutoff=1E-8, maxdim=200, refs=0)
    ψ = get_mps(state)
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    # Kraus operator: Π = (I + (-1)^m * λ*M) / √(2(1+λ²))
    Π = (I + (-1)^m * λ*M) / sqrt(2*(1+λ^2))
    g = op(Π, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...)

    ψ_new = apply(g, ψ; cutoff=cutoff, maxdim=maxdim)
    return normalize_state(PureStateMPS(ψ_new))
end

"""
    measure(state::DiagonalStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; kwargs...)

Apply weak measurement with fixed outcome m to a diagonal state.
For diagonal states, apply Π†Π (squared Kraus operator).
"""
function measure(state::DiagonalStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; 
                 cutoff=1E-8, maxdim=200, refs=0)
    ψ = get_mps(state)
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M) / sqrt(2*(1+λ^2))
    g = op(Π*Π, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...)

    ψ_new = apply(g, ψ; cutoff=cutoff, maxdim=maxdim)
    return normalize_state(DiagonalStateMPS(ψ_new))
end

"""
    measure(state::MixedStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; kwargs...)

Apply weak measurement with fixed outcome m to a mixed state.
Apply Kraus operator to both bra and ket parts of the doubled MPS.
"""
function measure(state::MixedStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; 
                 cutoff=1E-8, maxdim=200, refs=0)
    ρ = get_mps(state)
    sites = siteinds(ρ)
    M_width = Int(log2(size(M)[1]))
    L = length(sites)÷2 - refs

    Π = (I + (-1)^m * λ*M) / sqrt(2*(1+λ^2))
    g1 = op(Π, [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]...)
    g2 = op(Π, [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]...)

    ρ_new = apply([g1, g2], ρ; cutoff=cutoff, maxdim=maxdim)
    return normalize_state(MixedStateMPS(ρ_new))
end

# ========================================
# Measurements with Sampled Outcome
# ========================================

"""
    measure_with_outcome(state::StateMPS, M::AbstractMatrix, λ::Float64, pos::Int; kwargs...)

Perform weak measurement and sample the outcome based on Born rule.

# Returns
- `new_state`: State after measurement
- `outcome::Bool`: Sampled measurement outcome (false=0, true=1)
- `expval`: Expectation value before measurement
"""
function measure_with_outcome(state::StateMPS, M::AbstractMatrix, λ::Float64, pos::Int; 
                              kwargs...)
    val = expval(state, M, pos; kwargs...)
    # Born rule probability: p(0) = (1 + 2λ/(1+λ²)⟨M⟩) / 2
    prob = (1 + 2λ/(1+λ^2)*real(val)) / 2
    
    if rand() < abs(prob)
        return measure(state, M, λ, pos, false; kwargs...), false, val
    else
        return measure(state, M, λ, pos, true; kwargs...), true, val
    end
end

# ========================================
# Batch Measurements
# ========================================

"""
    measure(state::StateMPS, M::AbstractMatrix, λ::Float64, positions::AbstractVector, 
            outcomes::Vector{Bool}; kwargs...)

Apply measurements with fixed outcomes at multiple positions.
"""
function measure(state::StateMPS, M::AbstractMatrix, λ::Float64, 
                 positions::AbstractVector, outcomes::Vector{Bool}; kwargs...)
    for (i, pos) in enumerate(positions)
        state = measure(state, M, λ, pos, outcomes[i]; kwargs...)
    end
    return state
end

"""
    measure_with_outcome(state::StateMPS, M::AbstractMatrix, λ::Float64, 
                        positions::AbstractVector; kwargs...)

Apply measurements and sample outcomes at multiple positions.

# Returns
- `final_state`: State after all measurements
- `outcomes::Vector{Bool}`: Sampled measurement outcomes
- `expvals::Vector{ComplexF64}`: Expectation values before each measurement
"""
function measure_with_outcome(state::StateMPS, M::AbstractMatrix, λ::Float64,
                              positions::AbstractVector; kwargs...)
    outcomes = zeros(Bool, length(positions))
    expvals = zeros(ComplexF64, length(positions))

    for (i, pos) in enumerate(positions)
        state, outcome, val = measure_with_outcome(state, M, λ, pos; kwargs...)
        outcomes[i] = outcome
        expvals[i] = val
    end
    
    return state, outcomes, expvals
end

# ========================================
# Operator Construction Helper
# ========================================

"""
    op(M::AbstractMatrix, sites::Index...)

Create an ITensor operator from a matrix and site indices.
"""
function op(M::AbstractMatrix, sites::Index...)
    # Use ITensors built-in operator construction
    n = length(sites)
    d = size(M, 1)
    @assert d == 2^n "Matrix dimension must match number of sites"
    
    # Create operator tensor
    T = ITensor(M, [s' for s in sites]..., sites...)
    return T
end
