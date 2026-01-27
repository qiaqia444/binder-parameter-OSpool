"""
State construction utilities for quantum many-body simulations.

This module provides functions for creating initial MPS states in various configurations.
Supports multiple state representations: PureStateMPS, DiagonalStateMPS, MixedStateMPS.
"""

using ITensors, ITensorMPS

include("types.jl")
using .Main: PureStateMPS, DiagonalStateMPS, MixedStateMPS, StateMPS

export create_up_state_mps, create_plus_state_mps, create_product_mps
export zero_state, ghz_state, neel_state

"""
    create_up_state_mps(L::Int) -> (MPS, Vector{Index})

Create a product state MPS with all spins pointing up (|↑↑...↑⟩).

# Arguments
- `L::Int`: System size (number of sites)

# Returns
- `ψ`: MPS in the all-up state
- `sites`: Site indices for the MPS
"""
function create_up_state_mps(L::Int)
    sites = siteinds("S=1/2", L)
    ψ = productMPS(sites, fill("Up", L))
    return ψ, sites
end

"""
    create_plus_state_mps(L::Int) -> (MPS, Vector{Index})

Create a product state MPS with all spins in |+⟩ state (equal superposition).

# Arguments
- `L::Int`: System size (number of sites)

# Returns
- `ψ`: MPS in the all-plus state
- `sites`: Site indices for the MPS
"""
function create_plus_state_mps(L::Int)
    sites = siteinds("S=1/2", L)
    ψ = productMPS(sites, fill("+", L))
    return ψ, sites
end

"""
    create_product_mps(L::Int, state::String) -> (MPS, Vector{Index})

Create a product state MPS with all spins in the specified state.

# Arguments
- `L::Int`: System size (number of sites)
- `state::String`: Initial state for each site (e.g., "Up", "Dn", "+", "-")

# Returns
- `ψ`: MPS in the specified product state
- `sites`: Site indices for the MPS
"""

# ========================================
# Type-Safe State Constructors
# ========================================

"""
    zero_state(::Type{T}, L::Int) where {T<:StateMPS}

Create a state with all spins in |0⟩ (down) state.
"""
function zero_state(::Type{PureStateMPS}, L::Int)
    sites = siteinds("S=1/2", L)
    ψ = productMPS(sites, fill("Dn", L))
    return PureStateMPS(ψ)
end

function zero_state(::Type{DiagonalStateMPS}, L::Int)
    sites = siteinds("S=1/2", L)
    ψ = productMPS(sites, fill("Dn", L))
    return DiagonalStateMPS(ψ)
end

function zero_state(::Type{MixedStateMPS}, L::Int)
    sites = siteinds("S=1/2", 2L)
    ψ = productMPS(sites, fill("Dn", 2L))
    return MixedStateMPS(ψ)
end

"""
    ghz_state(::Type{T}, L::Int; ref=false) where {T<:StateMPS}

Create a GHZ state: (|00...0⟩ + |11...1⟩)/√2
"""
function ghz_state(::Type{PureStateMPS}, L::Int; ref=false)
    N = L + ref
    sites = siteinds("S=1/2", N)
    ψ0 = productMPS(sites, fill("Dn", N))
    ψ1 = productMPS(sites, fill("Up", N))
    return PureStateMPS((ψ0 + ψ1) / sqrt(2))
end

function ghz_state(::Type{MixedStateMPS}, L::Int; ref=false)
    N = 2L + 2ref
    sites = siteinds("S=1/2", N)
    ρ00 = productMPS(sites, fill("Dn", N))
    ρ11 = productMPS(sites, fill("Up", N))
    ρ01 = productMPS(sites, [iseven(i) ? "Dn" : "Up" for i in 1:N])
    ρ10 = productMPS(sites, [iseven(i) ? "Up" : "Dn" for i in 1:N])
    return MixedStateMPS((ρ00 + ρ01 + ρ10 + ρ11) / 2)
end

"""
    neel_state(::Type{T}, L::Int) where {T<:StateMPS}

Create a Neel state with alternating up/down spins: |↑↓↑↓...⟩
"""
function neel_state(::Type{PureStateMPS}, L::Int)
    sites = siteinds("S=1/2", L)
    state_pattern = [iseven(i) ? "Dn" : "Up" for i in 1:L]
    ψ = productMPS(sites, state_pattern)
    return PureStateMPS(ψ)
end

function neel_state(::Type{DiagonalStateMPS}, L::Int)
    sites = siteinds("S=1/2", L)
    state_pattern = [iseven(i) ? "Dn" : "Up" for i in 1:L]
    ψ = productMPS(sites, state_pattern)
    return DiagonalStateMPS(ψ)
end

function neel_state(::Type{MixedStateMPS}, L::Int)
    sites = siteinds("S=1/2", 2L)
    state_pattern = [iseven(i) ? "Dn" : "Up" for i in 1:2L]
    ψ = productMPS(sites, state_pattern)
    return MixedStateMPS(ψ)
end
function create_product_mps(L::Int, state::String)
    sites = siteinds("S=1/2", L)
    ψ = productMPS(sites, fill(state, L))
    return ψ, sites
end
