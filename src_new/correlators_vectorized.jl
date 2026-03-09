"""
ITensorCorrelators - Vectorized density matrix correlators for Z operators

Specialized correlator module for computing ⟨Z_i Z_j ... Z_k⟩ expectation values
from vectorized density matrix MPS representation.

Key features:
- Optimized for Z operators only (no fermion logic needed)
- Exploits Z^2 = I to reduce operator applications
- Canonicalizes tuples using commutativity of Z operators
- Computes Tr(Z_i Z_j ... ρ) via ⟨1|Z_i Z_j ... |ρ⟩

Vectorization convention:
- Each site has dimension d^2 (d=2 for qubits)
- Index mapping: α = i + (j-1)*d where i=row, j=column of ρ_ij
- Trace: ⟨1|ρ⟩ where ⟨1| = sum over diagonal entries α = i + (i-1)*d
"""
module ITensorCorrelators

using ITensors, ITensorMPS
using Statistics

export correlator, tr_rho, compute_M2_M4_vectorized, binder_from_trials

# ==============================================================================
# Vectorized Z operator for d^2-dimensional site indices
# ==============================================================================
"""
    make_vectorized_Z_op(s::Index)

Create Z operator for vectorized density matrix representation.

For vectorized |ρ⟩ with site index dimension d^2:
- Applies Z ⊗ I (left multiplication): |ρ⟩ → |Zρ⟩
- Matrix elements: (Z⊗I)_{αβ} = Z_{ik} δ_{jl}
- Where α = i + (j-1)d and β = k + (l-1)d

# Arguments
- `s::Index`: Site index with dimension d^2 (d=2 for qubits)

# Returns
- `ITensor`: Z operator acting on vectorized density matrix
"""
function make_vectorized_Z_op(s::Index)
    d = 2  # Qubits only
    T = ITensor(prime(s), s)

    # (Z ⊗ I)_{αβ} with α=(i,j), β=(k,l): = Z_{ik} δ_{jl}
    for i in 1:d, j in 1:d, k in 1:d, l in 1:d
        α = i + (j - 1) * d
        β = k + (l - 1) * d
        if j == l
            # Z matrix elements: diag(1, -1)
            z_ik = (i == k) ? (i == 1 ? 1.0 : -1.0) : 0.0
            if z_ik != 0.0
                T[prime(s) => α, s => β] = z_ik
            end
        end
    end
    return T
end

# Register Z operator for "DM" (density matrix) site type
ITensors.op(::OpName"Z", ::SiteType"DM", s::Index) = make_vectorized_Z_op(s)

# ==============================================================================
# Trace bra tensor: ⟨1| for computing Tr(ρ) = ⟨1|ρ⟩
# ==============================================================================

"""
    onebra_site_tensor(s::Index)

Create local ⟨1| bra tensor for trace computation.

Contracts with vec(I), extracting diagonal entries: α = i + (i-1)*d
This gives ⟨1|ρ⟩ = Tr(ρ)

# Arguments
- `s::Index`: Vectorized site index with dimension d^2

# Returns
- `ITensor`: Bra tensor for trace computation
"""
function onebra_site_tensor(s::Index)
    dim_s = dim(s)
    d = Int(round(sqrt(dim_s)))
    
    if d * d != dim_s
        error("Site index dim = $dim_s is not a perfect square; expected d^2 for vectorized |ρ⟩.")
    end

    b = ITensor(dag(s))
    
    # Set diagonal entries only: α = i + (i-1)*d
    for i in 1:d
        α = i + (i - 1) * d
        b[dag(s) => α] = 1.0
    end
    
    return b
end

"""
    make_onebra(psi::MPS)

Create ⟨1| bra for all sites of vectorized density matrix MPS.

# Arguments
- `psi::MPS`: Vectorized density matrix MPS

# Returns
- `Vector{ITensor}`: Array of bra tensors, one per site
"""
make_onebra(psi::MPS) = [onebra_site_tensor(siteinds(psi)[n]) for n in 1:length(psi)]

# ==============================================================================
# Trace computation
# ==============================================================================

"""
    tr_rho(psi::MPS)

Compute trace of vectorized density matrix: Tr(ρ) = ⟨1|ρ⟩

# Arguments
- `psi::MPS`: Vectorized density matrix MPS

# Returns
- `ComplexF64`: Trace value (should be real and ≈1 for normalized ρ)
"""
function tr_rho(psi::MPS)
    onebra = make_onebra(psi)
    L = ITensor(1.0)
    
    for n in 1:length(psi)
        L *= psi[n] * onebra[n]
    end
    
    return scalar(L)
end

# ==============================================================================
# Main correlator function
# ==============================================================================

"""
    correlator(psi::MPS, cor_ops::Tuple, op_sites::Vector{<:Tuple})

Compute Z-operator correlators: Tr(Z_i Z_j ... ρ) = ⟨1|Z_i Z_j ...|ρ⟩

Specialized for Z operators only. Exploits:
- Commutativity: Z_i Z_j = Z_j Z_i, so tuple order doesn't matter
- Idempotency: Z^2 = I, so only parity of multiplicity matters

# Arguments
- `psi::MPS`: Vectorized density matrix MPS
- `cor_ops::Tuple`: Tuple of operator names (must all be "Z")
- `op_sites::Vector{<:Tuple}`: Vector of site index tuples
  Example: [(1,2), (2,3), (1,1)] for ⟨Z_1 Z_2⟩, ⟨Z_2 Z_3⟩, ⟨Z_1 Z_1⟩

# Returns
- `Dict{Tuple{Vararg{Int}},ComplexF64}`: Dictionary mapping site tuples to expectation values

# Example
```julia
z2 = correlator(ρ_vec, ("Z", "Z"), [(1,2), (2,2), (3,1)])
z4 = correlator(ρ_vec, ("Z", "Z", "Z", "Z"), [(1,2,3,4), (1,1,1,1)])
```
"""
function correlator(psi::MPS, cor_ops::Tuple, op_sites::Vector{<:Tuple})
    # Validate that all operators are Z
    if !all(op == "Z" for op in cor_ops)
        error("This specialized correlator only supports Z operators. Got: $cor_ops")
    end
    
    # Handle empty case
    if isempty(op_sites)
        return Dict{Tuple{Vararg{Int}},ComplexF64}(() => tr_rho(psi))
    end

    N = length(cor_ops)
    
    # Validate tuple lengths
    if !all(length(st) == N for st in op_sites)
        error("Every tuple in op_sites must have length = length(cor_ops) = $N.")
    end

    # Canonicalize by sorting (Z operators commute)
    canonical_sites = [tuple(sort(collect(st))...) for st in op_sites]
    unique_sites = unique(canonical_sites)

    # Compute once per unique canonical tuple
    Ccanon = correlator_recursive_compact_zonly(psi, unique_sites)

    # Map back to original requested tuples
    C = Dict{Tuple{Vararg{Int}},ComplexF64}()
    for (orig, can) in zip(op_sites, canonical_sites)
        C[orig] = Ccanon[can]
    end
    
    return C
end

# ==============================================================================
# Recursive batch correlator for sorted site tuples
# ==============================================================================

"""
    correlator_recursive_compact_zonly(psi::MPS, sites::Vector{<:Tuple})

Recursively compute Z-correlators for multiple sorted site tuples.

Uses cached left environment contractions to efficiently evaluate many
correlators in a single pass through the MPS.

# Arguments
- `psi::MPS`: Vectorized density matrix MPS
- `sites::Vector{<:Tuple}`: Sorted site index tuples (lexicographic order)

# Returns
- `Dict{Tuple{Vararg{Int}},ComplexF64}`: Correlator values for each tuple
"""
function correlator_recursive_compact_zonly(psi::MPS, sites::Vector{<:Tuple})
    if isempty(sites)
        return Dict{Tuple{Vararg{Int}},ComplexF64}(() => tr_rho(psi))
    end

    sites = sort(sites)  # Lexicographic sort
    N = length(sites[1])

    # Validate all tuples have same length
    if !all(length(st) == N for st in sites)
        error("All site tuples must have the same length. Got lengths: $(length.(sites))")
    end

    C = Dict{Tuple{Vararg{Int}},ComplexF64}()
    onebra = make_onebra(psi)
    s = siteinds(psi)

    # For each tuple, determine first site and repeat count
    inds_ord = [st[1] for st in sites]
    repeats = [count(==(st[1]), st) - 1 for st in sites]

    op_inds = unique([(inds_ord[idx], repeats[idx]) for idx in eachindex(sites)])
    op_inds = sort(op_inds, by = x -> x[1])

    element = zeros(Int, N)
    L0 = ITensor(1.0)

    add_operator_zonly!(op_inds, sites, L0, 1, element, N, s, psi, onebra, C)
    
    return C
end

# ==============================================================================
# Recursive helper for cached contraction
# ==============================================================================

"""
    add_operator_zonly!(op_inds, sites_ind_prev, L_prev, counter, element, N, s, psi, onebra, C)

Recursive helper for computing Z-correlators with cached environments.

Exploits Z^2 = I: only applies Z if multiplicity is odd.

# Arguments
- `op_inds`: Vector of (site, repeat_count) tuples
- `sites_ind_prev`: Site tuples to evaluate
- `L_prev`: Cached left environment tensor
- `counter::Int`: Current position in operator string
- `element::Vector{Int}`: Current site tuple being built
- `N::Int`: Total number of operators
- `s`: Site indices
- `psi::MPS`: Vectorized density matrix MPS
- `onebra`: Trace bra tensors
- `C::Dict`: Output dictionary (mutated)
"""
function add_operator_zonly!(
    op_inds,
    sites_ind_prev,
    L_prev,
    counter::Int,
    element::Vector{Int},
    N::Int,
    s,
    psi::MPS,
    onebra,
    C::Dict,
)
    for (a, op_info) in enumerate(op_inds)
        site   = op_info[1]
        repeat = op_info[2]

        element[counter:(counter + repeat)] .= site

        # Initialize left environment for first site:
        # Contract sites 1:(site-1) with identity (trace bra only)
        if counter == 1
            L_prev = ITensor(1.0)
            if site > 1
                for n in 1:(site - 1)
                    L_prev = L_prev * psi[n] * onebra[n]
                end
            end
        end

        L = copy(L_prev)

        # Apply local operator: Z only if multiplicity is odd
        # Exploits Z^2 = I property
        op_psi = psi[site]
        if isodd(repeat + 1)
            op_psi = apply(op("Z", s[site]), op_psi)
        end

        L = L * op_psi * onebra[site]

        # Base case: all operator positions consumed
        if counter + repeat == N
            # Contract remaining sites to the right
            if site < length(psi)
                for n in (site + 1):length(psi)
                    L = L * psi[n] * onebra[n]
                end
            end
            
            # Store result
            C[tuple(element...)] = scalar(L)
            
        else
            # Recursive case: process next operator position
            
            # Keep only tuples consistent with current choice
            sites_ind = sites_ind_prev[
                findall(st -> st[counter + repeat] == site, sites_ind_prev)
            ]

            # Remove tuples whose next entry is same site (already handled in repeat)
            deleteat!(
                sites_ind,
                findall(st -> st[counter + repeat + 1] == site, sites_ind)
            )

            # Compute repeat counts for next position
            repeat_next = [
                count(==(st[counter + repeat + 1]), st) - 1
                for st in sites_ind
            ]

            inds_ord_next = [st[counter + repeat + 1] for st in sites_ind]

            op_inds_next = unique([
                (inds_ord_next[idx], repeat_next[idx]) for idx in eachindex(sites_ind)
            ])
            op_inds_next = sort(op_inds_next, by = x -> x[1])

            next_site = op_inds_next[1][1]

            # Contract intermediate sites up to next operator site
            if next_site > site + 1
                for n in (site + 1):(next_site - 1)
                    L = L * psi[n] * onebra[n]
                end
            end

            # Recurse to next operator position
            add_operator_zonly!(
                op_inds_next,
                sites_ind,
                L,
                counter + repeat + 1,
                element,
                N,
                s,
                psi,
                onebra,
                C,
            )
        end

        # Advance cached L_prev for next branch in this loop
        if a < length(op_inds)
            next_site = op_inds[a + 1][1]
            for n in site:(next_site - 1)
                L_prev = L_prev * psi[n] * onebra[n]
            end
        end
    end
end

# ==============================================================================
# High-level functions for Binder parameter calculation
# ==============================================================================

"""
    compute_M2_M4_vectorized(ρ_vec::MPS, L::Int)

Compute M₂ and M₄ for Edwards-Anderson Binder parameter.

    M₂ = (1/L²) ∑ᵢⱼ |⟨ZᵢZⱼ⟩|²
    M₄ = (1/L⁴) ∑ᵢⱼₖₗ |⟨ZᵢZⱼZₖZₗ⟩|²

Uses abs2() to correctly compute |⟨O⟩|² = |Tr(O ρ)|².

# Arguments
- `ρ_vec::MPS`: Vectorized density matrix MPS
- `L::Int`: Number of physical sites

# Returns
- `(M2::Float64, M4::Float64)`: M₂ and M₄ values

# Performance
- Computes L² two-point correlators
- Computes L⁴ four-point correlators
- Total cost: O(L³χ² + L⁵χ²) where χ is bond dimension
"""
function compute_M2_M4_vectorized(ρ_vec::MPS, L::Int)
    # Generate all pairs and quads
    pairs = [(i, j) for i in 1:L for j in 1:L]
    quads = [(i, j, k, l) for i in 1:L for j in 1:L for k in 1:L for l in 1:L]

    # Compute correlators
    z2 = correlator(ρ_vec, ("Z", "Z"), pairs)
    z4 = correlator(ρ_vec, ("Z", "Z", "Z", "Z"), quads)

    # Compute M₂ and M₄ with abs2() for |⟨O⟩|²
    M2 = sum(abs2(z2[p]) for p in pairs) / L^2
    M4 = sum(abs2(z4[q]) for q in quads) / L^4

    return M2, M4
end

"""
    binder_from_trials(M2_trials::AbstractVector, M4_trials::AbstractVector)

Compute Edwards-Anderson Binder parameter from trajectory-resolved observables.

    B_EA = 1 - ⟨M₄⟩ / (3⟨M₂⟩²)

This is the CORRECT way to compute the Binder parameter: average M₂ and M₄
FIRST, then compute the ratio. Do NOT compute B per trajectory then average.

# Arguments
- `M2_trials::AbstractVector`: M₂ values for each trajectory
- `M4_trials::AbstractVector`: M₄ values for each trajectory

# Returns
- `Float64`: Edwards-Anderson Binder parameter

# Notes
Correct:   B = 1 - mean(M₄) / (3 * mean(M₂)²)  ✓
Incorrect: B = mean(1 - M₄ / (3 * M₂²))        ✗
"""
function binder_from_trials(M2_trials::AbstractVector, M4_trials::AbstractVector)
    M2bar = mean(M2_trials)
    M4bar = mean(M4_trials)
    return 1.0 - M4bar / (3.0 * M2bar^2 + eps(Float64))
end

end # module
