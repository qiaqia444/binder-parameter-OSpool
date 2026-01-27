"""
Measurement operators and protocols for weak quantum measurements.

This module handles the creation of weak measurement operators and their
probabilistic application to quantum states.
"""

using Random
using ITensors, ITensorMPS, ITensorCorrelators

export create_weak_measurement_operators
export sample_and_apply
export compute_correlators

"""
    create_weak_measurement_operators(sites, lambda_x::Float64, lambda_zz::Float64)

Create weak measurement operators for X and ZZ measurements.

# Arguments
- `sites`: Site indices for the MPS
- `lambda_x::Float64`: Measurement strength for X measurements (0 ≤ λ ≤ 1)
- `lambda_zz::Float64`: Measurement strength for ZZ measurements (0 ≤ λ ≤ 1)

# Returns
- `WEAK_X_0`: Dict mapping site index to X measurement operator (outcome 0)
- `WEAK_X_1`: Dict mapping site index to X measurement operator (outcome 1)
- `WEAK_ZZ_0`: Dict mapping bond to ZZ measurement operator (outcome 0)
- `WEAK_ZZ_1`: Dict mapping bond to ZZ measurement operator (outcome 1)

# Notes
The weak measurement operators are defined as:
    K₀ = (I + λM) / √(2(1 + λ²))
    K₁ = (I - λM) / √(2(1 + λ²))
where M is the measurement observable (X or ZZ).
"""
function create_weak_measurement_operators(sites, lambda_x::Float64, lambda_zz::Float64)
    WEAK_X_0 = Dict{Int,ITensor}()
    WEAK_X_1 = Dict{Int,ITensor}()
    WEAK_ZZ_0 = Dict{Tuple{Int,Int},ITensor}()
    WEAK_ZZ_1 = Dict{Tuple{Int,Int},ITensor}()

    norm_x  = sqrt(2 * (1 + lambda_x^2))
    norm_zz = sqrt(2 * (1 + lambda_zz^2))

    # Create weak X measurement operators
    for i in 1:length(sites)
        Id_i = op("Id", sites[i])
        # Note: ITensor Sx = 0.5 * Pauli_X, so 2*Sx = Pauli_X
        X_i  = 2 * op("Sx", sites[i])
        WEAK_X_0[i] = (Id_i + lambda_x * X_i) / norm_x     # outcome = 0
        WEAK_X_1[i] = (Id_i - lambda_x * X_i) / norm_x     # outcome = 1
    end
    
    # Create weak ZZ measurement operators
    for i in 1:(length(sites)-1)
        # Note: ITensor Sz = 0.5 * Pauli_Z, so 2*Sz = Pauli_Z
        Z_i = 2 * op("Sz", sites[i])
        Z_j = 2 * op("Sz", sites[i+1])
        II  = op("Id", sites[i]) * op("Id", sites[i+1])
        ZZ  = Z_i * Z_j
        WEAK_ZZ_0[(i,i+1)] = (II + lambda_zz * ZZ) / norm_zz   # outcome = 0
        WEAK_ZZ_1[(i,i+1)] = (II - lambda_zz * ZZ) / norm_zz   # outcome = 1
    end
    
    return WEAK_X_0, WEAK_X_1, WEAK_ZZ_0, WEAK_ZZ_1
end

"""
    sample_and_apply(ψ::MPS, K0::ITensor, K1::ITensor, which::Vector{Int}; kwargs...)

Apply a measurement operator by sampling the measurement outcome.

# Arguments
- `ψ::MPS`: Input quantum state
- `K0::ITensor`: Kraus operator for outcome 0
- `K1::ITensor`: Kraus operator for outcome 1
- `which::Vector{Int}`: Sites to apply the measurement to
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `rng`: Random number generator

# Returns
- `ϕ::MPS`: State after measurement with sampled outcome
"""
function sample_and_apply(ψ::MPS, K0::ITensor, K1::ITensor, which::Vector{Int};
                          maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    # Apply both Kraus operators
    ϕ0 = product(K0, ψ, which; maxdim=maxdim, cutoff=cutoff)
    p0 = max(real(inner(dag(ϕ0), ϕ0)), 0.0)
    
    ϕ1 = product(K1, ψ, which; maxdim=maxdim, cutoff=cutoff)
    p1 = max(real(inner(dag(ϕ1), ϕ1)), 0.0)
    
    # Sample outcome based on Born rule probabilities
    tot = p0 + p1
    ϕ = (tot <= 0) ? ϕ0 : (rand(rng) < p0/tot ? ϕ0 : ϕ1)
    
    # Normalize the resulting state
    orthogonalize!(ϕ, which[end])
    normalize!(ϕ)
    return ϕ
end

"""
    compute_correlators(ψ::MPS, sites; operator="Z")

Compute 2-point and 4-point correlation functions.

# Arguments
- `ψ::MPS`: Quantum state
- `sites`: Site indices
- `operator::String="Z"`: Operator to measure (default: "Z")

# Returns
- `M2sq::Float64`: Sum of squared 2-point correlators divided by L²
- `M4sq::Float64`: Sum of squared 4-point correlators divided by L⁴
"""
function compute_correlators(ψ::MPS, sites; operator="Z")
    L_total = length(sites)
    idx = 1:L_total
    n = length(idx)

    # Generate all pairs and quads
    pairs = [(i,j) for i in idx for j in idx]
    quads = [(i,j,k,l) for i in idx for j in idx for k in idx for l in idx]

    # Compute correlation functions
    z2 = correlator(ψ, (operator, operator), pairs)
    z4 = correlator(ψ, (operator, operator, operator, operator), quads)

    # Sum squared correlators
    sum2_sq = 0.0
    @inbounds for (i,j) in pairs
        v = real(z2[(i,j)])
        sum2_sq += v*v
    end
    
    sum4_sq = 0.0
    @inbounds for (i,j,k,l) in quads
        v = real(z4[(i,j,k,l)])
        sum4_sq += v*v
    end

    M2sq = sum2_sq / n^2
    M4sq = sum4_sq / n^4
    
    return M2sq, M4sq
end
