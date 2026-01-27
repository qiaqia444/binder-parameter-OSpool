"""
Time evolution and quantum channel operations.

This module implements quantum dynamics including weak measurements,
dephasing channels, and unitary evolution.
"""

using Random
using ITensors, ITensorMPS

export evolve_one_trial, evolve_one_trial_dephasing
export apply_single_site_gate, apply_two_site_gate
export apply_x_dephasing, apply_zz_dephasing

# Import measurement functions
include("measurements.jl")
using .Main: create_weak_measurement_operators, sample_and_apply

# Import state constructors
include("state_constructors.jl")
using .Main: create_up_state_mps

"""
    apply_single_site_gate(ψ::MPS, site::Int, gate::ITensor; kwargs...)

Apply a single-site gate to an MPS.

# Arguments
- `ψ::MPS`: Input quantum state
- `site::Int`: Site to apply gate to
- `gate::ITensor`: Gate operator
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff

# Returns
- `ϕ::MPS`: State after gate application
"""
function apply_single_site_gate(ψ::MPS, site::Int, gate::ITensor;
                                maxdim::Int=256, cutoff::Float64=1e-12)
    ϕ = copy(ψ)
    s = siteind(ϕ, site)
    ϕ[site] = gate * ϕ[site]
    ϕ[site] = replaceinds(ϕ[site], s' => s)
    orthogonalize!(ϕ, site)
    normalize!(ϕ)
    return ϕ
end

"""
    apply_two_site_gate(ψ::MPS, i::Int, j::Int, gate::ITensor; kwargs...)

Apply a two-site gate to an MPS (sites must be adjacent).

# Arguments
- `ψ::MPS`: Input quantum state
- `i::Int`: First site
- `j::Int`: Second site (must be adjacent to i)
- `gate::ITensor`: Two-site gate operator
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff

# Returns
- `ϕ::MPS`: State after gate application
"""
function apply_two_site_gate(ψ::MPS, i::Int, j::Int, gate::ITensor;
                             maxdim::Int=256, cutoff::Float64=1e-12)
    @assert abs(i-j) == 1 "Sites must be adjacent"
    left_site = min(i, j)
    right_site = max(i, j)
    
    ϕ = copy(ψ)
    orthogonalize!(ϕ, left_site)
    
    s_left = siteind(ϕ, left_site)
    s_right = siteind(ϕ, right_site)
    
    # Contract two-site tensor
    wf = ϕ[left_site] * ϕ[right_site]
    # Apply gate
    wf = gate * wf
    # Fix indices
    wf = replaceinds(wf, [s_left', s_right'] => [s_left, s_right])
    
    # SVD decomposition
    U, S, V = svd(wf, uniqueinds(ϕ[left_site], ϕ[right_site]); 
                  maxdim=maxdim, cutoff=cutoff)
    ϕ[left_site] = U
    ϕ[right_site] = S * V
    
    orthogonalize!(ϕ, right_site)
    normalize!(ϕ)
    return ϕ
end

"""
    apply_x_dephasing(ψ::MPS, site::Int, P_x::Float64, sites; kwargs...)

Apply X dephasing channel to a single site.

The X dephasing channel is defined as:
    ρ → (1-P_x)ρ + P_x·X·ρ·X

For pure state trajectories, this is implemented by:
- With probability P_x: apply X (flip the spin)
- With probability (1-P_x): do nothing

# Arguments
- `ψ::MPS`: Input quantum state
- `site::Int`: Site to apply dephasing to
- `P_x::Float64`: Dephasing probability (0 ≤ P ≤ 1)
- `sites`: Site indices
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `rng`: Random number generator

# Returns
- `ψ::MPS`: State after dephasing channel
"""
function apply_x_dephasing(ψ::MPS, site::Int, P_x::Float64, sites; 
                           maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    P_x <= 0 && return ψ
    
    if rand(rng) < P_x
        # Apply Kraus operator K₁ = √P_x · X
        X_gate = 2 * op("Sx", sites[site])  # Pauli X = 2*Sx
        ψ_new = apply_single_site_gate(ψ, site, X_gate; maxdim=maxdim, cutoff=cutoff)
        normalize!(ψ_new)
        return ψ_new
    else
        # Apply Kraus operator K₀ = √(1-P_x) · I (identity)
        return ψ
    end
end

"""
    apply_zz_dephasing(ψ::MPS, i::Int, j::Int, P_zz::Float64, sites; kwargs...)

Apply ZZ dephasing channel to adjacent sites.

The ZZ dephasing channel is defined as:
    ρ → (1-P_zz)ρ + P_zz·(Z⊗Z)·ρ·(Z⊗Z)

For pure state trajectories, this is implemented by:
- With probability P_zz: apply Z⊗Z
- With probability (1-P_zz): do nothing

# Arguments
- `ψ::MPS`: Input quantum state
- `i::Int`: First site
- `j::Int`: Second site (must be adjacent)
- `P_zz::Float64`: Dephasing probability (0 ≤ P ≤ 1)
- `sites`: Site indices
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `rng`: Random number generator

# Returns
- `ψ::MPS`: State after dephasing channel
"""
function apply_zz_dephasing(ψ::MPS, i::Int, j::Int, P_zz::Float64, sites;
                            maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    P_zz <= 0 && return ψ
    if abs(i-j) != 1
        return ψ
    end
    
    if rand(rng) < P_zz
        # Apply Kraus operator K₁ = √P_zz · (Z⊗Z)
        Z_i = 2 * op("Sz", sites[i])  # Pauli Z = 2*Sz
        Z_j = 2 * op("Sz", sites[j])
        ZZ = Z_i * Z_j
        ψ_new = apply_two_site_gate(ψ, i, j, ZZ; maxdim=maxdim, cutoff=cutoff)
        normalize!(ψ_new)
        return ψ_new
    else
        # Apply Kraus operator K₀ = √(1-P_zz) · I (identity)
        return ψ
    end
end

"""
    evolve_one_trial(L::Int; lambda_x, lambda_zz, maxdim=256, cutoff=1e-12, rng)

Evolve a single quantum trajectory with weak measurements only (no dephasing).

The protocol consists of T_max = 2L time steps, where at each step:
1. Apply weak X measurements to all sites
2. Apply weak ZZ measurements to all adjacent bonds

# Arguments
- `L::Int`: System size
- `lambda_x::Float64`: X measurement strength
- `lambda_zz::Float64`: ZZ measurement strength
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `rng`: Random number generator

# Returns
- `ψ::MPS`: Final quantum state after evolution
- `sites`: Site indices
"""
function evolve_one_trial(L::Int; lambda_x::Float64, lambda_zz::Float64,
                          maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    ψ, sites = create_up_state_mps(L)
    KX0, KX1, KZZ0, KZZ1 = create_weak_measurement_operators(sites, lambda_x, lambda_zz)
    T_max = 2L
    
    for _ in 1:T_max
        # Weak X measurements on all sites
        for i in 1:L
            ψ = sample_and_apply(ψ, KX0[i], KX1[i], [i]; maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
        
        # Weak ZZ measurements on all adjacent bonds
        for i in 1:(L-1)
            ψ = sample_and_apply(ψ, KZZ0[(i,i+1)], KZZ1[(i,i+1)], [i,i+1];
                                 maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
    end
    
    return ψ, sites
end

"""
    evolve_one_trial_dephasing(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)

Evolve a single quantum trajectory with weak measurements and dephasing channels.

The protocol consists of T_max = 2L time steps, where at each step:
1. Apply weak X measurements to all sites
2. Apply X dephasing to all sites
3. Apply weak ZZ measurements to all adjacent bonds
4. Apply ZZ dephasing to all adjacent bonds

# Arguments
- `L::Int`: System size
- `lambda_x::Float64`: X measurement strength
- `lambda_zz::Float64`: ZZ measurement strength
- `P_x::Float64`: X dephasing probability
- `P_zz::Float64`: ZZ dephasing probability
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `rng`: Random number generator

# Returns
- `ψ::MPS`: Final quantum state after evolution
- `sites`: Site indices
"""
function evolve_one_trial_dephasing(L::Int; lambda_x::Float64, lambda_zz::Float64,
                                   P_x::Float64, P_zz::Float64,
                                   maxdim::Int=256, cutoff::Float64=1e-12, 
                                   rng=Random.GLOBAL_RNG)
    ψ, sites = create_up_state_mps(L)
    KX0, KX1, KZZ0, KZZ1 = create_weak_measurement_operators(sites, lambda_x, lambda_zz)
    T_max = 2L
    
    for _ in 1:T_max
        # First: Apply weak X measurements on all sites
        for i in 1:L
            ψ = sample_and_apply(ψ, KX0[i], KX1[i], [i]; maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
        # Second: Apply X dephasing on all sites
        for i in 1:L
            ψ = apply_x_dephasing(ψ, i, P_x, sites; maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
        
        # Third: Apply weak ZZ measurements on all adjacent bonds
        for i in 1:(L-1)
            ψ = sample_and_apply(ψ, KZZ0[(i,i+1)], KZZ1[(i,i+1)], [i,i+1];
                                 maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
        # Fourth: Apply ZZ dephasing on all adjacent bonds
        for i in 1:(L-1)
            ψ = apply_zz_dephasing(ψ, i, i+1, P_zz, sites; maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
    end
    
    return ψ, siteinds(ψ)
end
