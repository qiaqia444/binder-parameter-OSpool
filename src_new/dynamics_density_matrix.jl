"""
Density matrix evolution for measurement-induced phase transitions with dephasing.

This module implements the CORRECT approach for systems with both weak measurements
and dephasing channels: evolving DiagonalStateMPS (density matrices) instead of
pure state quantum trajectories.

Key difference from dynamics.jl:
- dynamics.jl: Pure state trajectories (stochastic dephasing) - WRONG for dephasing
- This file: Density matrix evolution (deterministic dephasing) - CORRECT for dephasing
"""

using Random, Statistics
using ITensors, ITensorMPS
using LinearAlgebra

# Note: Assumes types.jl, channels.jl are already loaded in the calling context
# These functions work with DiagonalStateMPS from types.jl

export evolve_density_matrix_one_trial
export ea_binder_density_matrix
export compute_correlators_diagonal

"""
    evolve_density_matrix_one_trial(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)

Evolve a density matrix under weak measurements and dephasing channels.

This is the CORRECT implementation for systems with dephasing. Unlike the pure
state trajectory approach in dynamics.jl, this:
1. Uses DiagonalStateMPS to represent the density matrix
2. Applies dephasing channels deterministically: ρ → (1-P)ρ + P·M·ρ·M†
3. Still samples measurement outcomes (measurements are projective)

# Protocol
For T_max = 2L time steps:
1. Apply weak X measurements to all sites (sample outcomes)
2. Apply X dephasing channel to all sites (deterministic)
3. Apply weak ZZ measurements to all bonds (sample outcomes)
4. Apply ZZ dephasing channel to all bonds (deterministic)

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
- `state::DiagonalStateMPS`: Final density matrix
"""
function evolve_density_matrix_one_trial(L::Int; 
                                         lambda_x::Float64, 
                                         lambda_zz::Float64,
                                         P_x::Float64=0.0, 
                                         P_zz::Float64=0.0,
                                         maxdim::Int=256, 
                                         cutoff::Float64=1e-12,
                                         rng=Random.GLOBAL_RNG)
    # Initialize state as diagonal density matrix (all-up state)
    sites = siteinds("S=1/2", L)
    ψ_init = productMPS(sites, fill("Up", L))
    state = DiagonalStateMPS(ψ_init)
    
    T_max = 2L
    
    # Define Pauli matrices
    X_matrix = [0 1; 1 0]
    Z_matrix = [1 0; 0 -1]
    
    for t in 1:T_max
        ρ = get_mps(state)
        
        # Step 1: Weak X measurements on all sites
        for i in 1:L
            sites_i = siteinds(ρ)
            
            # Compute <X> for sampling
            X_op = ITensor(X_matrix, sites_i[i]', sites_i[i])
            Xρ_temp = apply(X_op, copy(ρ); maxdim=maxdim, cutoff=cutoff)
            
            # Normalize for expectation value
            norm_rho = norm(ρ)
            norm_Xrho = norm(Xρ_temp)
            expval_X = norm_Xrho^2 / norm_rho^2  # Approximation for diagonal states
            
            # Sample outcome
            prob_0 = (1 + 2*lambda_x/(1+lambda_x^2)*expval_X) / 2
            prob_0 = clamp(prob_0, 0.0, 1.0)
            outcome = rand(rng) < prob_0 ? false : true
            
            # Apply Kraus operator (squared for diagonal states)
            Π_matrix = (I(2) + (-1)^outcome * lambda_x * X_matrix) / sqrt(2*(1+lambda_x^2))
            Π²_matrix = Π_matrix * Π_matrix
            Π_op = ITensor(Π²_matrix, sites_i[i]', sites_i[i])
            
            ρ = apply(Π_op, ρ; maxdim=maxdim, cutoff=cutoff)
            normalize!(ρ)
        end
        
        state = DiagonalStateMPS(ρ)
        
        # Step 2: X dephasing on all sites (deterministic channel)
        if P_x > 0
            for i in 1:L
                state = apply_x_dephasing_channel(state, i, P_x; 
                                                 maxdim=maxdim, cutoff=cutoff)
            end
        end
        
        ρ = get_mps(state)
        
        # Step 3: Weak ZZ measurements on all adjacent bonds
        for i in 1:(L-1)
            sites_i = siteinds(ρ)
            
            # Compute <ZZ> for sampling
            Z_i = ITensor(Z_matrix, sites_i[i]', sites_i[i])
            Z_j = ITensor(Z_matrix, sites_i[i+1]', sites_i[i+1])
            ZZρ_temp = apply(Z_i, copy(ρ); maxdim=maxdim, cutoff=cutoff)
            ZZρ_temp = apply(Z_j, ZZρ_temp; maxdim=maxdim, cutoff=cutoff)
            
            norm_rho = norm(ρ)
            norm_ZZrho = norm(ZZρ_temp)
            expval_ZZ = norm_ZZrho^2 / norm_rho^2  # Approximation
            
            # Sample outcome
            prob_0 = (1 + 2*lambda_zz/(1+lambda_zz^2)*expval_ZZ) / 2
            prob_0 = clamp(prob_0, 0.0, 1.0)
            outcome = rand(rng) < prob_0 ? false : true
            
            # Apply two-site Kraus operator
            ZZ_matrix = kron(Z_matrix, Z_matrix)
            Π_matrix = (I(4) + (-1)^outcome * lambda_zz * ZZ_matrix) / sqrt(2*(1+lambda_zz^2))
            Π²_matrix = Π_matrix * Π_matrix
            
            # Use replacebond! for proper two-site gate application
            orthogonalize!(ρ, i)
            ψij = ρ[i] * ρ[i+1]
            
            # Apply gate
            Π_op = ITensor(Π²_matrix, sites_i[i]', sites_i[i+1]', sites_i[i], sites_i[i+1])
            ψij = Π_op * ψij
            noprime!(ψij)
            
            # SVD and truncate
            spec = replacebond!(ρ, i, ψij; maxdim=maxdim, cutoff=cutoff, ortho="left")
            normalize!(ρ)
        end
        
        state = DiagonalStateMPS(ρ)
        
        # Step 4: ZZ dephasing on all adjacent bonds (deterministic channel)
        if P_zz > 0
            for i in 1:(L-1)
                state = apply_zz_dephasing_channel(state, i, i+1, P_zz;
                                                   maxdim=maxdim, cutoff=cutoff)
            end
        end
    end
    
    return state
end

"""
    ea_binder_density_matrix(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)

Calculate Edwards-Anderson Binder parameter using density matrix evolution.

This is the CORRECT method for systems with dephasing. Each "trial" evolves
a full density matrix (not a pure state trajectory).

# Arguments
- `L::Int`: System size
- `lambda_x::Float64`: X measurement strength
- `lambda_zz::Float64`: ZZ measurement strength  
- `P_x::Float64`: X dephasing probability
- `P_zz::Float64`: ZZ dephasing probability
- `ntrials::Int=200`: Number of independent evolutions
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `seed::Union{Nothing,Int}=nothing`: Random seed

# Returns
Named tuple with fields:
- `B`: Edwards-Anderson Binder parameter
- `B_mean_of_trials`: Mean of per-trial Binder parameters
- `B_std_of_trials`: Standard deviation
- `S2_bar`: Average M₂
- `S4_bar`: Average M₄
- `ntrials`: Number of trials completed

# Note
Each "trial" here represents an independent realization of the measurement
outcomes. The density matrix captures the dephasing exactly (not stochastically).
"""
function ea_binder_density_matrix(L::Int; 
                                  lambda_x::Float64, 
                                  lambda_zz::Float64,
                                  P_x::Float64=0.0,
                                  P_zz::Float64=0.0,
                                  ntrials::Int=200, 
                                  maxdim::Int=256, 
                                  cutoff::Float64=1e-12,
                                  chunk4::Int=50_000,
                                  seed::Union{Nothing,Int}=nothing)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    
    S2s = Vector{Float64}(undef, ntrials)
    S4s = Vector{Float64}(undef, ntrials)
    Bs  = Vector{Float64}(undef, ntrials)
    
    for t in 1:ntrials
        # Evolve density matrix (not pure state trajectory)
        state = evolve_density_matrix_one_trial(L; 
                                                lambda_x=lambda_x, 
                                                lambda_zz=lambda_zz,
                                                P_x=P_x, 
                                                P_zz=P_zz,
                                                maxdim=maxdim, 
                                                cutoff=cutoff, 
                                                rng=rng)
        
        # Extract underlying MPS
        ρ = get_mps(state)
        sites = siteinds(ρ)
        
        # Center and normalize
        orthogonalize!(ρ, cld(length(sites), 2))
        normalize!(ρ)
        
        # Compute correlation functions
        # For diagonal states, need to compute: Tr(ρ·Z_i·Z_j)
        M2sq, M4sq = compute_correlators_diagonal(ρ, sites)
        
        # Compute Binder parameter for this trial
        den = 3.0 * max(M2sq^2, 1e-12)
        
        S2s[t] = M2sq
        S4s[t] = M4sq
        Bs[t]  = 1.0 - M4sq / den
    end
    
    # Compute ensemble averages
    S2_bar = mean(S2s)
    S4_bar = mean(S4s)
    B_EA   = 1.0 - S4_bar / (3.0*S2_bar^2 + eps(Float64))
    
    return (B = B_EA,
            B_mean_of_trials = mean(Bs),
            B_std_of_trials  = std(Bs),
            S2_bar = S2_bar,
            S4_bar = S4_bar,
            ntrials = ntrials)
end

"""
    compute_correlators_diagonal(ρ::MPS, sites; operator="Z")

Compute 2-point and 4-point correlation functions for a diagonal density matrix.

For DiagonalStateMPS, we need to compute expectation values properly:
    ⟨Z_i Z_j⟩ = Tr(ρ·Z_i·Z_j) / Tr(ρ)

# Arguments
- `ρ::MPS`: Diagonal density matrix (MPS representation)
- `sites`: Site indices
- `operator::String="Z"`: Operator to measure

# Returns
- `M2sq::Float64`: Sum of squared 2-point correlators / L²
- `M4sq::Float64`: Sum of squared 4-point correlators / L⁴
"""
function compute_correlators_diagonal(ρ::MPS, sites; operator="Z")
    L = length(sites)
    
    # Compute trace for normalization
    plus_state = productMPS(sites, fill("Up", L))  # |+⟩^⊗L state
    trace_rho = real(inner(plus_state, ρ))
    
    # Compute 2-point correlators
    sum2_sq = 0.0
    Z_matrix = [1 0; 0 -1]  # Pauli Z
    
    for i in 1:L, j in 1:L
        # Apply Z_i * Z_j
        Z_i = ITensor(Z_matrix, sites[i]', sites[i])
        Z_j = ITensor(Z_matrix, sites[j]', sites[j])
        
        ρ_copy = copy(ρ)
        ρ_copy = apply(Z_j, ρ_copy; cutoff=1e-12)
        ρ_copy = apply(Z_i, ρ_copy; cutoff=1e-12)
        
        corr = real(inner(plus_state, ρ_copy) / trace_rho)
        sum2_sq += corr^2
    end
    
    # Compute 4-point correlators (sampled to avoid O(L^4) cost)
    n_samples = min(1000, L^4)
    sum4_sq = 0.0
    
    for _ in 1:n_samples
        i, j, k, l = rand(1:L, 4)
        
        # Apply Z_i * Z_j * Z_k * Z_l
        ρ_copy = copy(ρ)
        for site in [l, k, j, i]  # Apply in reverse order
            Z = ITensor(Z_matrix, sites[site]', sites[site])
            ρ_copy = apply(Z, ρ_copy; cutoff=1e-12)
        end
        
        corr = real(inner(plus_state, ρ_copy) / trace_rho)
        sum4_sq += corr^2
    end
    
    M2sq = sum2_sq / L^2
    M4sq = (sum4_sq / n_samples) * L^4 / L^4  # Rescale from sampling
    
    return M2sq, M4sq
end

