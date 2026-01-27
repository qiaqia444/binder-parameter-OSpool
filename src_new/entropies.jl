"""
Entanglement entropy calculations for quantum many-body systems.

This module provides functions for computing bipartite entanglement entropy
and related information-theoretic quantities.
"""

using Random, Statistics
using ITensors, ITensorMPS

export calculate_bipartite_entropy, weak_bipartite_entropy

"""
    calculate_bipartite_entropy(ψ::MPS, cut::Int) -> Float64

Calculate bipartite entanglement entropy across a cut.

The von Neumann entropy is computed as:
    S = -Σᵢ λᵢ² log₂(λᵢ²)
where λᵢ are the Schmidt values across the cut.

# Arguments
- `ψ::MPS`: Quantum state
- `cut::Int`: Position of the bipartite cut (between sites cut and cut+1)

# Returns
- `entropy::Float64`: von Neumann entanglement entropy in bits
"""
function calculate_bipartite_entropy(ψ::MPS, cut::Int)
    L = length(ψ)
    if cut <= 0 || cut >= L
        return 0.0
    end
    
    try
        # Create copy and move orthogonality center to cut position
        ψ_copy = copy(ψ)
        orthogonalize!(ψ_copy, cut)
        
        # Get bond between sites cut and cut+1
        if cut < L
            bond_idx = commonind(ψ_copy[cut], ψ_copy[cut+1])
            if bond_idx !== nothing
                # Perform SVD: M = USV†
                U, S, V = svd(ψ_copy[cut], (bond_idx,))
                
                # Extract Schmidt values (singular values)
                schmidt_vals = Float64[]
                for i in 1:dim(S, 1)
                    val = S[i,i]
                    if abs(val) > 1e-12
                        push!(schmidt_vals, abs(val))
                    end
                end
                
                if isempty(schmidt_vals)
                    return 0.0
                end
                
                # Normalize Schmidt values
                schmidt_vals = schmidt_vals ./ norm(schmidt_vals)
                
                # Calculate von Neumann entropy: S = -Σ λᵢ² log₂(λᵢ²)
                entropy = -sum(s^2 * log2(s^2 + 1e-16) for s in schmidt_vals)
                return entropy
            end
        end
        
        return 0.0
    catch e
        @warn "Error calculating entropy at cut $cut: $e"
        return 0.0
    end
end

"""
    weak_bipartite_entropy(L::Int, lambda_x, lambda_zz, P_x, P_zz; kwargs...)

Compute bipartite entanglement entropy for weak measurement protocol.

This function evolves multiple quantum trajectories and computes the average
bipartite entanglement entropy across the center of the chain.

# Arguments
- `L::Int`: System size
- `lambda_x::Float64`: X measurement strength
- `lambda_zz::Float64`: ZZ measurement strength
- `P_x::Float64`: X dephasing probability
- `P_zz::Float64`: ZZ dephasing probability
- `ntrials::Int=100`: Number of trajectories to average over
- `maxdim::Int=128`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `seed::Union{Nothing,Int}=nothing`: Random seed

# Returns
Named tuple with fields:
- `entropy_mean`: Mean bipartite entropy
- `entropy_std`: Standard deviation of entropy
- `entropy_sem`: Standard error of the mean
- `ntrials`: Number of trials completed
"""
function weak_bipartite_entropy(L::Int; lambda_x::Float64, lambda_zz::Float64,
                                P_x::Float64=0.0, P_zz::Float64=0.0,
                                ntrials::Int=100, maxdim::Int=128, cutoff::Float64=1e-12,
                                seed::Union{Nothing,Int}=nothing)
    # Import dynamics function
    include("dynamics.jl")
    using .Main: evolve_one_trial_dephasing
    
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    entropies = Vector{Float64}(undef, ntrials)
    
    # Position of bipartite cut (middle of chain)
    cut_position = div(L, 2)
    
    for t in 1:ntrials
        # Evolve one quantum trajectory
        if P_x == 0.0 && P_zz == 0.0
            # Use standard evolution without dephasing
            include("dynamics.jl")
            using .Main: evolve_one_trial
            ψ, sites = evolve_one_trial(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                        maxdim=maxdim, cutoff=cutoff, rng=rng)
        else
            # Use evolution with dephasing
            ψ, sites = evolve_one_trial_dephasing(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                                  P_x=P_x, P_zz=P_zz,
                                                  maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
        
        # Compute bipartite entropy
        entropy = calculate_bipartite_entropy(ψ, cut_position)
        entropies[t] = entropy
    end
    
    return (entropy_mean = mean(entropies),
            entropy_std = std(entropies),
            entropy_sem = std(entropies) / sqrt(ntrials),
            ntrials = ntrials)
end
