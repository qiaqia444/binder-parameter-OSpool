"""
Edwards-Anderson Binder parameter calculations.

This module computes the EA Binder parameter using Monte Carlo sampling
over quantum trajectories.
"""

using Random, Statistics
using ITensors, ITensorMPS

export ea_binder_mc, ea_binder_mc_dephasing

# Import required functions
include("dynamics.jl")
using .Main: evolve_one_trial, evolve_one_trial_dephasing

include("measurements.jl")
using .Main: compute_correlators

"""
    ea_binder_mc(L::Int; lambda_x, lambda_zz, ntrials=200, kwargs...)

Calculate Edwards-Anderson Binder parameter via Monte Carlo sampling.

The Binder parameter is defined as:
    B_EA = 1 - ⟨M₄⟩ / (3⟨M₂⟩²)
where:
    M₂ = (1/L²) Σᵢⱼ |⟨ZᵢZⱼ⟩|²
    M₄ = (1/L⁴) Σᵢⱼₖₗ |⟨ZᵢZⱼZₖZₗ⟩|²

# Arguments
- `L::Int`: System size
- `lambda_x::Float64`: X measurement strength
- `lambda_zz::Float64`: ZZ measurement strength
- `ntrials::Int=200`: Number of Monte Carlo trajectories
- `maxdim::Int=256`: Maximum bond dimension
- `cutoff::Float64=1e-12`: Truncation cutoff
- `chunk4::Int=50_000`: Chunk size for 4-point correlators (unused, for compatibility)
- `seed::Union{Nothing,Int}=nothing`: Random seed

# Returns
Named tuple with fields:
- `B`: Edwards-Anderson Binder parameter
- `B_mean_of_trials`: Mean of per-trajectory Binder parameters
- `B_std_of_trials`: Standard deviation of per-trajectory Binder parameters
- `S2_bar`: Average M₂
- `S4_bar`: Average M₄
- `ntrials`: Number of trials completed
"""
function ea_binder_mc(L::Int; lambda_x::Float64, lambda_zz::Float64,
                      ntrials::Int=200, maxdim::Int=256, cutoff::Float64=1e-12,
                      chunk4::Int=50_000, seed::Union{Nothing,Int}=nothing)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    S2s = Vector{Float64}(undef, ntrials)
    S4s = Vector{Float64}(undef, ntrials)
    Bs  = Vector{Float64}(undef, ntrials)

    for t in 1:ntrials
        # Evolve one quantum trajectory
        ψ, sites = evolve_one_trial(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                    maxdim=maxdim, cutoff=cutoff, rng=rng)
        
        # Center and normalize the state
        orthogonalize!(ψ, cld(length(sites),2))
        normalize!(ψ)

        # Compute correlation functions
        M2sq, M4sq = compute_correlators(ψ, sites; operator="Z")

        # Compute Binder parameter for this trajectory
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
    ea_binder_mc_dephasing(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)

Calculate Edwards-Anderson Binder parameter with dephasing channels.

Same as `ea_binder_mc` but includes environmental dephasing noise.

# Additional Arguments
- `P_x::Float64`: X dephasing probability
- `P_zz::Float64`: ZZ dephasing probability

# Returns
Same named tuple as `ea_binder_mc`
"""
function ea_binder_mc_dephasing(L::Int; lambda_x::Float64, lambda_zz::Float64,
                                P_x::Float64, P_zz::Float64,
                                ntrials::Int=1000, maxdim::Int=256, cutoff::Float64=1e-12,
                                chunk4::Int=50_000, seed::Union{Nothing,Int}=nothing)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    S2s = Vector{Float64}(undef, ntrials)
    S4s = Vector{Float64}(undef, ntrials)
    Bs  = Vector{Float64}(undef, ntrials)

    for t in 1:ntrials
        # Evolve one quantum trajectory with dephasing
        ψ, sites = evolve_one_trial_dephasing(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                              P_x=P_x, P_zz=P_zz,
                                              maxdim=maxdim, cutoff=cutoff, rng=rng)
        
        # Center and normalize the state
        orthogonalize!(ψ, cld(length(sites),2))
        normalize!(ψ)

        # Compute correlation functions
        M2sq, M4sq = compute_correlators(ψ, sites; operator="Z")

        # Compute Binder parameter for this trajectory
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
