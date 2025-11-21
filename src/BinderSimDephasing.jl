module BinderSimDephasing

using Random, Statistics
using ITensors, ITensorMPS, ITensorCorrelators

export ea_binder_mc_dephasing, evolve_one_trial_dephasing
export apply_x_dephasing, apply_zz_dephasing
export apply_single_site_gate, apply_two_site_gate

# Import shared functions from BinderSim
include("BinderSim.jl")
using .BinderSim: create_up_state_mps, create_weak_measurement_operators, sample_and_apply

# ========================================
# Gate Application Functions
# ========================================

function apply_single_site_gate(ψ::MPS, site::Int, gate::ITensor;
                                maxdim::Int=256, cutoff::Float64=1e-12)
    """Apply single-site gate to MPS"""
    ϕ = copy(ψ)
    s = siteind(ϕ, site)
    ϕ[site] = gate * ϕ[site]
    ϕ[site] = replaceinds(ϕ[site], s' => s)
    orthogonalize!(ϕ, site)
    normalize!(ϕ)
    return ϕ
end

function apply_two_site_gate(ψ::MPS, i::Int, j::Int, gate::ITensor;
                             maxdim::Int=256, cutoff::Float64=1e-12)
    """Apply two-site gate to MPS"""
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

# ========================================
# Dephasing Channel Functions
# ========================================

function apply_x_dephasing(ψ::MPS, site::Int, P_x::Float64, sites; 
                           maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    """
    Apply X dephasing channel: ρ → (1-P_x)ρ + P_x·X·ρ·X
    
    For pure state trajectories:
    - With probability (1-P_x): apply √(1-P_x)·I (keep state)
    - With probability P_x: apply √P_x·X (flip state)
    
    This samples from the quantum channel mixture.
    """
    P_x <= 0 && return ψ
    
    # Sample which Kraus operator to apply
    if rand(rng) < P_x
        # Apply Kraus operator K₁ = √P_x · X
        X_gate = 2 * op("Sx", sites[site])  # Pauli X = 2*Sx
        ψ_new = apply_single_site_gate(ψ, site, X_gate; maxdim=maxdim, cutoff=cutoff)
        # Renormalize (the √P_x factor is handled by the probabilistic sampling)
        normalize!(ψ_new)
        return ψ_new
    else
        # Apply Kraus operator K₀ = √(1-P_x) · I (identity)
        return ψ
    end
end

function apply_zz_dephasing(ψ::MPS, i::Int, j::Int, P_zz::Float64, sites;
                            maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    """
    Apply ZZ dephasing channel: ρ → (1-P_zz)ρ + P_zz·(Z⊗Z)·ρ·(Z⊗Z)
    
    For pure state trajectories:
    - With probability (1-P_zz): apply √(1-P_zz)·I (keep state)
    - With probability P_zz: apply √P_zz·(Z⊗Z) (apply ZZ)
    
    This samples from the quantum channel mixture.
    """
    P_zz <= 0 && return ψ
    if abs(i-j) != 1
        return ψ
    end
    
    # Sample which Kraus operator to apply
    if rand(rng) < P_zz
        # Apply Kraus operator K₁ = √P_zz · (Z⊗Z)
        Z_i = 2 * op("Sz", sites[i])  # Pauli Z = 2*Sz
        Z_j = 2 * op("Sz", sites[j])
        ZZ = Z_i * Z_j
        ψ_new = apply_two_site_gate(ψ, i, j, ZZ; maxdim=maxdim, cutoff=cutoff)
        # Renormalize (the √P_zz factor is handled by the probabilistic sampling)
        normalize!(ψ_new)
        return ψ_new
    else
        # Apply Kraus operator K₀ = √(1-P_zz) · I (identity)
        return ψ
    end
end

# ========================================
# Evolution and Binder Parameter Calculation
# ========================================

function evolve_one_trial_dephasing(L::Int; lambda_x::Float64, lambda_zz::Float64,
                                   P_x::Float64, P_zz::Float64,
                                   maxdim::Int=256, cutoff::Float64=1e-12, 
                                   rng=Random.GLOBAL_RNG)
    """Evolve one trial with weak measurements + dephasing channels, T_max = 2L"""
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
    
    # Return MPS and the original sites (not modified by gate applications)
    return ψ, siteinds(ψ)
end

function ea_binder_mc_dephasing(L::Int; lambda_x::Float64, lambda_zz::Float64,
                                P_x::Float64, P_zz::Float64,
                                ntrials::Int=1000, maxdim::Int=256, cutoff::Float64=1e-12,
                                chunk4::Int=50_000, seed::Union{Nothing,Int}=nothing)
    """
    Calculate Binder parameter with dephasing channels
    Uses same correlation measurement approach as standard simulations
    """
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    S2s = Vector{Float64}(undef, ntrials)
    S4s = Vector{Float64}(undef, ntrials)
    Bs  = Vector{Float64}(undef, ntrials)

    for t in 1:ntrials
        ψ, sites = evolve_one_trial_dephasing(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                              P_x=P_x, P_zz=P_zz,
                                              maxdim=maxdim, cutoff=cutoff, rng=rng)
        
        # Center and normalize
        orthogonalize!(ψ, cld(length(sites),2))
        normalize!(ψ)

        # Use ALL sites for correlation calculation
        L_total = length(sites)
        idx = 1:L_total
        n = length(idx)

        pairs = [(i,j) for i in idx for j in idx]
        quads = [(i,j,k,l) for i in idx for j in idx for k in idx for l in idx]

        z2 = correlator(ψ, ("Z","Z"), pairs)
        z4 = correlator(ψ, ("Z","Z","Z","Z"), quads)

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
        den  = 3.0 * max(M2sq^2, 1e-12)
        
        S2s[t] = M2sq
        S4s[t] = M4sq
        Bs[t]  = 1.0 - M4sq / den
    end

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

end
