"""
Density matrix evolution with MixedStateMPS (full density matrix).

This is the CORRECT implementation for dephasing. Uses full density matrix representation
in the doubled (Choi-Jamiolkowski) MPS form.

Key: MixedStateMPS represents ρ as a "doubled" MPS where indices come in pairs (bra, ket).
Sites 2i-1 = bra, sites 2i = ket. This allows representing general mixed states.

Based on src_1/MPS implementation patterns.
"""

using Random, Statistics
using ITensors, ITensorMPS
using LinearAlgebra

export evolve_density_matrix_one_trial
export ea_binder_density_matrix

# Pauli matrices
const σx = Float64[0 1; 1 0]
const σz = Float64[1 0; 0 -1]

"""
    evolve_density_matrix_one_trial(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)

Evolve density matrix using MixedStateMPS (doubled MPS representation).

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
- `state::MixedStateMPS`: Final density matrix
"""
function evolve_density_matrix_one_trial(L::Int; 
                                         lambda_x::Float64, 
                                         lambda_zz::Float64,
                                         P_x::Float64=0.0, 
                                         P_zz::Float64=0.0,
                                         maxdim::Int=256, 
                                         cutoff::Float64=1e-12,
                                         rng=Random.GLOBAL_RNG)
    # Initialize as |↑↑...↑⟩⟨↑↑...↑| in doubled representation
    # Sites: 1,2 = site 1 (bra, ket), 3,4 = site 2 (bra, ket), ...
    sites = siteinds("Qubit", 2L)
    ρ = MPS(sites, _ -> "Up")  # All up state
    
    T_max = 2L
    
    for t in 1:T_max
        # Step 1: Weak X measurements on all sites
        for i in 1:L
            bra_idx = 2i - 1
            ket_idx = 2i
            
            # Sample measurement outcome based on ⟨X⟩
            # For MixedState: ⟨X⟩ = Tr(ρ·X)
            X_bra = op(σx, sites[bra_idx])
            Xρ = apply(X_bra, copy(ρ); cutoff=cutoff, maxdim=maxdim)
            expval_X = real(inner(ρ, Xρ) / inner(ρ, ρ))
            
            prob_0 = (1 + 2*lambda_x/(1+lambda_x^2)*expval_X) / 2
            prob_0 = clamp(prob_0, 0.0, 1.0)
            outcome = rand(rng) < prob_0 ? 0 : 1
            
            # Apply Kraus operator Π = (I + (-1)^m λ X) / √(2(1+λ²))
            # to BOTH bra and ket
            Π = (I(2) + (-1)^outcome * lambda_x * σx) / sqrt(2*(1+lambda_x^2))
            Π_bra = op(Π, sites[bra_idx])
            Π_ket = op(Π, sites[ket_idx])
            
            ρ = apply([Π_bra, Π_ket], ρ; cutoff=cutoff, maxdim=maxdim)
            normalize!(ρ)
        end
        
        # Step 2: X dephasing on all sites (deterministic channel)
        # ρ → (1-P_x)ρ + P_x·X·ρ·X†
        if P_x > 0
            for i in 1:L
                bra_idx = 2i - 1
                ket_idx = 2i
                
                # Build gate: (1-p)I⊗I + p·X⊗X
                gate = (1-P_x)*op(I(2), sites[bra_idx])*op(I(2), sites[ket_idx]) + 
                       P_x*op(σx, sites[bra_idx])*op(σx, sites[ket_idx])
                
                ρ = apply(gate, ρ; cutoff=cutoff, maxdim=maxdim)
                normalize!(ρ)
            end
        end
        
        # Step 3: Weak ZZ measurements on adjacent bonds
        for i in 1:(L-1)
            j = i + 1
            bra_i = 2i - 1
            ket_i = 2i
            bra_j = 2j - 1
            ket_j = 2j
            
            # Sample based on ⟨ZZ⟩
            Z_bra_i = op(σz, sites[bra_i])
            Z_bra_j = op(σz, sites[bra_j])
            ZZρ = apply([Z_bra_i, Z_bra_j], copy(ρ); cutoff=cutoff, maxdim=maxdim)
            expval_ZZ = real(inner(ρ, ZZρ) / inner(ρ, ρ))
            
            prob_0 = (1 + 2*lambda_zz/(1+lambda_zz^2)*expval_ZZ) / 2
            prob_0 = clamp(prob_0, 0.0, 1.0)
            outcome = rand(rng) < prob_0 ? 0 : 1
            
            # Apply Π to each site (both bra and ket)
            Π = (I(2) + (-1)^outcome * lambda_zz * σz) / sqrt(2*(1+lambda_zz^2))
            gates = [op(Π, sites[bra_i]), op(Π, sites[ket_i]),
                    op(Π, sites[bra_j]), op(Π, sites[ket_j])]
            
            ρ = apply(gates, ρ; cutoff=cutoff, maxdim=maxdim)
            normalize!(ρ)
        end
        
        # Step 4: ZZ dephasing on adjacent bonds
        # ρ → (1-P_zz)ρ + P_zz·(Z⊗Z)·ρ·(Z⊗Z)†
        if P_zz > 0
            for i in 1:(L-1)
                j = i + 1
                bra_i = 2i - 1
                ket_i = 2i
                bra_j = 2j - 1
                ket_j = 2j
                
                # Build ZZ gate on bra and ket
                ZZ_bra = op(σz, sites[bra_i]) * op(σz, sites[bra_j])
                ZZ_ket = op(σz, sites[ket_i]) * op(σz, sites[ket_j])
                II_bra = op(I(2), sites[bra_i]) * op(I(2), sites[bra_j])
                II_ket = op(I(2), sites[ket_i]) * op(I(2), sites[ket_j])
                
                gate = (1-P_zz)*(II_bra*II_ket) + P_zz*(ZZ_bra*ZZ_ket)
                
                ρ = apply(gate, ρ; cutoff=cutoff, maxdim=maxdim)
                normalize!(ρ)
            end
        end
    end
    
    return MixedStateMPS(ρ)
end

"""
    ea_binder_density_matrix(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)

Calculate Edwards-Anderson Binder parameter using MixedStateMPS.

Each trial evolves a full density matrix. Disorder/measurement averaging over trials.
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
        # Evolve density matrix
        state = evolve_density_matrix_one_trial(L; 
                                               lambda_x=lambda_x, 
                                               lambda_zz=lambda_zz,
                                               P_x=P_x, 
                                               P_zz=P_zz,
                                               maxdim=maxdim, 
                                               cutoff=cutoff, 
                                               rng=rng)
        
        # Extract MPS and normalize
        ρ = get_mps(state)
        normalize!(ρ)
        
        # Compute correlators for mixed state
        M2sq, M4sq = compute_correlators_mixed(ρ, L; cutoff=cutoff)
        
        # Compute Binder parameter
        den = 3.0 * max(M2sq^2, 1e-12)
        
        S2s[t] = M2sq
        S4s[t] = M4sq
        Bs[t]  = 1.0 - M4sq / den
    end
    
    # Ensemble averages
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
    compute_correlators_mixed(ρ::MPS, L::Int; cutoff)

Compute correlation functions for MixedStateMPS (doubled representation).

For density matrix ρ in doubled form:
    ⟨Z_i Z_j⟩ = Tr(ρ·Z_i·Z_j)
    
Apply Z operators to bra part only (trace automatically from ket part).
"""
function compute_correlators_mixed(ρ::MPS, L::Int; cutoff=1e-12)
    sites = siteinds(ρ)
    
    # Compute trace for normalization
    tr = inner(ρ, ρ)
    
    # 2-point correlators: ⟨Z_i Z_j⟩ = Tr(ρ Z_i Z_j)
    # In doubled rep: apply Z to bra indices
    sum2_sq = 0.0
    
    for i in 1:L, j in 1:L
        bra_i = 2i - 1
        bra_j = 2j - 1
        
        # Apply Z to bra part
        Z_i = op(σz, sites[bra_i])
        Z_j = op(σz, sites[bra_j])
        
        ρ_temp = apply([Z_i, Z_j], copy(ρ); cutoff=cutoff)
        corr = real(inner(ρ, ρ_temp) / tr)
        sum2_sq += corr^2
    end
    
    # 4-point correlators (sampled for efficiency)
    n_samples = min(1000, L^4)
    sum4_sq = 0.0
    
    for _ in 1:n_samples
        i, j, k, l = rand(1:L, 4)
        
        gates = [op(σz, sites[2*site-1]) for site in [i, j, k, l]]
        ρ_temp = apply(gates, copy(ρ); cutoff=cutoff)
        
        corr = real(inner(ρ, ρ_temp) / tr)
        sum4_sq += corr^2
    end
    
    M2sq = sum2_sq / L^2
    M4sq = (sum4_sq / n_samples) * L^4 / L^4
    
    return M2sq, M4sq
end
