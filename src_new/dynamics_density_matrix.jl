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
    M_bra(sites, M, pos; refs=0)

Create "doubled" MPS for operator M at position pos, identity elsewhere.
From src_1/MPS/tools.jl pattern.
"""
function M_bra(sites::Vector{<:Index}, M::AbstractMatrix, pos::Int; refs=0)
    M_width = Int(log2(size(M)[1]))
    L = length(sites)÷2 - refs
    
    bra = bell(sites)
    # Apply M to ket indices (even: 2,4,6,...)
    bra = apply(op(M, [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]...), bra)
    return bra
end

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
            
            # Sample measurement outcome based on ⟨X⟩ = Tr(ρ·X) / Tr(ρ)
            # Use M_bra / doubledtrace pattern from src_1
            X_bra_state = M_bra(sites, σx, i)
            tr = doubledtrace(ρ)
            expval_X = real(inner(X_bra_state, ρ) / tr)
            
            prob_0 = (1 + 2*lambda_x/(1+lambda_x^2)*expval_X) / 2
            prob_0 = clamp(prob_0, 0.0, 1.0)
            outcome = rand(rng) < prob_0 ? 0 : 1
            
            # Apply Kraus operator Π = (I + (-1)^m λ X) / √(2(1+λ²))
            # to BOTH bra and ket
            Π = (I(2) + (-1)^outcome * lambda_x * σx) / sqrt(2*(1+lambda_x^2))
            Π_bra = op(Π, sites[bra_idx])
            Π_ket = op(Π, sites[ket_idx])
            
            ρ = apply([Π_bra, Π_ket], ρ; cutoff=cutoff, maxdim=maxdim)
            # Renormalize by trace after measurement
            ρ = ρ / doubledtrace(ρ)
        end
        
        # Step 2: X dephasing on all sites (deterministic channel)
        # ρ → (1-P_x)ρ + P_x·X·ρ·X†
        # In vectorized form: |ρ⟩ → [(1-p)I + p(X⊗X*)] |ρ⟩
        # CPTP channel preserves trace - do NOT renormalize
        if P_x > 0
            for i in 1:L
                bra_idx = 2i - 1
                ket_idx = 2i
                
                # Apply CPTP map as linear combination of MPOs
                # Identity superoperator: I⊗I acting on (bra,ket)
                ρ_copy = copy(ρ)
                
                # Conjugation superoperator: X⊗X* = X on ket, X on bra
                X_bra = op(σx, sites[bra_idx])
                X_ket = op(σx, sites[ket_idx])
                ρ_dephased = apply([X_bra, X_ket], ρ; cutoff=cutoff, maxdim=maxdim)
                
                # Linear combination: (1-p)|ρ⟩ + p(X⊗X)|ρ⟩
                ρ = (1-P_x) * ρ_copy + P_x * ρ_dephased
                # No renormalization - CPTP preserves trace
            end
        end
        
        # Step 3: Weak ZZ measurements on adjacent bonds
        # Use proper 2-site Kraus operator: K_m = (I₄ + (-1)^m λ_zz (Z⊗Z)) / √(2(1+λ_zz²))
        for i in 1:(L-1)
            j = i + 1
            bra_i = 2i - 1
            ket_i = 2i
            bra_j = 2j - 1
            ket_j = 2j
            
            # Sample based on ⟨Z_i Z_j⟩ = Tr(ρ·Z_i·Z_j) / Tr(ρ)
            ZZ = kron(σz, σz)
            ZZ_bra_state = M_bra(sites, ZZ, i)
            tr = doubledtrace(ρ)
            expval_ZZ = real(inner(ZZ_bra_state, ρ) / tr)
            
            prob_0 = (1 + 2*lambda_zz/(1+lambda_zz^2)*expval_ZZ) / 2
            prob_0 = clamp(prob_0, 0.0, 1.0)
            outcome = rand(rng) < prob_0 ? 0 : 1
            
            # Build 2-site Kraus operator: K_m = (I₄ + (-1)^m λ_zz (Z⊗Z)) / √(2(1+λ_zz²))
            I4 = Matrix{Float64}(I, 4, 4)
            ZZ_2site = kron(σz, σz)
            K_m = (I4 + (-1)^outcome * lambda_zz * ZZ_2site) / sqrt(2*(1+lambda_zz^2))
            
            # Apply K_m to ket pair (i,j) and K_m to bra pair (since K_m is real)
            K_ket = op(K_m, sites[ket_i], sites[ket_j])
            K_bra = op(K_m, sites[bra_i], sites[bra_j])
            
            ρ = apply([K_bra, K_ket], ρ; cutoff=cutoff, maxdim=maxdim)
            # Renormalize by trace after measurement
            ρ = ρ / doubledtrace(ρ)
        end
        
        # Step 4: ZZ dephasing on adjacent bonds
        # ρ → (1-P_zz)ρ + P_zz·(Z⊗Z)·ρ·(Z⊗Z)†
        # In vectorized form: |ρ⟩ → [(1-p)I + p(ZZ⊗ZZ*)] |ρ⟩
        if P_zz > 0
            for i in 1:(L-1)
                j = i + 1
                bra_i = 2i - 1
                ket_i = 2i
                bra_j = 2j - 1
                ket_j = 2j
                
                # Apply CPTP map as linear combination of MPOs
                ρ_copy = copy(ρ)
                
                # Conjugation superoperator: (ZZ)⊗(ZZ)* = ZZ on ket, ZZ on bra
                ZZ_bra = op(σz, sites[bra_i]) * op(σz, sites[bra_j])
                ZZ_ket = op(σz, sites[ket_i]) * op(σz, sites[ket_j])
                ρ_dephased = apply([ZZ_bra, ZZ_ket], ρ; cutoff=cutoff, maxdim=maxdim)
                
                # Linear combination: (1-p)|ρ⟩ + p(ZZ⊗ZZ)|ρ⟩
                ρ = (1-P_zz) * ρ_copy + P_zz * ρ_dephased
                # No renormalization - CPTP preserves trace
            end
        end
    end
    
    return MixedStateMPS(ρ), sites
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
        state, _ = evolve_density_matrix_one_trial(L; 
                                               lambda_x=lambda_x, 
                                               lambda_zz=lambda_zz,
                                               P_x=P_x, 
                                               P_zz=P_zz,
                                               maxdim=maxdim, 
                                               cutoff=cutoff, 
                                               rng=rng)
        
        # Extract MPS and normalize by trace (not HS norm)
        ρ = get_mps(state)
        ρ = ρ / doubledtrace(ρ)
        
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
    bell(sites)

Create bell state for doubled MPS trace computation.
Uses src_1/MPS pattern: |00⟩+|11⟩ for each (bra,ket) pair.
"""
function id_mps(s1::Index, s2::Index)
    return MPS([1.0; 0.0; 0.0; 1.0], [s1, s2])
end

function bell(sites::Vector{<:Index})
    N = length(sites)
    @assert iseven(N) "Sites must have even length for doubled MPS"
    
    tensors = ITensor[]
    for i in 1:2:N
        a, b = id_mps(sites[i], sites[i+1])
        push!(tensors, a, b)
    end
    return MPS(tensors)
end

"""
    doubledtrace(ρ)

Compute trace of density matrix in doubled MPS representation.
"""
function doubledtrace(ρ::MPS)
    return real(inner(bell(siteinds(ρ)), ρ))
end

"""
    compute_correlators_mixed(ρ, L; cutoff)

Compute correlators for EA Binder parameter.
CORRECTED version based on bug guide.
"""
function corr_ZZ(ρ::MPS, i::Int, j::Int)
    sites = siteinds(ρ)
    tr = doubledtrace(ρ)
    
    if i == j
        # ⟨Z_i Z_i⟩ = ⟨Z²⟩ = ⟨I⟩ = 1
        return 1.0
    elseif abs(i - j) == 1
        # Nearest neighbors: can use M_bra with kron(σz, σz)
        pos = min(i, j)
        M = kron(σz, σz)
        bra = M_bra(sites, M, pos)
        return real(inner(bra, ρ) / tr)
    else
        # Non-nearest neighbors: apply two single-site Z operators
        # ⟨Z_i Z_j⟩ = Tr(ρ · Z_i · Z_j) = inner(bell, Z_i · Z_j · ρ)
        # Apply Z to ket indices (even)
        bra = bell(sites)
        Z_i = op(σz, sites[2*i])
        Z_j = op(σz, sites[2*j])
        bra = apply([Z_i, Z_j], bra; cutoff=1e-14)
        return real(inner(bra, ρ) / tr)
    end
end

function compute_correlators_mixed(ρ::MPS, L::Int; cutoff=1e-12)
    sites = siteinds(ρ)
    tr = doubledtrace(ρ)
    
    # 2-point correlators: compute all C_ij = ⟨Z_i Z_j⟩
    sum2_sq = 0.0
    for i in 1:L, j in 1:L
        C_ij = corr_ZZ(ρ, i, j)
        sum2_sq += C_ij^2
    end
    
    # 4-point correlators: ⟨Z_i Z_j Z_k Z_l⟩
    # Always apply all four Z operators sequentially (duplicates cancel: Z²=I)
    n_samples = min(1000, L^4)
    sum4_sq = 0.0
    
    for _ in 1:n_samples
        indices = rand(1:L, 4)
        
        # Apply four Z operators sequentially to ket legs
        bra = bell(sites)
        for idx in indices
            Z_gate = op(σz, sites[2*idx])
            bra = apply(Z_gate, bra; cutoff=1e-14)
        end
        corr = real(inner(bra, ρ) / tr)
        sum4_sq += corr^2
    end
    
    M2sq = sum2_sq / L^2
    M4sq = (sum4_sq / n_samples) * L^4 / L^4
    
    return M2sq, M4sq
end
