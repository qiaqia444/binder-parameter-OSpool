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
using Base.Threads

# Import types for MixedStateMPS
include("types.jl")
using .Main: MixedStateMPS, get_mps

export evolve_density_matrix_one_trial
export ea_binder_density_matrix

# Import vectorized correlators module
include("correlators_vectorized.jl")
using .ITensorCorrelators: correlator

# Pauli matrices
const σx = Float64[0 1; 1 0]
const σy = ComplexF64[0 -im; im 0]
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

"""
TRAJECTORY SAMPLING VERSION (COMMENTED OUT - requires lambda_x > 0 for entanglement)
This version samples jump vs no-jump for dephasing channels.

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
    ρ = MPS(sites, _ -> "Up")

    T_max = 50 * L

    # Save state after X noise on final timestep (for Binder calculation)
    ρ_after_X_noise = nothing

    I2 = Matrix{Float64}(I, 2, 2)
    I4 = Matrix{Float64}(I, 4, 4)
    ZZ_2site = kron(σz, σz)

    for t in 1:T_max
        # --------------------------------------------------
        # Step 1: Weak X measurements on all sites
        # --------------------------------------------------
        if lambda_x > 0
            for i in 1:L
                bra_idx = 2i - 1
                ket_idx = 2i

                # Sample measurement outcome from <X> = Tr(ρ X) / Tr(ρ)
                X_bra_state = M_bra(sites, σx, i)
                tr = doubledtrace(ρ)
                expval_X = real(inner(X_bra_state, ρ) / tr)

                prob_0 = (1 + 2 * lambda_x / (1 + lambda_x^2) * expval_X) / 2
                prob_0 = clamp(prob_0, 0.0, 1.0)
                outcome = rand(rng) < prob_0 ? 0 : 1

                # Π_m = (I + (-1)^m λ X) / sqrt(2(1+λ²))
                Π = (I2 + (-1)^outcome * lambda_x * σx) / sqrt(2 * (1 + lambda_x^2))
                Π_bra = op(Π, sites[bra_idx])
                Π_ket = op(Π, sites[ket_idx])

                ρ = apply([Π_bra, Π_ket], ρ; cutoff=cutoff, maxdim=maxdim)

                # Renormalize after sampled measurement
                ρ = ρ / doubledtrace(ρ)
            end
        end

        # --------------------------------------------------
        # Step 2: X dephasing on all sites
        # TRAJECTORY VERSION:
        #   with prob 1-P_x: do nothing
        #   with prob P_x:   ρ -> X ρ X
        # --------------------------------------------------
        if P_x > 0
            for i in 1:L
                if rand(rng) < P_x
                    bra_idx = 2i - 1
                    ket_idx = 2i

                    X_bra = op(σx, sites[bra_idx])
                    X_ket = op(σx, sites[ket_idx])

                    ρ = apply([X_bra, X_ket], ρ; cutoff=cutoff, maxdim=maxdim)

                    # No renormalization needed:
                    # X ρ X is unitary conjugation, so Tr(ρ) is preserved.
                end
            end
        end

        # Save state after X noise on final timestep
        if t == T_max
            ρ_after_X_noise = deepcopy(ρ)
        end

        # --------------------------------------------------
        # Step 3: Weak ZZ measurements on adjacent bonds
        # --------------------------------------------------
        if lambda_zz > 0
            for i in 1:(L - 1)
                j = i + 1
                bra_i = 2i - 1
                ket_i = 2i
                bra_j = 2j - 1
                ket_j = 2j

                # Sample measurement outcome from <Z_i Z_j> = Tr(ρ Z_i Z_j) / Tr(ρ)
                ZZ_bra_state = M_bra(sites, ZZ_2site, i)
                tr = doubledtrace(ρ)
                expval_ZZ = real(inner(ZZ_bra_state, ρ) / tr)

                prob_0 = (1 + 2 * lambda_zz / (1 + lambda_zz^2) * expval_ZZ) / 2
                prob_0 = clamp(prob_0, 0.0, 1.0)
                outcome = rand(rng) < prob_0 ? 0 : 1

                # K_m = (I₄ + (-1)^m λ_zz (Z⊗Z)) / sqrt(2(1+λ_zz²))
                K_m = (I4 + (-1)^outcome * lambda_zz * ZZ_2site) / sqrt(2 * (1 + lambda_zz^2))

                K_bra = op(K_m, sites[bra_i], sites[bra_j])
                K_ket = op(K_m, sites[ket_i], sites[ket_j])

                ρ = apply([K_bra, K_ket], ρ; cutoff=cutoff, maxdim=maxdim)

                # Renormalize after sampled measurement
                ρ = ρ / doubledtrace(ρ)
            end
        end

        # --------------------------------------------------
        # Step 4: ZZ dephasing on adjacent bonds
        # TRAJECTORY VERSION:
        #   with prob 1-P_zz: do nothing
        #   with prob P_zz:   ρ -> (ZZ) ρ (ZZ)
        # --------------------------------------------------
        if P_zz > 0
            for i in 1:(L - 1)
                j = i + 1

                if rand(rng) < P_zz
                    bra_i = 2i - 1
                    ket_i = 2i
                    bra_j = 2j - 1
                    ket_j = 2j

                    ZZ_bra = op(ZZ_2site, sites[bra_i], sites[bra_j])
                    ZZ_ket = op(ZZ_2site, sites[ket_i], sites[ket_j])

                    ρ = apply([ZZ_bra, ZZ_ket], ρ; cutoff=cutoff, maxdim=maxdim)

                    # No renormalization needed:
                    # ZZ ρ ZZ is unitary conjugation, so Tr(ρ) is preserved.
                end
            end
        end
    end

    # Return state after X noise on final timestep
    ρ_return = isnothing(ρ_after_X_noise) ? ρ : ρ_after_X_noise
    return MixedStateMPS(ρ_return), sites
end

"""

"""
CHANNEL AVERAGING VERSION (ACTIVE - works with lambda_x = 0)
This version applies dephasing as linear combination: ρ → (1-P)ρ + P·U·ρ·U†
Creates mixed states even without measurements (lambda_x = 0).
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
    
    T_max = 50 * L  # Time evolution
    
    # Save state after X noise on final timestep (for Binder calculation)
    ρ_after_X_noise = nothing
    
    for t in 1:T_max
        # Step 1: Weak X measurements on all sites
        if lambda_x > 0
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
        end
        
        # Step 2: X dephasing on all sites
        # ρ → (1-P_x)ρ + P_x·X·ρ·X
        # In vectorized form: |ρ⟩ → [(1-p)I + p(X⊗X*)] |ρ⟩
        if P_x > 0
            for i in 1:L
                bra_idx = 2i - 1
                ket_idx = 2i
                
                # Apply CPTP map as linear combination
                ρ_copy = copy(ρ)
                
                # Apply X·ρ·X (conjugation superoperator: X⊗X*)
                X_bra = op(σx, sites[bra_idx])
                X_ket = op(σx, sites[ket_idx])
                ρ_X = apply([X_bra, X_ket], ρ_copy; cutoff=cutoff, maxdim=maxdim)
                
                # Linear combination: (1-p)ρ + p·X·ρ·X
                ρ = (1 - P_x) * ρ_copy + P_x * ρ_X
                # No renormalization - CPTP preserves trace
            end
        end
        
        # Save state after X noise on final timestep (this is the strobe point)
        if t == T_max
            ρ_after_X_noise = deepcopy(ρ)
        end
        
        # Step 3: Weak ZZ measurements on adjacent bonds
        # Use proper 2-site Kraus operator: K_m = (I₄ + (-1)^m λ_zz (Z⊗Z)) / √(2(1+λ_zz²))
        if lambda_zz > 0
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
    
    # Return state after X noise on final timestep (most thermalized strobe point)
    ρ_return = isnothing(ρ_after_X_noise) ? ρ : ρ_after_X_noise
    return MixedStateMPS(ρ_return), sites
end
"""
    ea_binder_density_matrix(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)

Calculate Edwards-Anderson Binder parameter using MixedStateMPS.

Each trial evolves a full density matrix. Uses vectorized correlator for computing
expectation values via <1|O|ρ> = Tr(O·ρ).
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
        
        # Convert MixedStateMPS (doubled) to vectorized density matrix
        ρ_vec = mixed_to_vectorized(ρ, L)
        
        # Compute correlators using vectorized correlator <1|O|ρ>
        M2sq, M4sq = compute_correlators_vectorized(ρ_vec, L)
        
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
    mixed_to_vectorized(ρ_mixed::MPS, L::Int) -> MPS

Convert MixedStateMPS (doubled representation) to vectorized density matrix.

MixedStateMPS has 2L sites with pairs (bra, ket) for each physical site.
Vectorized MPS has L sites, each with dim=d² encoding the density matrix.

For qubits: each physical site (i,j) element → vectorized index α = i + (j-1)*d
"""
function mixed_to_vectorized(ρ_mixed::MPS, L::Int)
    d = 2  # Qubits
    
    # Get site indices from mixed state
    sites_mixed = siteinds(ρ_mixed)
    @assert length(sites_mixed) == 2*L "MixedStateMPS must have 2L sites"
    
    # Create vectorized site indices (dim = d²)
    sites_vec = [Index(d^2, "DM,n=$n") for n in 1:L]
    
    # Build vectorized MPS tensors
    vec_tensors = Vector{ITensor}(undef, L)
    
    # For each physical site, extract density matrix elements
    for n in 1:L
        bra_site = sites_mixed[2n-1]  # Odd index: bra
        ket_site = sites_mixed[2n]    # Even index: ket
        
        # Get the tensor(s) from mixed state for this physical site
        # The MPS structure may have bond indices connecting sites
        T_bra = ρ_mixed[2n-1]
        T_ket = ρ_mixed[2n]
        
        # Get bond indices
        all_inds_bra = inds(T_bra)
        all_inds_ket = inds(T_ket)
        bond_inds_bra = setdiff(all_inds_bra, [bra_site])
        bond_inds_ket = setdiff(all_inds_ket, [ket_site])
        
        # Contract the bra and ket tensors to get the density matrix elements
        # This is complex - we need to extract ρᵢⱼ = ⟨i|ρ|j⟩
        
        # For simplicity, contract the two tensors for this site pair
        T_pair = T_bra * T_ket
        
        # Get remaining bond indices (left and right)
        pair_inds = inds(T_pair)
        site_inds = [bra_site, ket_site]
        bond_inds = setdiff(pair_inds, site_inds)
        
        # Create vectorized tensor
        if length(bond_inds) == 0
            # No bond dimension - product state
            T_vec = ITensor(sites_vec[n])
            for i in 1:d, j in 1:d
                α = i + (j-1)*d
                coeff = T_pair[bra_site => i, ket_site => j]
                T_vec[sites_vec[n] => α] = coeff
            end
        elseif length(bond_inds) == 2
            # Two bond indices (typical case)
            bond_left = bond_inds[1]
            bond_right = bond_inds[2]
            dim_left = dim(bond_left)
            dim_right = dim(bond_right)
            
            T_vec = ITensor(sites_vec[n], bond_left, bond_right)
            
            for i in 1:d, j in 1:d, bl in 1:dim_left, br in 1:dim_right
                α = i + (j-1)*d
                coeff = T_pair[bra_site => i, ket_site => j, bond_left => bl, bond_right => br]
                T_vec[sites_vec[n] => α, bond_left => bl, bond_right => br] = coeff
            end
        else
            # One bond index (first or last site)
            bond = bond_inds[1]
            dim_bond = dim(bond)
            
            T_vec = ITensor(sites_vec[n], bond)
            
            for i in 1:d, j in 1:d, b in 1:dim_bond
                α = i + (j-1)*d
                coeff = T_pair[bra_site => i, ket_site => j, bond => b]
                T_vec[sites_vec[n] => α, bond => b] = coeff
            end
        end
        
        vec_tensors[n] = T_vec
    end
    
    return MPS(vec_tensors)
end

"""
Define vectorized Z operator for d²-dimensional site indices.
"""
function make_vectorized_Z_op(s::Index)
    d = 2  # Qubits
    T = ITensor(s, s')
    # For Z ⊗ I (acting on left): (Z⊗I) |ρ⟩ = |Zρ⟩
    # α = i + (j-1)*d, β = k + (l-1)*d
    # (Z⊗I)_{αβ} with Zᵢₖ δⱼₗ where Z = diag(1,-1)
    for i in 1:d, j in 1:d, k in 1:d, l in 1:d
        α = i + (j-1)*d
        β = k + (l-1)*d
        if j == l  # j must equal l
            # Z matrix elements: diag(1, -1)
            Z_ik = (i == k) ? (i == 1 ? 1.0 : -1.0) : 0.0
            if Z_ik != 0
                T[s=>α, s'=>β] = Z_ik
            end
        end
    end
    return T
end

# Register Z operator for vectorized density matrices
ITensors.op(::OpName"Z", ::SiteType"DM", s::Index) = make_vectorized_Z_op(s)

"""
    compute_correlators_vectorized(ρ_vec::MPS, L::Int) -> (M2sq, M4sq)

Compute M₂ and M₄ using vectorized density matrix correlator function.

Uses correlator from correlators_vectorized.jl which computes <1|O|ρ> for vectorized |ρ⟩.
Parallelized using multi-threading for faster computation.
"""
function compute_correlators_vectorized(ρ_vec::MPS, L::Int)
    # Generate all pairs and quads
    pairs = [(i,j) for i in 1:L for j in 1:L]
    quads = [(i,j,k,l) for i in 1:L for j in 1:L for k in 1:L for l in 1:L]
    
    # Compute correlation functions using vectorized correlator from local module
    # This is the bottleneck - computing L⁴ correlators
    z2 = correlator(ρ_vec, ("Z", "Z"), pairs)
    z4 = correlator(ρ_vec, ("Z", "Z", "Z", "Z"), quads)
    
    # Parallelize summation over threads (small gain but helps)
    sum2_sq = zeros(Float64, nthreads())
    @threads for idx in 1:length(pairs)
        (i,j) = pairs[idx]
        sum2_sq[threadid()] += abs2(z2[(i,j)])
    end
    
    sum4_sq = zeros(Float64, nthreads())
    @threads for idx in 1:length(quads)
        (i,j,k,l) = quads[idx]
        sum4_sq[threadid()] += abs2(z4[(i,j,k,l)])
    end
    
    M2sq = sum(sum2_sq) / L^2
    M4sq = sum(sum4_sq) / L^4
    
    return M2sq, M4sq
end