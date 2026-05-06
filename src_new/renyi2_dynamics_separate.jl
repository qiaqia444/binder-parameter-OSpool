# ============================================================
# SEPARATE Rényi-2 DYNAMICS PIPELINE
#
# This file provides INDEPENDENT evolution functions for Rényi-2 Binder.
# It does NOT modify the original EA Binder or old dynamics functions.
#
# Key difference from original:
#  - Independent strobe options for intermediate measurement recording
#  - Otherwise identical physics to evolve_density_matrix_one_trial_new
# ============================================================

using Random
using Statistics
using ITensors, ITensorMPS
using LinearAlgebra

# Import types
include("types.jl")
using .Main: MixedStateMPS, get_mps

export evolve_renyi2_density_matrix_one_trial_separate
export renyi2_binder_density_matrix_separate

# Pauli matrices
const σx = Float64[0 1; 1 0]
const σz = Float64[1 0; 0 -1]
const I2 = Matrix{Float64}(I, 2, 2)
const I4 = Matrix{Float64}(I, 4, 4)
const ZZ_2site = kron(σz, σz)


"""
    M_bra(sites, M, pos; refs=0)

Create "doubled" MPS for operator M at position pos, identity elsewhere.
From original dynamics_density_matrix.jl.
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
    evolve_renyi2_density_matrix_one_trial_separate(
        L::Int;
        lambda_x::Float64,
        lambda_zz::Float64,
        P_x::Float64=0.0,
        P_zz::Float64=0.0,
        maxdim::Int=256,
        cutoff::Float64=1e-12,
        rng::AbstractRNG=MersenneTwister(),
        strobe::Symbol=:after_full_layer,
    )

Separate Rényi-2 evolution (identical physics to evolve_density_matrix_one_trial_new).

Strobe options (new, only for this pipeline):
  :after_x_measurement    - Record state after X measurements
  :after_x_noise          - Record state after X dephasing
  :after_zz_measurement   - Record state after ZZ measurements
  :after_full_layer       - Record state after full evolution (default)

Returns:
  (state, strobe_record)
  state: MixedStateMPS at final time
  strobe_record: Dict of intermediate states
"""
function evolve_renyi2_density_matrix_one_trial_separate(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    P_x::Float64=0.0,
    P_zz::Float64=0.0,
    maxdim::Int=256,
    cutoff::Float64=1e-12,
    rng::AbstractRNG=MersenneTwister(),
    strobe::Symbol=:after_full_layer,
)
    @assert strobe in (:after_x_measurement, :after_x_noise, :after_zz_measurement, :after_full_layer) "Invalid strobe option"
    
    # Initialize as |↑↑...↑⟩⟨↑↑...↑| in doubled representation
    sites = siteinds("Qubit", 2L)
    ρ = MPS(sites, _ -> "Up")  # All up state
    
    T_max = 2 * L  # Time evolution
    
    strobe_record = Dict{Symbol, MixedStateMPS}()

    for t in 1:T_max
        # Step 1: Weak X measurements on all sites
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
        
        if strobe === :after_x_measurement
            strobe_record[:after_x_measurement] = MixedStateMPS(deepcopy(ρ))
        end

        # Step 2: X dephasing on all sites
        # Channel averaging: ρ → (1-P_x)ρ + P_x·X·ρ·X
        if P_x > 0
            for i in 1:L
                bra_idx = 2i - 1
                ket_idx = 2i
                
                # Copy for unperturbed part
                ρ_copy = copy(ρ)
                
                # Apply X to both bra and ket: X·ρ·X
                X_bra = op(σx, sites[bra_idx])
                X_ket = op(σx, sites[ket_idx])
                ρ_X = apply([X_bra, X_ket], ρ_copy; cutoff=cutoff, maxdim=maxdim)
                
                # Linear combination: (1-P_x)ρ + P_x·X·ρ·X
                ρ = (1 - P_x) * ρ + P_x * ρ_X
                # No renormalization - CPTP preserves trace
            end
        end
        
        if strobe === :after_x_noise
            strobe_record[:after_x_noise] = MixedStateMPS(deepcopy(ρ))
        end

        # Step 3: Weak ZZ measurements on adjacent bonds
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
        
        if strobe === :after_zz_measurement
            strobe_record[:after_zz_measurement] = MixedStateMPS(deepcopy(ρ))
        end

        # Step 4: ZZ dephasing on adjacent bonds
        # Channel averaging: ρ → (1-P_zz)ρ + P_zz·(ZZ)·ρ·(ZZ)†
        if P_zz > 0
            for i in 1:(L-1)
                j = i + 1
                bra_i = 2i - 1
                ket_i = 2i
                bra_j = 2j - 1
                ket_j = 2j
                
                # Copy for unperturbed part
                ρ_copy = copy(ρ)
                
                # Apply ZZ to both bra and ket: (ZZ)·ρ·(ZZ)
                ZZ_bra = op(ZZ_2site, sites[bra_i], sites[bra_j])
                ZZ_ket = op(ZZ_2site, sites[ket_i], sites[ket_j])
                ρ_ZZ = apply([ZZ_bra, ZZ_ket], ρ_copy; cutoff=cutoff, maxdim=maxdim)
                
                # Linear combination: (1-P_zz)ρ + P_zz·(ZZ)·ρ·(ZZ)
                ρ = (1 - P_zz) * ρ + P_zz * ρ_ZZ
                # No renormalization - CPTP preserves trace
            end
        end
    end

    return MixedStateMPS(ρ), strobe_record
end


"""
    renyi2_binder_density_matrix_separate(L::Int; kwargs...)

Separate Rényi-2 Binder calculation (identical physics to original).

Uses:
  - evolve_renyi2_density_matrix_one_trial_separate() for evolution with strobe
  - MPO method for moment computation on doubled MPS
  - Ensemble averaging: B = 1 - <M4> / (3 <M2>²)

Returns ensemble-averaged results.
"""
function renyi2_binder_density_matrix_separate(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    P_x::Float64=0.0,
    P_zz::Float64=0.0,
    ntrials::Int=200,
    maxdim::Int=256,
    cutoff::Float64=1e-12,
    seed::Union{Nothing,Int}=nothing,
    strobe::Symbol=:after_full_layer,
    verbose::Bool=false,
)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)

    M2s = fill(NaN, ntrials)
    M4s = fill(NaN, ntrials)
    Bs = fill(NaN, ntrials)
    purities = fill(NaN, ntrials)
    max_linkdims = fill(0, ntrials)

    for t in 1:ntrials
        # Evolve one trial with separate Rényi-2 dynamics
        state, _ = evolve_renyi2_density_matrix_one_trial_separate(
            L;
            lambda_x=lambda_x,
            lambda_zz=lambda_zz,
            P_x=P_x,
            P_zz=P_zz,
            maxdim=maxdim,
            cutoff=cutoff,
            rng=rng,
            strobe=strobe,
        )

        # Extract doubled MPS and normalize
        ρ = get_mps(state)
        tr_val = doubledtrace(ρ)
        ρ_norm = ρ / tr_val

        # Purity = Tr(ρ²)
        purity = real(inner(ρ_norm, ρ_norm))
        purities[t] = purity

        # Max link dimension
        max_linkdims[t] = max_linkdim_mps(ρ_norm)

        # Use renyi2_binder_correlator_doubled from renyi2_binder.jl
        res = renyi2_binder_correlator_doubled(ρ_norm, L; layer=:bra)

        M2s[t] = res.M2
        M4s[t] = res.M4
        Bs[t] = res.B

        if verbose && (t == 1 || t == ntrials || t % max(1, ntrials ÷ 10) == 0)
            println("trial $t/$ntrials: B=$(res.B), M2=$(res.M2), M4=$(res.M4), purity=$purity")
        end
    end

    valid = isfinite.(M2s) .& isfinite.(M4s) .& isfinite.(Bs)
    n_valid = count(valid)
    n_invalid = ntrials - n_valid

    if n_valid == 0
        return (
            B = NaN,
            B_mean_of_trials = NaN,
            B_std_of_trials = NaN,
            M2_bar = NaN,
            M4_bar = NaN,
            purity_bar = NaN,
            M2s = M2s,
            M4s = M4s,
            Bs = Bs,
            purities = purities,
            ntrials = ntrials,
            n_valid = 0,
            n_invalid = n_invalid,
            max_linkdim = maximum(max_linkdims),
        )
    end

    M2_bar = mean(M2s[valid])
    M4_bar = mean(M4s[valid])
    purity_bar = mean(purities[valid])

    # Ensemble Binder: B = 1 - <M4> / (3 <M2>²)
    B_ensemble = 1.0 - M4_bar / (3.0 * M2_bar^2)

    return (
        B = B_ensemble,
        B_mean_of_trials = mean(Bs[valid]),
        B_std_of_trials = std(Bs[valid]),
        M2_bar = M2_bar,
        M4_bar = M4_bar,
        purity_bar = purity_bar,
        M2s = M2s,
        M4s = M4s,
        Bs = Bs,
        purities = purities,
        ntrials = ntrials,
        n_valid = n_valid,
        n_invalid = n_invalid,
        max_linkdim = maximum(max_linkdims),
    )
end


"""
    max_linkdim_mps(ψ::MPS)

Safe helper for bond-dimension diagnostics.
"""
function max_linkdim_mps(ψ::MPS)
    N = length(ψ)
    N <= 1 && return 1

    dims = Int[]
    for j in 1:(N - 1)
        l = commonind(ψ[j], ψ[j + 1])
        if l !== nothing
            push!(dims, dim(l))
        end
    end
    return isempty(dims) ? 1 : maximum(dims)
end
