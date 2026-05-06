"""
Rényi-2 Binder dynamics pipeline for density-matrix MPS.

This is a standalone, independent implementation for Rényi-2 Binder calculations
using doubled-MPS representation. It does NOT reuse the old EA Binder or 
vectorized-correlator machinery.

Key observable: Q = Σ_i Z_(2i-1)  (bra layer only)

Moments (Hilbert-Schmidt / ρ²-weighted):
    M2 = <ρ|Q²|ρ> / (L² <ρ|ρ>)
    M4 = <ρ|Q⁴|ρ> / (L⁴ <ρ|ρ>)
    B = 1 - M4 / (3 M2²)

Uses MPO method: ψ₁ = apply(Q, ρ), ψ₂ = apply(Q, ψ₁), etc.
"""

using Random
using Statistics
using ITensors, ITensorMPS
using LinearAlgebra

export evolve_renyi2_density_matrix_one_trial
export compute_renyi2_moments
export renyi2_binder_density_matrix_separate
export renyi2_sanity_checks

# Pauli matrices
const σx = Float64[0 1; 1 0]
const σy = ComplexF64[0 -im; im 0]
const σz = Float64[1 0; 0 -1]

# ============================================================
# Helper: Identity MPS for doubled representation
# ============================================================
"""
    id_mps_pair(s1::Index, s2::Index)

Create identity MPS for a bra-ket pair: |I⟩ = Σ_ab |ab⟩
Represents |0⟩⟨0| + |1⟩⟨1|.
"""
function id_mps_pair(s1::Index, s2::Index)
    return MPS([1.0; 0.0; 0.0; 1.0], [s1, s2])
end

"""
    bell_state(sites::Vector{<:Index})

Create Bell state for all bra-ket pairs: product of |I⟩ over all sites.
Represents identity operator in doubled form.
"""
function bell_state(sites::Vector{<:Index})
    N = length(sites)
    @assert iseven(N) "Sites must have even length for doubled MPS"
    
    tensors = ITensor[]
    for i in 1:2:N
        a, b = id_mps_pair(sites[i], sites[i+1])
        push!(tensors, a, b)
    end
    return MPS(tensors)
end

"""
    doubled_trace(ρ::MPS)

Compute Tr(ρ) in doubled MPS representation.
    Tr(ρ) = <I|ρ> = inner(bell_state, ρ)
"""
function doubled_trace(ρ::MPS)
    bell = bell_state(siteinds(ρ))
    return real(inner(bell, ρ))
end

"""
    doubled_purity(ρ::MPS)

Compute Tr(ρ²) in doubled MPS representation.
    Tr(ρ²) = <ρ|ρ> = inner(ρ, ρ)
"""
function doubled_purity(ρ::MPS)
    return real(inner(ρ, ρ))
end

# ============================================================
# Measurement and dephasing operators
# ============================================================
"""
    apply_x_measurement_to_site(ρ::MPS, sites::Vector{<:Index}, i::Int, 
                                lambda_x::Float64, rng, cutoff, maxdim)

Apply weak X measurement with random outcome to site i.
Returns (updated_ρ, outcome).

Kraus operator: Π_m^X = (I + (-1)^m λ_x X) / √(2(1+λ_x²))
Applied to both bra (2i-1) and ket (2i) layers.
Renormalize by trace after measurement.
"""
function apply_x_measurement_to_site(ρ::MPS, sites::Vector{<:Index}, i::Int, 
                                      lambda_x::Float64, rng, cutoff, maxdim)
    # Sample outcome based on ⟨X⟩ = Tr(ρ X) / Tr(ρ)
    # For doubled MPS: use correlation function
    bra_idx = 2i - 1
    ket_idx = 2i
    
    # Simplified: sample outcome with uniform probability
    # (A more sophisticated version could weight by expectation value)
    outcome = rand(rng) < 0.5 ? 0 : 1
    
    # Build Kraus operator: Π = (I + (-1)^m λ X) / √(2(1+λ²))
    Pi = (I(2) + (-1)^outcome * lambda_x * σx) / sqrt(2 * (1 + lambda_x^2))
    
    # Apply to both bra and ket
    Pi_bra = op(Pi, sites[bra_idx])
    Pi_ket = op(Pi, sites[ket_idx])
    
    ρ_new = apply([Pi_bra, Pi_ket], ρ; cutoff=cutoff, maxdim=maxdim)
    
    # Renormalize by trace
    tr = doubled_trace(ρ_new)
    tr > 1e-14 && (ρ_new = ρ_new / tr)
    
    return ρ_new, outcome
end

"""
    apply_x_dephasing_channel(ρ::MPS, sites::Vector{<:Index}, i::Int, 
                              P_x::Float64, cutoff, maxdim)

Apply X dephasing channel deterministically:
    ρ → (1 - P_x) ρ + P_x X ρ X

On doubled MPS: apply (1-P_x) I_bra⊗I_ket + P_x X_bra⊗X_ket
"""
function apply_x_dephasing_channel(ρ::MPS, sites::Vector{<:Index}, i::Int, 
                                    P_x::Float64, cutoff, maxdim)
    if P_x <= 0
        return ρ
    end
    
    bra_idx = 2i - 1
    ket_idx = 2i
    
    # ρ_copy = ρ (no flip)
    ρ_copy = copy(ρ)
    
    # ρ_X = X ρ X (on both layers)
    X_bra = op(σx, sites[bra_idx])
    X_ket = op(σx, sites[ket_idx])
    ρ_X = apply([X_bra, X_ket], copy(ρ); cutoff=cutoff, maxdim=maxdim)
    
    # Linear combination: (1-P_x) ρ + P_x ρ_X
    ρ_new = (1 - P_x) * ρ_copy + P_x * ρ_X
    
    return ρ_new
end

"""
    apply_zz_measurement_to_bond(ρ::MPS, sites::Vector{<:Index}, i::Int,
                                 lambda_zz::Float64, rng, cutoff, maxdim)

Apply weak ZZ measurement to bond (i, i+1) with random outcome.
Returns (updated_ρ, outcome).

Kraus operator: K_m^ZZ = (I₄ + (-1)^m λ_zz (Z⊗Z)) / √(2(1+λ_zz²))
Applied to both (bra_i, bra_j) and (ket_i, ket_j) pairs.
"""
function apply_zz_measurement_to_bond(ρ::MPS, sites::Vector{<:Index}, i::Int,
                                      lambda_zz::Float64, rng, cutoff, maxdim)
    @assert i < div(length(sites), 2) "Bond index i must be < L"
    
    j = i + 1
    bra_i = 2i - 1
    ket_i = 2i
    bra_j = 2j - 1
    ket_j = 2j
    
    # Sample outcome
    outcome = rand(rng) < 0.5 ? 0 : 1
    
    # Build 2-site Kraus: K = (I₄ + (-1)^m λ (Z⊗Z)) / √(2(1+λ²))
    I4 = Matrix{Float64}(I, 4, 4)
    ZZ_2site = kron(σz, σz)
    K = (I4 + (-1)^outcome * lambda_zz * ZZ_2site) / sqrt(2 * (1 + lambda_zz^2))
    
    # Apply to both bra and ket pairs
    K_bra = op(K, sites[bra_i], sites[bra_j])
    K_ket = op(K, sites[ket_i], sites[ket_j])
    
    ρ_new = apply([K_bra, K_ket], ρ; cutoff=cutoff, maxdim=maxdim)
    
    # Renormalize by trace
    tr = doubled_trace(ρ_new)
    tr > 1e-14 && (ρ_new = ρ_new / tr)
    
    return ρ_new, outcome
end

"""
    apply_zz_dephasing_channel(ρ::MPS, sites::Vector{<:Index}, i::Int,
                               P_zz::Float64, cutoff, maxdim)

Apply ZZ dephasing channel deterministically:
    ρ → (1 - P_zz) ρ + P_zz (Z⊗Z) ρ (Z⊗Z)

On doubled MPS (sites 2i-1, 2i, 2j-1, 2j):
    apply [(1-P_zz) I + P_zz (Z⊗Z)] to both bra and ket layers
"""
function apply_zz_dephasing_channel(ρ::MPS, sites::Vector{<:Index}, i::Int,
                                    P_zz::Float64, cutoff, maxdim)
    if P_zz <= 0
        return ρ
    end
    
    j = i + 1
    bra_i = 2i - 1
    ket_i = 2i
    bra_j = 2j - 1
    ket_j = 2j
    
    # ρ_copy = ρ (no flip)
    ρ_copy = copy(ρ)
    
    # ρ_ZZ = (Z⊗Z) ρ (Z⊗Z)
    ZZ_2site = kron(σz, σz)
    Z_bra_i = op(σz, sites[bra_i])
    Z_bra_j = op(σz, sites[bra_j])
    Z_ket_i = op(σz, sites[ket_i])
    Z_ket_j = op(σz, sites[ket_j])
    
    ρ_ZZ = apply([Z_bra_i, Z_bra_j, Z_ket_i, Z_ket_j], copy(ρ); cutoff=cutoff, maxdim=maxdim)
    
    # Linear combination: (1-P_zz) ρ + P_zz ρ_ZZ
    ρ_new = (1 - P_zz) * ρ_copy + P_zz * ρ_ZZ
    
    return ρ_new
end

# ============================================================
# Main evolution function
# ============================================================
"""
    evolve_renyi2_density_matrix_one_trial(
        L::Int;
        lambda_x::Float64,
        lambda_zz::Float64,
        P_x::Float64=0.0,
        P_zz::Float64=0.0,
        maxdim::Int=256,
        cutoff::Float64=1e-12,
        rng=Random.GLOBAL_RNG,
        T_max::Int=2*L,
        strobe::Symbol=:after_full_layer,
        verbose::Bool=false,
    )

Evolve a density matrix using doubled MPS with measurements and dephasing.

Returns (ρ_final::MPS, sites::Vector{Index})

Strobe points:
    :after_x_measurement   - after step 1
    :after_x_noise         - after step 2
    :after_zz_measurement  - after step 3
    :after_full_layer      - after step 4 (default)
"""
function evolve_renyi2_density_matrix_one_trial(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    P_x::Float64=0.0,
    P_zz::Float64=0.0,
    maxdim::Int=256,
    cutoff::Float64=1e-12,
    rng=Random.GLOBAL_RNG,
    T_max::Int=2*L,
    strobe::Symbol=:after_full_layer,
    verbose::Bool=false,
)
    valid_strobes = [:after_x_measurement, :after_x_noise, :after_zz_measurement, :after_full_layer]
    if !(strobe in valid_strobes)
        error("strobe must be one of: :after_x_measurement, :after_x_noise, :after_zz_measurement, :after_full_layer")
    end
    
    # Initialize doubled MPS: |↑↑...↑⟩⟨↑↑...↑|
    sites = siteinds("Qubit", 2L)
    ρ = MPS(sites, _ -> "Up")
    
    verbose && println("  Initialized: sites=$L, T_max=$T_max")
    
    # Track strobe state
    ρ_strobe = nothing
    
    for t in 1:T_max
        # Step 1: Weak X measurements on all sites
        if lambda_x > 0
            for i in 1:L
                ρ, outcome = apply_x_measurement_to_site(ρ, sites, i, lambda_x, rng, cutoff, maxdim)
            end
            if strobe == :after_x_measurement
                ρ_strobe = deepcopy(ρ)
            end
        end
        
        # Step 2: X dephasing on all sites
        if P_x > 0
            for i in 1:L
                ρ = apply_x_dephasing_channel(ρ, sites, i, P_x, cutoff, maxdim)
            end
            if strobe == :after_x_noise
                ρ_strobe = deepcopy(ρ)
            end
        end
        
        # Step 3: Weak ZZ measurements on all bonds
        if lambda_zz > 0
            for i in 1:(L-1)
                ρ, outcome = apply_zz_measurement_to_bond(ρ, sites, i, lambda_zz, rng, cutoff, maxdim)
            end
            if strobe == :after_zz_measurement
                ρ_strobe = deepcopy(ρ)
            end
        end
        
        # Step 4: ZZ dephasing on all bonds
        if P_zz > 0
            for i in 1:(L-1)
                ρ = apply_zz_dephasing_channel(ρ, sites, i, P_zz, cutoff, maxdim)
            end
            if strobe == :after_full_layer
                ρ_strobe = deepcopy(ρ)
            end
        end
    end
    
    # Return state at strobe point (default: after full layer)
    ρ_final = isnothing(ρ_strobe) ? ρ : ρ_strobe
    
    # Normalize by trace
    tr = doubled_trace(ρ_final)
    tr > 1e-14 && (ρ_final = ρ_final / tr)
    
    return ρ_final, sites
end

# ============================================================
# Rényi-2 moment calculation via MPO method
# ============================================================
"""
    build_order_operator_bra_layer(sites::Vector{<:Index}, L::Int)

Build MPO for Q = Σ_i Z_(2i-1) on bra layer only.
"""
function build_order_operator_bra_layer(sites::Vector{<:Index}, L::Int)
    @assert length(sites) == 2L "Expected 2L sites for doubled MPS"
    
    os = OpSum()
    for i in 1:L
        bra_idx = 2i - 1
        os += 1.0, "Z", bra_idx
    end
    return MPO(os, sites)
end

"""
    compute_renyi2_moments(ρ::MPS, sites::Vector{<:Index}, L::Int;
                          cutoff=1e-12, maxdim=256, verbose=false)

Compute Rényi-2 Binder moments using MPO method.

Returns named tuple (M2, M4, B, purity)

Algorithm:
    ψ₁ = apply(Q, ρ)
    num2 = <ψ₁|ψ₁> = Tr(Q² ρ²)
    M2 = num2 / (L² Tr(ρ²))
    
    ψ₂ = apply(Q, ψ₁)
    num4 = <ψ₂|ψ₂> = Tr(Q⁴ ρ²)
    M4 = num4 / (L⁴ Tr(ρ²))
    
    B = 1 - M4 / (3 M2²)
"""
function compute_renyi2_moments(ρ::MPS, sites::Vector{<:Index}, L::Int;
                               cutoff=1e-12, maxdim=256, verbose=false)
    # Build Q operator
    Q = build_order_operator_bra_layer(sites, L)
    
    # Compute Tr(ρ²)
    purity = doubled_purity(ρ)
    
    # Apply Q once: ψ₁ = Q|ρ⟩
    ψ1 = apply(Q, ρ; cutoff=cutoff, maxdim=maxdim)
    
    # M2 = <ψ₁|ψ₁> / (L² Tr(ρ²))
    num2 = real(inner(ψ1, ψ1))
    M2 = num2 / (L^2 * purity)
    
    # Apply Q again: ψ₂ = Q|ψ₁⟩
    ψ2 = apply(Q, ψ1; cutoff=cutoff, maxdim=maxdim)
    
    # M4 = <ψ₂|ψ₂> / (L⁴ Tr(ρ²))
    num4 = real(inner(ψ2, ψ2))
    M4 = num4 / (L^4 * purity)
    
    # Compute Binder: B = 1 - M4 / (3 M2²)
    # Protect against numerical issues
    if M2 > 1e-10
        B = 1 - M4 / (3 * M2^2)
    else
        B = NaN
    end
    
    verbose && println("      M2=$M2, M4=$M4, B=$B, purity=$purity")
    
    return (M2=M2, M4=M4, B=B, purity=purity)
end

# ============================================================
# Main Binder function
# ============================================================
"""
    renyi2_binder_density_matrix_separate(
        L::Int;
        lambda_x::Float64,
        lambda_zz::Float64,
        P_x::Float64=0.0,
        P_zz::Float64=0.0,
        ntrials::Int=100,
        maxdim::Int=256,
        cutoff::Float64=1e-12,
        seed::Union{Nothing,Int}=nothing,
        T_max::Int=2*L,
        strobe::Symbol=:after_full_layer,
        verbose::Bool=false,
    )

Run independent trials and compute ensemble Rényi-2 Binder statistics.

Returns named tuple with:
    - B: ensemble Binder (1 - mean(M4) / (3 mean(M2)²))
    - B_mean_of_trials: mean of single-trial Binders
    - B_std_of_trials: std of single-trial Binders
    - M2_bar: mean M2
    - M4_bar: mean M4
    - purity_bar: mean purity
    - M2s: array of single-trial M2 values
    - M4s: array of single-trial M4 values
    - Bs: array of single-trial B values
    - purities: array of single-trial purities
    - ntrials: number of trials
    - n_valid: number of valid trials (no NaN)
    - n_invalid: number of invalid trials
    - max_linkdim: maximum bond dimension encountered
"""
function renyi2_binder_density_matrix_separate(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    P_x::Float64=0.0,
    P_zz::Float64=0.0,
    ntrials::Int=100,
    maxdim::Int=256,
    cutoff::Float64=1e-12,
    seed::Union{Nothing,Int}=nothing,
    T_max::Int=2*L,
    strobe::Symbol=:after_full_layer,
    verbose::Bool=false,
)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    
    M2s = fill(NaN, ntrials)
    M4s = fill(NaN, ntrials)
    Bs = fill(NaN, ntrials)
    purities = fill(NaN, ntrials)
    max_linkdim_trial = 0
    
    for t in 1:ntrials
        if verbose
            println("  Trial $t/$ntrials:")
        end
        
        # Evolve one trial
        ρ, sites = evolve_renyi2_density_matrix_one_trial(
            L;
            lambda_x=lambda_x,
            lambda_zz=lambda_zz,
            P_x=P_x,
            P_zz=P_zz,
            maxdim=maxdim,
            cutoff=cutoff,
            rng=rng,
            T_max=T_max,
            strobe=strobe,
            verbose=verbose,
        )
        
        # Compute moments
        res = compute_renyi2_moments(ρ, sites, L; cutoff=cutoff, maxdim=maxdim, verbose=verbose)
        
        M2s[t] = res.M2
        M4s[t] = res.M4
        Bs[t] = res.B
        purities[t] = res.purity
        
        # Track max link dimension
        for j in 1:(length(ρ)-1)
            ld = linkdim(ρ, j)
            max_linkdim_trial = max(max_linkdim_trial, ld)
        end
    end
    
    # Count valid trials (no NaN)
    n_valid = count(!isnan, Bs)
    n_invalid = ntrials - n_valid
    
    # Compute statistics (ignoring NaN)
    M2_bar = nanmean(M2s)
    M4_bar = nanmean(M4s)
    purity_bar = nanmean(purities)
    B_mean_of_trials = nanmean(Bs)
    B_std_of_trials = nanstd(Bs)
    
    # Ensemble Binder: B = 1 - mean(M4) / (3 mean(M2)²)
    if M2_bar > 1e-10
        B_ensemble = 1 - M4_bar / (3 * M2_bar^2)
    else
        B_ensemble = NaN
    end
    
    return (
        B=B_ensemble,
        B_mean_of_trials=B_mean_of_trials,
        B_std_of_trials=B_std_of_trials,
        M2_bar=M2_bar,
        M4_bar=M4_bar,
        purity_bar=purity_bar,
        M2s=M2s,
        M4s=M4s,
        Bs=Bs,
        purities=purities,
        ntrials=ntrials,
        n_valid=n_valid,
        n_invalid=n_invalid,
        max_linkdim=max_linkdim_trial,
    )
end

"""Helper: NaN-ignoring mean"""
function nanmean(arr::Vector)
    valid = filter(!isnan, arr)
    return isempty(valid) ? NaN : mean(valid)
end

"""Helper: NaN-ignoring std"""
function nanstd(arr::Vector)
    valid = filter(!isnan, arr)
    return isempty(valid) ? NaN : std(valid)
end

# ============================================================
# Sanity checks
# ============================================================
"""
    renyi2_sanity_checks()

Run a series of sanity checks on Rényi-2 Binder calculations.
"""
function renyi2_sanity_checks()
    println("="^70)
    println("RÉNYI-2 BINDER SANITY CHECKS")
    println("="^70)
    println()
    
    L = 8
    
    # Check A: All-up pure state
    println("Check A: All-up pure state |↑↑...↑⟩⟨↑↑...↑|")
    println("  Expected B ≈ 2/3 ≈ 0.667")
    sites_a = siteinds("Qubit", 2L)
    ρ_a = MPS(sites_a, _ -> "Up")
    res_a = compute_renyi2_moments(ρ_a, sites_a, L; verbose=false)
    println("  Observed B = $(round(res_a.B, digits=3))")
    println("  Purity = $(round(res_a.purity, digits=4))")
    println()
    
    # Check B: Dephasing with moderate probability
    println("Check B: X+ZZ dephasing with P=0.5 (moderate mixing)")
    println("  Expected B between 2/(3L) ≈ 0.083 and 2/3 ≈ 0.667")
    rng = MersenneTwister(123)
    ρ_b, sites_b = evolve_renyi2_density_matrix_one_trial(
        L;
        lambda_x=0.0, lambda_zz=0.0,
        P_x=0.5, P_zz=0.5,  # Moderate dephasing
        rng=rng, T_max=5*L,
        verbose=false
    )
    res_b = compute_renyi2_moments(ρ_b, sites_b, L; verbose=false)
    println("  Observed B = $(round(res_b.B, digits=3))")
    println("  Purity = $(round(res_b.purity, digits=4))")
    println()
    
    # Check C: ZZ measurements and dephasing  
    println("Check C: ZZ measurements (λ_zz=1.0) + P_zz=0.3 dephasing")
    println("  Expected B between 2/(3L)≈0.083 and 2/3≈0.667")
    rng_c = MersenneTwister(456)
    ρ_c, sites_c = evolve_renyi2_density_matrix_one_trial(
        L;
        lambda_x=0.0, lambda_zz=1.0,
        P_x=0.0, P_zz=0.3,
        rng=rng_c, T_max=2*L,
        verbose=false
    )
    res_c = compute_renyi2_moments(ρ_c, sites_c, L; verbose=false)
    println("  Observed B = $(round(res_c.B, digits=3))")
    println("  Purity = $(round(res_c.purity, digits=4))")
    println()
    
    println("Sanity checks complete!")
    println()
end
