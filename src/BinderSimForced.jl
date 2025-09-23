module BinderSimForced

using Random, Statistics
using ITensors, ITensorMPS, ITensorCorrelators

export ea_binder_mc_forced, create_weak_measurement_operators, create_up_state_mps, evolve_one_trial_forced

function create_weak_measurement_operators(sites, lambda_x::Float64, lambda_zz::Float64)
    WEAK_X_0 = Dict{Int,ITensor}()
    WEAK_X_1 = Dict{Int,ITensor}()
    WEAK_ZZ_0 = Dict{Tuple{Int,Int},ITensor}()
    WEAK_ZZ_1 = Dict{Tuple{Int,Int},ITensor}()

    norm_x  = sqrt(2 * (1 + lambda_x^2))
    norm_zz = sqrt(2 * (1 + lambda_zz^2))

    for i in 1:length(sites)
        Id_i = op("Id", sites[i])
        # Match sparse matrix exactly: X = [[0,1],[1,0]]
        # ITensor Sx = 0.5 * Pauli_X, so 2*Sx = Pauli_X
        X_i  = 2 * op("Sx", sites[i])
        # Match sparse: (I + (-1)^outcome * lam * X) / norm
        WEAK_X_0[i] = (Id_i + lambda_x * X_i) / norm_x     # outcome = 0
        WEAK_X_1[i] = (Id_i - lambda_x * X_i) / norm_x     # outcome = 1
    end
    for i in 1:(length(sites)-1)
        # Match sparse matrix exactly: Z = [[1,0],[0,-1]]  
        # ITensor Sz = 0.5 * Pauli_Z, so 2*Sz = Pauli_Z
        Z_i = 2 * op("Sz", sites[i])
        Z_j = 2 * op("Sz", sites[i+1])
        II  = op("Id", sites[i]) * op("Id", sites[i+1])
        ZZ  = Z_i * Z_j
        # Match sparse: (I⊗I + (-1)^outcome * lam * ZZ) / norm
        WEAK_ZZ_0[(i,i+1)] = (II + lambda_zz * ZZ) / norm_zz   # outcome = 0
        WEAK_ZZ_1[(i,i+1)] = (II - lambda_zz * ZZ) / norm_zz   # outcome = 1
    end
    return WEAK_X_0, WEAK_X_1, WEAK_ZZ_0, WEAK_ZZ_1
end

function create_up_state_mps(L::Int)
    sites = siteinds("S=1/2", L)
    ψ = productMPS(sites, fill("Up", L))
    return ψ, sites
end

"""
    apply_forced_measurement(ψ::MPS, K1::ITensor, which::Vector{Int}; kwargs...)

Apply a forced +1 measurement outcome. Always uses the K1 (outcome = 1) operator,
corresponding to the (-1)^1 = -1 case in the weak measurement formulation.
"""
function apply_forced_measurement(ψ::MPS, K1::ITensor, which::Vector{Int};
                                 maxdim::Int=256, cutoff::Float64=1e-12)
    # Always apply the +1 outcome operator (K1)
    ϕ = product(K1, ψ, which; maxdim=maxdim, cutoff=cutoff)
    orthogonalize!(ϕ, which[end])
    normalize!(ϕ)
    return ϕ
end

"""
    evolve_one_trial_forced(L::Int; lambda_x, lambda_zz, maxdim, cutoff)

Evolution with forced +1 measurement outcomes for all weak measurements.
This creates a deterministic (non-random) trajectory where all measurements
give the +1 result.
"""
function evolve_one_trial_forced(L::Int; lambda_x::Float64, lambda_zz::Float64,
                                maxdim::Int=256, cutoff::Float64=1e-12)
    ψ, sites = create_up_state_mps(L)
    KX0, KX1, KZZ0, KZZ1 = create_weak_measurement_operators(sites, lambda_x, lambda_zz)
    T_max = 2L
    
    println("  Evolving with forced +1 measurements for $T_max time steps...")
    
    for t in 1:T_max
        # weak X on all sites - FORCED +1 outcomes
        for i in 1:L
            ψ = apply_forced_measurement(ψ, KX1[i], [i]; maxdim=maxdim, cutoff=cutoff)
        end
        
        # weak ZZ simultaneous: all adjacent bonds - FORCED +1 outcomes (physically correct)
        for i in 1:(L-1)
            ψ = apply_forced_measurement(ψ, KZZ1[(i,i+1)], [i,i+1]; maxdim=maxdim, cutoff=cutoff)
        end
        
        # Progress indicator
        if t % (T_max ÷ 4) == 0
            bond_dim = maxlinkdim(ψ)
            println("    Time step $t/$T_max, max bond dimension: $bond_dim")
        end
    end
    return ψ, sites
end

"""
    ea_binder_mc_forced(L::Int; kwargs...)

Edwards-Anderson Binder parameter calculation with forced +1 measurement outcomes.
Since all trajectories are now deterministic (same forced outcomes), we only need
ntrials to check reproducibility - they should all give identical results.
"""
function ea_binder_mc_forced(L::Int; lambda_x::Float64, lambda_zz::Float64,
                            ntrials::Int=1, maxdim::Int=256, cutoff::Float64=1e-12,
                            chunk4::Int=50_000, seed::Union{Nothing,Int}=nothing)
    
    println("Edwards-Anderson Binder parameter with FORCED +1 measurements")
    println("L = $L, λₓ = $lambda_x, λ_zz = $lambda_zz")
    println("Note: All measurements forced to +1 outcome (deterministic)")
    
    S2s = Vector{Float64}(undef, ntrials)
    S4s = Vector{Float64}(undef, ntrials)
    Bs  = Vector{Float64}(undef, ntrials)

    for t in 1:ntrials
        println("\\nTrial $t/$ntrials:")
        
        ψ, sites = evolve_one_trial_forced(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                          maxdim=maxdim, cutoff=cutoff)
        
        # Use the working approach from the notebook
        # --- put the center in the middle and ensure norm 1 ---
        orthogonalize!(ψ, cld(length(sites),2))
        normalize!(ψ)

        # Use ALL sites for correlation calculation (not just central fraction)
        L_total = length(sites)
        idx = 1:L_total
        n = length(idx)

        println("  Computing correlators on ALL sites: 1 to $L_total ($n sites)")

        pairs = [(i,j) for i in idx for j in idx]
        quads = [(i,j,k,l) for i in idx for j in idx for k in idx for l in idx]

        println("  Computing $(length(pairs)) 2-point and $(length(quads)) 4-point correlators...")

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
        
        println("  Trial $t results: M2sq = $M2sq, M4sq = $M4sq, B = $(Bs[t])")
    end

    S2_bar = mean(S2s)
    S4_bar = mean(S4s)
    B_EA   = 1.0 - S4_bar / (3.0*S2_bar^2 + 1e-30)

    println("\\nFinal ensemble results:")
    println("  Ensemble B_EA = $B_EA")
    println("  Mean of trials = $(mean(Bs))")
    println("  Std of trials = $(std(Bs))")

    return (B = B_EA,
            B_mean_of_trials = mean(Bs),
            B_std_of_trials  = std(Bs),
            S2_bar = S2_bar,
            S4_bar = S4_bar,
            ntrials = ntrials)
end

end