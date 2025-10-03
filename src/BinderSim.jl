module BinderSim

using Random, Statistics
using ITensors, ITensorMPS, ITensorCorrelators

export ea_binder_mc, create_weak_measurement_operators, create_up_state_mps, evolve_one_trial

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

function sample_and_apply(ψ::MPS, K0::ITensor, K1::ITensor, which::Vector{Int};
                          maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    ϕ0 = product(K0, ψ, which; maxdim=maxdim, cutoff=cutoff)
    p0 = max(real(inner(dag(ϕ0), ϕ0)), 0.0)
    ϕ1 = product(K1, ψ, which; maxdim=maxdim, cutoff=cutoff)
    p1 = max(real(inner(dag(ϕ1), ϕ1)), 0.0)
    # Check
    tot = p0 + p1
    ϕ = (tot <= 0) ? ϕ0 : (rand(rng) < p0/tot ? ϕ0 : ϕ1)
    orthogonalize!(ϕ, which[end])
    normalize!(ϕ)
    return ϕ
end

# One trajectory, open chain, T_max = 2L
function evolve_one_trial(L::Int; lambda_x::Float64, lambda_zz::Float64,
                          maxdim::Int=256, cutoff::Float64=1e-12, rng=Random.GLOBAL_RNG)
    ψ, sites = create_up_state_mps(L)
    KX0, KX1, KZZ0, KZZ1 = create_weak_measurement_operators(sites, lambda_x, lambda_zz)
    T_max = 2L
    for _ in 1:T_max
        # weak X on all sites - random sampling
        # check random
        for i in 1:L
            ψ = sample_and_apply(ψ, KX0[i], KX1[i], [i]; maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
        # weak ZZ simultaneous: all adjacent bonds - random sampling (physically correct)
        for i in 1:(L-1)
            ψ = sample_and_apply(ψ, KZZ0[(i,i+1)], KZZ1[(i,i+1)], [i,i+1];
                                 maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
    end
    return ψ, sites
end

function ea_binder_mc(L::Int; lambda_x::Float64, lambda_zz::Float64,
                      ntrials::Int=200, maxdim::Int=256, cutoff::Float64=1e-12,
                      chunk4::Int=50_000, seed::Union{Nothing,Int}=nothing)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    S2s = Vector{Float64}(undef, ntrials)
    S4s = Vector{Float64}(undef, ntrials)
    Bs  = Vector{Float64}(undef, ntrials)

    for t in 1:ntrials
        ψ, sites = evolve_one_trial(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                    maxdim=maxdim, cutoff=cutoff, rng=rng)
        
        # Use the working approach from the notebook
        # --- put the center in the middle and ensure norm 1 ---
        orthogonalize!(ψ, cld(length(sites),2))
        normalize!(ψ)

        # Use ALL sites for correlation calculation (not just central fraction)
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
