module BinderSimManual

using Random, Statistics
using ITensors, ITensorMPS
include("ManualCorrelators.jl")
using .ManualCorrelators

export ea_binder_mc_manual, create_weak_measurement_operators, create_up_state_mps, evolve_one_trial

function create_weak_measurement_operators(sites, lambda_x::Float64, lambda_zz::Float64)
    WEAK_X_0 = Dict{Int,ITensor}()
    WEAK_X_1 = Dict{Int,ITensor}()
    WEAK_ZZ_0 = Dict{Tuple{Int,Int},ITensor}()
    WEAK_ZZ_1 = Dict{Tuple{Int,Int},ITensor}()

    norm_x  = sqrt(2 * (1 + lambda_x^2))
    norm_zz = sqrt(2 * (1 + lambda_zz^2))

    for i in 1:length(sites)
        Id_i = op("Id", sites[i])
        X_i  = 2 * op("Sx", sites[i])
        WEAK_X_0[i] = (Id_i + lambda_x * X_i) / norm_x
        WEAK_X_1[i] = (Id_i - lambda_x * X_i) / norm_x
    end
    for i in 1:(length(sites)-1)
        Z_i = 2 * op("Sz", sites[i])
        Z_j = 2 * op("Sz", sites[i+1])
        II  = op("Id", sites[i]) * op("Id", sites[i+1])
        ZZ  = Z_i * Z_j
        WEAK_ZZ_0[(i,i+1)] = (II + lambda_zz * ZZ) / norm_zz
        WEAK_ZZ_1[(i,i+1)] = (II - lambda_zz * ZZ) / norm_zz
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
        # weak X on all sites
        for i in 1:L
            ψ = sample_and_apply(ψ, KX0[i], KX1[i], [i]; maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
        # weak ZZ simultaneous: all adjacent bonds (physically correct)
        for i in 1:(L-1)
            ψ = sample_and_apply(ψ, KZZ0[(i,i+1)], KZZ1[(i,i+1)], [i,i+1];
                                 maxdim=maxdim, cutoff=cutoff, rng=rng)
        end
    end
    return ψ, sites
end

"""
    ea_binder_mc_manual(L::Int; kwargs...)

Edwards-Anderson Binder parameter calculation using manual tensor contractions.
This bypasses ITensorCorrelators and should work for larger systems.
"""
function ea_binder_mc_manual(L::Int; lambda_x::Float64, lambda_zz::Float64,
                            ntrials::Int=200, maxdim::Int=256, cutoff::Float64=1e-12,
                            seed::Union{Nothing,Int}=nothing,
                            manual_maxdim::Int=64, manual_cutoff::Float64=1e-10,
                            chunk_size::Int=100)
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    
    println("Running Edwards-Anderson Binder parameter calculation with manual correlators")
    println("L = $L, ntrials = $ntrials, manual_maxdim = $manual_maxdim")
    
    Bs = Vector{Float64}(undef, ntrials)
    S2s = Vector{Float64}(undef, ntrials)
    S4s = Vector{Float64}(undef, ntrials)

    for t in 1:ntrials
        println("Trial $t/$ntrials")
        
        # Evolve one trial trajectory
        ψ, sites = evolve_one_trial(L; lambda_x=lambda_x, lambda_zz=lambda_zz,
                                    maxdim=maxdim, cutoff=cutoff, rng=rng)
        
        # Compute Binder parameter using manual correlators
        try
            if L <= 14
                # Use full calculation for smaller systems
                result = ManualCorrelators.manual_binder_parameter(ψ, sites; 
                                               maxdim=manual_maxdim, cutoff=manual_cutoff,
                                               chunk_size=chunk_size)
            else
                # Use chunked calculation for larger systems
                result = ManualCorrelators.chunked_manual_binder_parameter(ψ, sites; 
                                                       maxdim=manual_maxdim, cutoff=manual_cutoff,
                                                       max_4point_chunk=chunk_size)
            end
            
            Bs[t] = result.B
            S2s[t] = result.S2
            S4s[t] = result.S4
            
            println("  Trial $t: B = $(result.B)")
            
        catch e
            println("Error in trial $t: $e")
            Bs[t] = NaN
            S2s[t] = NaN
            S4s[t] = NaN
        end
    end

    # Compute ensemble averages
    valid_trials = .!isnan.(Bs)
    if sum(valid_trials) == 0
        println("Warning: No valid trials completed!")
        return (B = NaN, B_mean_of_trials = NaN, B_std_of_trials = NaN,
                S2_bar = NaN, S4_bar = NaN, ntrials = ntrials,
                ntrials_completed = 0)
    end
    
    B_mean = mean(Bs[valid_trials])
    B_std = std(Bs[valid_trials])
    S2_bar = mean(S2s[valid_trials])
    S4_bar = mean(S4s[valid_trials])
    B_EA = 1.0 - S4_bar / (3.0 * max(S2_bar^2, 1e-30))
    
    ntrials_completed = sum(valid_trials)
    
    println("Completed $ntrials_completed/$ntrials trials")
    println("Ensemble-averaged Binder parameter: $B_EA")
    
    return (B = B_EA,
            B_mean_of_trials = B_mean,
            B_std_of_trials = B_std,
            S2_bar = S2_bar,
            S4_bar = S4_bar,
            ntrials = ntrials,
            ntrials_completed = ntrials_completed)
end

end
