module BipartiteEntropy

using Random, Statistics
using ITensors, ITensorMPS
using JSON

export weak_bipartite_entropy

function create_weak_measurement_operators(sites, lambda_x::Float64, lambda_zz::Float64)
    WEAK_X_0 = Dict{Int,ITensor}()
    WEAK_X_1 = Dict{Int,ITensor}()
    WEAK_ZZ_0 = Dict{Tuple{Int,Int},ITensor}()
    WEAK_ZZ_1 = Dict{Tuple{Int,Int},ITensor}()

    # Normalization factors
    norm_x = sqrt(2 * (1 + lambda_x^2))
    norm_zz = sqrt(2 * (1 + lambda_zz^2))

    # Create weak X operators for each site
    for i in 1:length(sites)
        Id_i = op("Id", sites[i])
        X_i = op("X", sites[i])
        
        # K_0 = (I + λ_x X) / √(2(1 + λ_x²))
        WEAK_X_0[i] = (Id_i + lambda_x * X_i) / norm_x
        # K_1 = (I - λ_x X) / √(2(1 + λ_x²))  
        WEAK_X_1[i] = (Id_i - lambda_x * X_i) / norm_x
    end

    # Create weak ZZ operators for adjacent pairs
    for i in 1:(length(sites)-1)
        Z_i = op("Z", sites[i])
        Z_j = op("Z", sites[i+1])
        Id_i = op("Id", sites[i])
        Id_j = op("Id", sites[i+1])
        
        ZZ = Z_i * Z_j
        II = Id_i * Id_j
        
        # K_0 = (II + λ_zz ZZ) / √(2(1 + λ_zz²))
        WEAK_ZZ_0[(i, i+1)] = (II + lambda_zz * ZZ) / norm_zz
        # K_1 = (II - λ_zz ZZ) / √(2(1 + λ_zz²))
        WEAK_ZZ_1[(i, i+1)] = (II - lambda_zz * ZZ) / norm_zz
    end

    return WEAK_X_0, WEAK_X_1, WEAK_ZZ_0, WEAK_ZZ_1
end

function create_plus_state_mps(L::Int)
    sites = siteinds("S=1/2", L)
    state = ["+" for _ in 1:L]
    ψ = productMPS(sites, state)
    return ψ, sites
end

function apply_single_site_gate(ψ::MPS, site::Int, gate::ITensor; 
                                maxdim::Int=128, cutoff::Float64=1e-12)
    ψ_new = copy(ψ)
    ψ_new[site] = gate * ψ_new[site]
    noprime!(ψ_new[site])

    if maxlinkdim(ψ_new) > maxdim
        ψ_new = truncate(ψ_new; maxdim=maxdim, cutoff=cutoff)
    end
    
    return ψ_new
end

function apply_two_site_gate(ψ::MPS, site1::Int, site2::Int, gate::ITensor; 
                            maxdim::Int=128, cutoff::Float64=1e-12)
    if abs(site1 - site2) == 1
        s1, s2 = min(site1, site2), max(site1, site2)
        return product(gate, ψ, [s1, s2]; maxdim=maxdim, cutoff=cutoff)
    else
        return ψ
    end
end

function apply_weak_x_measurement(ψ::MPS, site::Int, prob::Float64,
                                 WEAK_X_0::Dict, WEAK_X_1::Dict;
                                 maxdim::Int=128, cutoff::Float64=1e-12)
    if rand() >= prob
        return ψ
    end
    
    ψ_0 = apply_single_site_gate(ψ, site, WEAK_X_0[site]; maxdim=maxdim, cutoff=cutoff)
    ψ_1 = apply_single_site_gate(ψ, site, WEAK_X_1[site]; maxdim=maxdim, cutoff=cutoff)
    
    prob_0 = real(inner(ψ_0, ψ_0))
    prob_1 = real(inner(ψ_1, ψ_1))
    
    prob_0 = max(prob_0, 0.0)
    prob_1 = max(prob_1, 0.0)
    total_prob = prob_0 + prob_1
    
    if total_prob < 1e-12
        return nothing
    end
    
    prob_0 /= total_prob
    prob_1 /= total_prob
    
    if rand() < prob_0
        if prob_0 < 1e-12
            return nothing
        end
        normalize!(ψ_0)
        return ψ_0
    else
        if prob_1 < 1e-12
            return nothing
        end
        normalize!(ψ_1)
        return ψ_1
    end
end

function apply_weak_zz_measurement(ψ::MPS, site1::Int, site2::Int, prob::Float64,
                                  WEAK_ZZ_0::Dict, WEAK_ZZ_1::Dict;
                                  maxdim::Int=128, cutoff::Float64=1e-12)
    if rand() >= prob
        return ψ
    end
    
    ψ_0 = apply_two_site_gate(ψ, site1, site2, WEAK_ZZ_0[(site1, site2)]; 
                             maxdim=maxdim, cutoff=cutoff)
    ψ_1 = apply_two_site_gate(ψ, site1, site2, WEAK_ZZ_1[(site1, site2)]; 
                             maxdim=maxdim, cutoff=cutoff)
    
    prob_0 = real(inner(ψ_0, ψ_0))
    prob_1 = real(inner(ψ_1, ψ_1))
    
    prob_0 = max(prob_0, 0.0)
    prob_1 = max(prob_1, 0.0)
    total_prob = prob_0 + prob_1
    
    if total_prob < 1e-12
        return nothing
    end
    
    prob_0 /= total_prob
    prob_1 /= total_prob
    
    if rand() < prob_0
        if prob_0 < 1e-12
            return nothing
        end
        normalize!(ψ_0)
        return ψ_0
    else
        if prob_1 < 1e-12
            return nothing
        end
        normalize!(ψ_1)
        return ψ_1
    end
end

function calculate_bipartite_entropy(ψ::MPS, cut::Int)
    """
    Calculate bipartite entanglement entropy across cut position
    """
    L = length(ψ)
    if cut <= 0 || cut >= L
        return 0.0
    end
    
    try
        # Create copy and move orthogonality center to cut position
        ψ_copy = copy(ψ)
        orthogonalize!(ψ_copy, cut)
        
        # Get bond between sites cut and cut+1
        if cut < L
            bond_idx = commonind(ψ_copy[cut], ψ_copy[cut+1])
            if bond_idx !== nothing
                # Perform SVD: M = USV†
                U, S, V = svd(ψ_copy[cut], (bond_idx,))
                
                # Extract Schmidt values (singular values)
                schmidt_vals = Float64[]
                for i in 1:dim(S, 1)
                    val = S[i,i]
                    if abs(val) > 1e-12
                        push!(schmidt_vals, abs(val))
                    end
                end
                
                if isempty(schmidt_vals)
                    return 0.0
                end
                
                # Normalize Schmidt values
                schmidt_vals = schmidt_vals ./ norm(schmidt_vals)
                
                # Calculate von Neumann entropy: S = -Σ λᵢ log₂(λᵢ)
                entropy = -sum(s^2 * log2(s^2 + 1e-16) for s in schmidt_vals)
                return entropy
            end
        end
        
        return 0.0
    catch e
        println("Warning: Entropy calculation failed for cut=$cut: $e")
        return 0.0
    end
end

function run_single_trial_entropy(seed::Int, L::Int, T_max::Int, lambda_x::Float64, lambda_zz::Float64, 
                                 cuts::Vector{Int}; maxdim::Int=128, cutoff::Float64=1e-12)
    """
    Run single trial and calculate bipartite entropy at specified cuts
    """
    Random.seed!(seed)
    
    ψ, sites = create_plus_state_mps(L)
    
    WEAK_X_0, WEAK_X_1, WEAK_ZZ_0, WEAK_ZZ_1 = create_weak_measurement_operators(sites, lambda_x, lambda_zz)
    
    for t in 1:T_max
        # Apply weak X measurements on all sites
        for i in 1:L
            ψ = apply_weak_x_measurement(ψ, i, 1.0, WEAK_X_0, WEAK_X_1; 
                                        maxdim=maxdim, cutoff=cutoff)
            if ψ === nothing
                return fill(NaN, length(cuts))
            end
        end
        
        # Apply weak ZZ measurements on adjacent pairs
        for i in 1:(L-1)
            ψ = apply_weak_zz_measurement(ψ, i, i+1, 1.0, WEAK_ZZ_0, WEAK_ZZ_1;
                                         maxdim=maxdim, cutoff=cutoff)
            if ψ === nothing
                return fill(NaN, length(cuts))
            end
        end
    end
    
    # Calculate bipartite entropies at all requested cuts
    entropies = [calculate_bipartite_entropy(ψ, cut) for cut in cuts]
    return entropies
end

function weak_bipartite_entropy(L::Int; lambda_x::Float64, lambda_zz::Float64, T_max::Int,
                               ntrials::Int=100, maxdim::Int=128, cutoff::Float64=1e-12,
                               seed::Union{Nothing,Int}=nothing)
    """
    Calculate average bipartite entanglement entropy across all cuts for given parameters
    """
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    
    # Calculate entropy at all possible bipartite cuts
    cuts = collect(1:(L-1))
    
    # Storage for results
    all_entropies = Vector{Vector{Float64}}(undef, ntrials)
    
    println("Running $ntrials trials for L=$L, T_max=$T_max, λ_x=$lambda_x, λ_zz=$lambda_zz")
    
    for t in 1:ntrials
        trial_seed = rand(rng, 1:1000000)
        entropies = run_single_trial_entropy(trial_seed, L, T_max, lambda_x, lambda_zz, cuts;
                                            maxdim=maxdim, cutoff=cutoff)
        all_entropies[t] = entropies
        
        if t % 10 == 0
            println("  Completed $t/$ntrials trials")
        end
    end
    
    # Process results
    n_cuts = length(cuts)
    avg_entropies = Vector{Float64}(undef, n_cuts)
    std_entropies = Vector{Float64}(undef, n_cuts)
    
    for i in 1:n_cuts
        # Extract entropy values for cut i across all trials
        cut_entropies = [all_entropies[t][i] for t in 1:ntrials]
        valid_entropies = filter(!isnan, cut_entropies)
        
        if !isempty(valid_entropies)
            avg_entropies[i] = mean(valid_entropies)
            std_entropies[i] = std(valid_entropies)
        else
            avg_entropies[i] = NaN
            std_entropies[i] = NaN
        end
    end
    
    # Calculate central entropy and area law coefficient
    central_cut = div(L, 2)
    central_entropy = avg_entropies[central_cut]
    
    # Simple area law analysis: S ~ a * log(L) + const
    # For 1D systems, "area" is just the boundary (constant), so S ~ const + corrections
    max_entropy = maximum(filter(!isnan, avg_entropies))
    
    return (
        L = L,
        T_max = T_max,
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        ntrials = ntrials,
        cuts = cuts,
        avg_entropies = avg_entropies,
        std_entropies = std_entropies,
        central_entropy = central_entropy,
        max_entropy = max_entropy,
        successful_trials = count(t -> !all(isnan.(all_entropies[t])), 1:ntrials)
    )
end

end # module