module ManualCorrelators

using ITensors, ITensorMPS
using Random, Statistics

export manual_binder_parameter, manual_correlator_2point, manual_correlator_4point

"""
    manual_correlator_2point(ψ::MPS, sites, i, j; maxdim=64, cutoff=1e-10)

Compute ⟨Sz_i Sz_j⟩ using direct tensor contractions instead of ITensorCorrelators.
This bypasses the scaling limitations of ITensorCorrelators.
"""
function manual_correlator_2point(ψ::MPS, sites, i, j; maxdim=64, cutoff=1e-10)
    # Make a copy to avoid modifying original
    ψ_work = copy(ψ)
    
    # Orthogonalize around the leftmost operator site
    left_site = min(i, j)
    orthogonalize!(ψ_work, left_site)
    
    # Apply Sz operators
    Sz_i = op("Sz", sites[i])
    Sz_j = op("Sz", sites[j])
    
    # Apply first operator
    ψ_work[i] = Sz_i * ψ_work[i]
    noprime!(ψ_work[i])
    
    # Apply second operator
    ψ_work[j] = Sz_j * ψ_work[j]
    noprime!(ψ_work[j])
    
    # Truncate if needed
    if maxlinkdim(ψ_work) > maxdim
        ψ_work = truncate(ψ_work; maxdim=maxdim, cutoff=cutoff)
    end
    
    # Compute inner product
    return real(inner(ψ, ψ_work))
end

"""
    manual_correlator_4point(ψ::MPS, sites, i, j, k, l; maxdim=64, cutoff=1e-10)

Compute ⟨Sz_i Sz_j Sz_k Sz_l⟩ using direct tensor contractions.
Applies operators sequentially with truncation between each application.
"""
function manual_correlator_4point(ψ::MPS, sites, i, j, k, l; maxdim=64, cutoff=1e-10)
    # Make a copy to avoid modifying original
    ψ_work = copy(ψ)
    
    # Sort sites for optimal orthogonalization
    sites_ordered = sort([i, j, k, l])
    left_site = sites_ordered[1]
    
    # Orthogonalize around the leftmost operator site
    orthogonalize!(ψ_work, left_site)
    
    # Apply operators in spatial order to minimize entanglement growth
    for site_idx in sites_ordered
        Sz = op("Sz", sites[site_idx])
        ψ_work[site_idx] = Sz * ψ_work[site_idx]
        noprime!(ψ_work[site_idx])
        
        # Truncate after each operator to control bond dimension
        if maxlinkdim(ψ_work) > maxdim
            ψ_work = truncate(ψ_work; maxdim=maxdim, cutoff=cutoff)
            normalize!(ψ_work)
        end
    end
    
    # Compute inner product
    return real(inner(ψ, ψ_work))
end

"""
    manual_binder_parameter(ψ::MPS, sites; maxdim=64, cutoff=1e-10, chunk_size=1000)

Compute the Edwards-Anderson Binder parameter using manual tensor contractions.
This should work for larger systems where ITensorCorrelators fails.

B = 1 - ⟨M⁴⟩/(3⟨M²⟩²) where M² = Σᵢⱼ ⟨Sz_i Sz_j⟩²
"""
function manual_binder_parameter(ψ::MPS, sites; maxdim=64, cutoff=1e-10, chunk_size=1000)
    L = length(sites)
    
    println("Computing manual Binder parameter for L=$L with maxdim=$maxdim")
    
    # Ensure MPS is normalized and well-conditioned
    ψ_work = copy(ψ)
    orthogonalize!(ψ_work, L÷2)
    normalize!(ψ_work)
    
    # Compute all 2-point correlators ⟨Sz_i Sz_j⟩
    println("Computing 2-point correlators...")
    sum2_sq = 0.0
    total_2point = L^2
    
    for i in 1:L, j in 1:L
        try
            c2_val = manual_correlator_2point(ψ_work, sites, i, j; maxdim=maxdim, cutoff=cutoff)
            sum2_sq += c2_val^2
        catch e
            println("Warning: 2-point correlator ($i,$j) failed: $e")
        end
    end
    
    # Compute 4-point correlators in chunks to manage memory
    println("Computing 4-point correlators...")
    sum4_sq = 0.0
    total_4point = L^4
    processed = 0
    
    # Process in chunks to avoid memory issues
    for i in 1:L, j in 1:L
        chunk_4point = Vector{Float64}()
        
        for k in 1:L, l in 1:L
            try
                c4_val = manual_correlator_4point(ψ_work, sites, i, j, k, l; 
                                                maxdim=maxdim, cutoff=cutoff)
                sum4_sq += c4_val^2
                processed += 1
                
                # Progress indicator
                if processed % chunk_size == 0
                    progress = 100 * processed / total_4point
                    println("  Progress: $(round(progress, digits=1))% ($processed/$total_4point)")
                end
                
            catch e
                println("Warning: 4-point correlator ($i,$j,$k,$l) failed: $e")
                processed += 1
            end
        end
    end
    
    println("Completed correlator calculations.")
    
    # Apply proper normalization like the notebook
    M2sq = sum2_sq / L^2
    M4sq = sum4_sq / L^4
    
    println("M2sq = $M2sq, M4sq = $M4sq")
    
    # Compute Binder parameter
    denominator = 3.0 * max(M2sq^2, 1e-30)
    B = 1.0 - M4sq / denominator
    
    println("Binder parameter B = $B")
    
    return (B = B, S2 = sum2_sq, S4 = sum4_sq, M2sq = M2sq, M4sq = M4sq)
end

"""
    chunked_manual_binder_parameter(ψ::MPS, sites; maxdim=64, cutoff=1e-10, 
                                   max_4point_chunk=100)

Memory-efficient version that processes 4-point correlators in smaller chunks.
Recommended for L ≥ 16.
"""
function chunked_manual_binder_parameter(ψ::MPS, sites; maxdim=64, cutoff=1e-10, 
                                       max_4point_chunk=100)
    L = length(sites)
    
    println("Computing chunked manual Binder parameter for L=$L")
    
    # Ensure MPS is normalized and well-conditioned
    ψ_work = copy(ψ)
    orthogonalize!(ψ_work, L÷2)
    normalize!(ψ_work)
    
    # Compute 2-point correlators
    println("Computing 2-point correlators...")
    sum2_sq = 0.0
    for i in 1:L, j in 1:L
        try
            c2_val = manual_correlator_2point(ψ_work, sites, i, j; maxdim=maxdim, cutoff=cutoff)
            sum2_sq += c2_val^2
        catch e
            println("Warning: 2-point correlator ($i,$j) failed: $e")
        end
    end
    
    # Compute 4-point correlators in very small chunks
    println("Computing 4-point correlators in chunks...")
    sum4_sq = 0.0
    total_processed = 0
    total_4point = L^4
    
    # Process 4-point correlators in manageable chunks
    for chunk_start in 1:max_4point_chunk:L^4
        chunk_end = min(chunk_start + max_4point_chunk - 1, L^4)
        
        for linear_idx in chunk_start:chunk_end
            # Convert linear index to (i,j,k,l)
            temp = linear_idx - 1
            i = (temp ÷ (L^3)) + 1
            temp = temp % (L^3)
            j = (temp ÷ (L^2)) + 1
            temp = temp % (L^2)
            k = (temp ÷ L) + 1
            l = (temp % L) + 1
            
            try
                c4_val = manual_correlator_4point(ψ_work, sites, i, j, k, l; 
                                                maxdim=maxdim, cutoff=cutoff)
                sum4_sq += c4_val^2
                
            catch e
                println("Warning: 4-point correlator ($i,$j,$k,$l) failed: $e")
            end
            
            total_processed += 1
        end
        
        # Progress update
        progress = 100 * total_processed / total_4point
        println("  Chunk complete. Progress: $(round(progress, digits=1))% ($total_processed/$total_4point)")
        
        # Force garbage collection between chunks
        Base.GC.gc()
    end
    
    # Apply proper normalization like the notebook
    M2sq = sum2_sq / L^2
    M4sq = sum4_sq / L^4
    
    # Compute Binder parameter
    denominator = 3.0 * max(M2sq^2, 1e-30)
    B = 1.0 - M4sq / denominator
    
    println("Final results: M2sq = $M2sq, M4sq = $M4sq, B = $B")
    
    return (B = B, S2 = sum2_sq, S4 = sum4_sq, M2sq = M2sq, M4sq = M4sq)
end

end
