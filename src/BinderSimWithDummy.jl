"""
BinderSimWithDummy.jl

Implementation of quantum trajectory simulation with a dummy site at the end.
The system has L+1 sites total, but we only focus on the first L-1 sites for physics.
The dummy site (site L) allows ZZ measurements on all bonds including (L-1, L).
"""

module BinderSimWithDummy

using ITensors, ITensorMPS
using LinearAlgebra
using Printf
using Statistics
using Random

export ea_binder_mc_with_dummy

function create_weak_measurement_operators_with_dummy(sites, lambda_x::Float64, lambda_zz::Float64)
    """
    Create weak measurement operators for system with dummy site.
    sites: Index array of length L+1 (includes dummy site)
    Only sites 1:(L-1) participate in X measurements
    ZZ measurements occur on bonds 1:(L-1) → includes (L-1,L) with dummy
    """
    L_total = length(sites)
    L_physical = L_total - 1  # Exclude dummy site from physical system
    
    WEAK_X_0 = Dict{Int,ITensor}()
    WEAK_X_1 = Dict{Int,ITensor}()
    WEAK_ZZ_0 = Dict{Tuple{Int,Int},ITensor}()
    WEAK_ZZ_1 = Dict{Tuple{Int,Int},ITensor}()

    norm_x  = sqrt(2 * (1 + lambda_x^2))
    norm_zz = sqrt(2 * (1 + lambda_zz^2))

    # X measurements only on physical sites (1 to L-1)
    for i in 1:L_physical
        Id_i = op("Id", sites[i])
        X_i  = 2 * op("Sx", sites[i])  # ITensor Sx = 0.5 * Pauli_X
        WEAK_X_0[i] = (Id_i + lambda_x * X_i) / norm_x     # outcome = 0
        WEAK_X_1[i] = (Id_i - lambda_x * X_i) / norm_x     # outcome = 1
    end
    
    # ZZ measurements on ALL bonds including dummy: (1,2), (2,3), ..., (L-1,L)
    for i in 1:(L_total-1)  # This includes the bond (L-1,L) with dummy site
        Z_i = 2 * op("Sz", sites[i])   # ITensor Sz = 0.5 * Pauli_Z
        Z_j = 2 * op("Sz", sites[i+1])
        II  = op("Id", sites[i]) * op("Id", sites[i+1])
        ZZ  = Z_i * Z_j
        WEAK_ZZ_0[(i,i+1)] = (II + lambda_zz * ZZ) / norm_zz   # outcome = 0
        WEAK_ZZ_1[(i,i+1)] = (II - lambda_zz * ZZ) / norm_zz   # outcome = 1
    end
    
    println("Created operators for L_physical=$L_physical, L_total=$L_total")
    println("X measurements on sites: 1:$L_physical")
    println("ZZ measurements on bonds: 1:$(L_total-1) (includes dummy)")
    
    return WEAK_X_0, WEAK_X_1, WEAK_ZZ_0, WEAK_ZZ_1
end

function create_initial_state_with_dummy(L_physical::Int)
    """
    Create initial state with L_physical + 1 sites (including dummy)
    All sites initialized to |↑⟩ state
    """
    L_total = L_physical + 1
    sites = siteinds("S=1/2", L_total)
    ψ = productMPS(sites, fill("Up", L_total))
    
    println("Created initial state: L_physical=$L_physical, L_total=$L_total")
    return ψ, sites
end

function sample_and_apply_with_dummy(ψ::MPS, K0::ITensor, K1::ITensor, which::Vector{Int};
                                    maxdim::Int=256, cutoff::Float64=1e-12)
    """Same sampling function but works with dummy site system"""
    # Calculate probabilities
    ϕ0 = product(K0, ψ, which; maxdim=maxdim, cutoff=cutoff)
    ϕ1 = product(K1, ψ, which; maxdim=maxdim, cutoff=cutoff)
    
    p0 = real(inner(ϕ0, ϕ0))
    p1 = real(inner(ϕ1, ϕ1))
    
    # Normalize probabilities
    total_p = p0 + p1
    p0 /= total_p
    p1 /= total_p
    
    # Sample outcome
    outcome = rand() < p0 ? 0 : 1
    
    # Apply corresponding operator
    if outcome == 0
        ϕ = ϕ0
        orthogonalize!(ϕ, which[end])
        normalize!(ϕ)
        return ϕ, outcome
    else
        ϕ = ϕ1
        orthogonalize!(ϕ, which[end])
        normalize!(ϕ)
        return ϕ, outcome
    end
end

function calculate_correlations_with_dummy(ψ, sites, L_physical::Int)
    """
    Calculate correlations only for the physical sites (1 to L_physical)
    Ignore the dummy site in correlation calculations
    """
    # Only use physical sites for correlation calculation
    corr_matrix = correlation_matrix(ψ, "Sz", "Sz")
    
    # Sum correlations only over physical sites
    total_corr = 0.0
    count = 0
    for i in 1:L_physical
        for j in 1:L_physical
            if i != j
                total_corr += real(corr_matrix[i, j])
                count += 1
            end
        end
    end
    
    avg_correlation = total_corr / count
    return avg_correlation
end

function ea_binder_mc_with_dummy(L_physical::Int; 
                                lambda_x::Float64=0.3, 
                                lambda_zz::Float64=0.7, 
                                T_max::Int=20,
                                maxdim::Int=256,
                                cutoff::Float64=1e-12,
                                seed::Union{Int,Nothing}=nothing)
    """
    Edwards-Anderson Binder parameter simulation with dummy site
    
    Args:
        L_physical: Number of physical sites (actual system size)
        lambda_x: X measurement strength
        lambda_zz: ZZ measurement strength  
        T_max: Number of timesteps
        maxdim: Maximum bond dimension
        cutoff: SVD cutoff
        seed: Random seed
        
    Returns:
        final_correlation: Average correlation over physical sites only
    """
    
    if seed !== nothing
        Random.seed!(seed)
    end
    
    println("=" ^ 60)
    println("EDWARDS-ANDERSON SIMULATION WITH DUMMY SITE")
    println("L_physical = $L_physical (+ 1 dummy site)")
    println("λₓ = $lambda_x, λ_zz = $lambda_zz")
    println("T_max = $T_max, maxdim = $maxdim")
    println("=" ^ 60)
    
    # Initialize system with dummy site
    ψ, sites = create_initial_state_with_dummy(L_physical)
    L_total = length(sites)
    
    # Create measurement operators
    KX0, KX1, KZZ0, KZZ1 = create_weak_measurement_operators_with_dummy(sites, lambda_x, lambda_zz)
    
    # Evolution loop
    for t in 1:T_max
        # X measurements on physical sites only (1 to L_physical)
        for i in 1:L_physical
            ψ, outcome_x = sample_and_apply_with_dummy(ψ, KX0[i], KX1[i], [i]; 
                                                      maxdim=maxdim, cutoff=cutoff)
        end
        
        # ZZ measurements on ALL bonds (including dummy): (1,2), (2,3), ..., (L_physical, L_total)
        for i in 1:(L_total-1)
            bond = (i, i+1)
            ψ, outcome_zz = sample_and_apply_with_dummy(ψ, KZZ0[bond], KZZ1[bond], [i, i+1]; 
                                                       maxdim=maxdim, cutoff=cutoff)
        end
        
        # Monitor bond dimension
        current_maxdim = maxlinkdim(ψ)
        if t % (T_max ÷ 4) == 0 || t == T_max
            correlation = calculate_correlations_with_dummy(ψ, sites, L_physical)
            println("Step $t/$T_max: maxdim = $current_maxdim, avg_corr = $(round(correlation, digits=6))")
        end
    end
    
    # Final correlation calculation (physical sites only)
    final_correlation = calculate_correlations_with_dummy(ψ, sites, L_physical)
    
    println("\nFINAL RESULTS:")
    println("Physical system size: $L_physical sites")
    println("Total system size: $L_total sites (including dummy)")
    println("Final average correlation: $(round(final_correlation, digits=8))")
    println("ZZ bonds measured: $(L_total-1) (includes dummy bond)")
    
    return final_correlation
end

# Test function to compare with and without dummy
function test_dummy_effect()
    """Test the effect of adding dummy site"""
    
    L_test = 6
    
    println("\nTESTING DUMMY SITE EFFECT")
    println("=" ^ 50)
    
    # Test with dummy site
    println("\nWith dummy site:")
    corr_with_dummy = ea_binder_mc_with_dummy(L_test; seed=12345)
    
    # For comparison, we could implement the standard version
    # But for now, just show the dummy version results
    
    println("\nSummary:")
    println("L_physical = $L_test")
    println("L_total = $(L_test + 1) (with dummy)")
    println("Average correlation = $(round(corr_with_dummy, digits=8))")
    
    return corr_with_dummy
end

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    # Run test when script is executed directly
    test_dummy_effect()
end

end  # module BinderSimWithDummy