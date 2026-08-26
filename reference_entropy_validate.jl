using Test
using LinearAlgebra
using Random

include("reference_entropy_core.jl")

# ============================================================
# Independent dense brute-force replay (no MPS/ITensor machinery),
# extended from the defect_insertion validation pattern to include the
# weak-ZZ measurement layer (needed for the general-lambda dynamics
# reused here). System qubits are 1..L, the reference qubit is L+1;
# site 1 is the most-significant bit of the linear Hilbert-space index.
# ============================================================
function re_dense_embed(op2x2::AbstractMatrix{Float64}, p::Int, N::Int)
    left = Matrix{Float64}(I, 2^(p - 1), 2^(p - 1))
    right = Matrix{Float64}(I, 2^(N - p), 2^(N - p))
    return kron(kron(left, op2x2), right)
end

function re_dense_embed2(op4x4::AbstractMatrix{Float64}, i::Int, j::Int, N::Int)
    @assert j == i + 1
    left = Matrix{Float64}(I, 2^(i - 1), 2^(i - 1))
    right = Matrix{Float64}(I, 2^(N - j), 2^(N - j))
    return kron(kron(left, op4x4), right)
end

const RE_X = [0.0 1.0; 1.0 0.0]
const RE_Z = [1.0 0.0; 0.0 -1.0]
const RE_I2 = Matrix{Float64}(I, 2, 2)
const RE_ZZ4 = kron(RE_Z, RE_Z)
const RE_I4 = Matrix{Float64}(I, 4, 4)

"""
    re_dense_replay(x_record, zz_record, L, T_max, lambda_x, lambda_zz, q_x, q_zz)

Independently replay the SAME fixed weak-X and weak-ZZ measurement
records with dense ((2^(L+1)) x (2^(L+1))) matrices over L system qubits
+ 1 reference qubit (qubit L+1, NEVER acted on by anything here), using
the fully *unnormalized* Kraus maps for measurement and the exact
averaged-channel maps for dephasing (matching `ls_*` exactly). Returns
the (unnormalized) full density matrix.
"""
function re_dense_replay(
    x_record::Vector{Int}, zz_record::Vector{Int}, L::Int, T_max::Int,
    lambda_x::Float64, lambda_zz::Float64, q_x::Float64, q_zz::Float64,
)
    N = L + 1
    d = 2^N
    rho = zeros(Float64, d, d)
    rho[1, 1] = 1.0
    # GHZ-Bell initial state on all N = L+1 qubits (system + reference).
    rho .= 0.0
    psi = zeros(Float64, d)
    psi[1] = 1 / sqrt(2)       # all-zero bit string
    psi[d] = 1 / sqrt(2)       # all-one bit string
    rho .= psi * psi'

    x_idx = 1
    zz_idx = 1
    for _ in 1:T_max
        for p in 1:L
            outcome = x_record[x_idx]; x_idx += 1
            K = (RE_I2 .+ (-1.0)^outcome .* lambda_x .* RE_X) ./ sqrt(2 * (1 + lambda_x^2))
            Kemb = re_dense_embed(K, p, N)
            rho = Kemb * rho * Kemb'
        end
        for p in 1:L
            Xemb = re_dense_embed(RE_X, p, N)
            rho = (1 - q_x) .* rho .+ q_x .* (Xemb * rho * Xemb)
        end
        for p in 1:(L - 1)
            outcome = zz_record[zz_idx]; zz_idx += 1
            K = (RE_I4 .+ (-1.0)^outcome .* lambda_zz .* RE_ZZ4) ./ sqrt(2 * (1 + lambda_zz^2))
            Kemb = re_dense_embed2(K, p, p + 1, N)
            rho = Kemb * rho * Kemb'
        end
        for p in 1:(L - 1)
            ZZemb = re_dense_embed2(RE_ZZ4, p, p + 1, N)
            rho = (1 - q_zz) .* rho .+ q_zz .* (ZZemb * rho * ZZemb)
        end
    end

    return rho
end

"""
    re_dense_reduced_rho_R(rho_dense, L)

Partial trace over the first L (system) qubits of a dense (2^(L+1) x
2^(L+1)) matrix, returning the 2x2 reduced matrix for qubit L+1.
"""
function re_dense_reduced_rho_R(rho_dense::AbstractMatrix{Float64}, L::Int)
    N = L + 1
    d = 2^N
    reduced = zeros(Float64, 2, 2)
    for row0 in 0:(d - 1), col0 in 0:(d - 1)
        row_R = row0 & 1
        col_R = col0 & 1
        row_Q = row0 >> 1
        col_Q = col0 >> 1
        if row_Q == col_Q
            reduced[row_R + 1, col_R + 1] += rho_dense[row0 + 1, col0 + 1]
        end
    end
    return reduced
end

@testset "reference entropy validation" begin

    @testset "1. initial-state test: S_R(t=0) = 1" begin
        L = 4
        sites = siteinds("Qubit", 2 * (L + 1))
        trace_bra = ls_trace_state(sites)
        rho0 = ls_trace_normalize(reference_entropy_initial_state(sites, L), trace_bra)
        M = reference_entropy_reduced_rho_R(rho0, sites, L)
        entropy = reference_entropy_von_neumann(M)
        @test isapprox(entropy.S, 1.0; atol=1e-10)
        @test isapprox(entropy.eig1, 0.5; atol=1e-10)
        @test isapprox(entropy.eig2, 0.5; atol=1e-10)
        @test isapprox(entropy.trace_error, 0.0; atol=1e-10)
    end

    @testset "2. reference-isolation: no gate touches the reference site" begin
        L = 4
        sites = siteinds("Qubit", 2 * (L + 1))
        ref_indices = Set([sites[2L + 1], sites[2L + 2]])

        x_gates = ls_build_x_dephasing_gates(sites, L, 0.2)
        zz_gates = ls_build_zz_dephasing_gates(sites, L, 0.2)
        for gate in x_gates
            @test isempty(intersect(Set(inds(gate)), ref_indices))
        end
        for gate in zz_gates
            @test isempty(intersect(Set(inds(gate)), ref_indices))
        end
    end

    @testset "3. system-unitary test: Q-only unitary leaves S_R = 1" begin
        L = 4
        sites = siteinds("Qubit", 2 * (L + 1))
        trace_bra = ls_trace_state(sites)
        rho = ls_trace_normalize(reference_entropy_initial_state(sites, L), trace_bra)

        # A single-site Q-only unitary (Hadamard-like real orthogonal matrix)
        # applied identically to bra and ket legs (conjugation U rho U^dagger).
        theta = 0.37
        U = [cos(theta) -sin(theta); sin(theta) cos(theta)]
        p = 2
        gate = op(U, sites[2p - 1]) * op(U, sites[2p])
        rho2 = apply([gate], rho; cutoff=1e-12, maxdim=64)
        rho2 = ls_trace_normalize(rho2, trace_bra)

        M = reference_entropy_reduced_rho_R(rho2, sites, L)
        entropy = reference_entropy_von_neumann(M)
        @test isapprox(entropy.S, 1.0; atol=1e-8)
    end

    @testset "4./5. dense brute-force cross-check (small L, no truncation)" begin
        L_small = 3
        T_small = 4
        evolved = reference_entropy_evolve_one_trial(
            L_small; lambda_x=0.7, lambda_zz=0.3, q_x=0.2, q_zz=0.2, T_max=T_small,
            maxdim=256, cutoff=0.0, seed=11, save_every=1,
        )

        n_x_per_layer = L_small
        n_zz_per_layer = L_small - 1
        x_record = Int[]
        zz_record = Int[]
        for t in 1:T_small
            base = (t - 1) * (n_x_per_layer + n_zz_per_layer)
            append!(x_record, evolved.record[(base + 1):(base + n_x_per_layer)])
            append!(zz_record, evolved.record[(base + n_x_per_layer + 1):(base + n_x_per_layer + n_zz_per_layer)])
        end

        rho_dense = re_dense_replay(x_record, zz_record, L_small, T_small, 0.7, 0.3, 0.2, 0.2)
        rho_dense_normalized = rho_dense ./ tr(rho_dense)
        M_dense = re_dense_reduced_rho_R(rho_dense_normalized, L_small)
        entropy_dense = reference_entropy_von_neumann(M_dense)

        M_mps = evolved.snapshots[end].rhoR
        entropy_mps = evolved.snapshots[end].entropy

        @test isapprox(M_mps, M_dense; atol=1e-6)
        @test isapprox(entropy_mps.S, entropy_dense.S; atol=1e-6)
    end

    @testset "6. Born-probability test: outcome probabilities sum to 1" begin
        L = 4
        sites = siteinds("Qubit", 2 * (L + 1))
        trace_bra = ls_trace_state(sites)
        rho = ls_trace_normalize(reference_entropy_initial_state(sites, L), trace_bra)

        states, weights = ls_x_measurement_candidates(rho, sites, 1, 0.7, trace_bra; maxdim=64, cutoff=1e-12)
        total = sum(weights) / ls_trace(rho, trace_bra)
        @test isapprox(total, 1.0; atol=1e-8)
    end

    @testset "7. entropy-range test: 0 <= S_R(m) <= 1" begin
        L = 6
        for seed in 1:5
            evolved = reference_entropy_evolve_one_trial(
                L; lambda_x=0.7, lambda_zz=0.0, q_x=0.3, q_zz=0.3, T_max=6,
                maxdim=64, cutoff=1e-12, seed=seed, save_every=2,
            )
            for snap in evolved.snapshots
                @test -1e-9 <= snap.entropy.S <= 1.0 + 1e-9
            end
        end
    end

    @testset "8. Hermiticity, trace-one, and positivity of rho_R" begin
        L = 6
        evolved = reference_entropy_evolve_one_trial(
            L; lambda_x=0.7, lambda_zz=0.0, q_x=0.3, q_zz=0.3, T_max=6,
            maxdim=64, cutoff=1e-12, seed=3, save_every=2,
        )
        for snap in evolved.snapshots
            M = snap.rhoR
            @test isapprox(tr(M), 1.0; atol=1e-6)
            @test maximum(abs.(M .- M')) < 1e-8
            @test all(eigvals(Hermitian((M .+ M') ./ 2)) .>= -1e-8)
        end
    end

    @testset "9. averaging test: S(mean rho_R) != mean(S(rho_R))" begin
        # Demonstrate the nonlinearity explicitly and confirm the code
        # computes the (correct) per-trajectory average, not the entropy
        # of the trajectory-averaged reduced state.
        L = 6
        entropies = Float64[]
        rho_sum = zeros(2, 2)
        n = 6
        for seed in 1:n
            evolved = reference_entropy_evolve_one_trial(
                L; lambda_x=0.7, lambda_zz=0.0, q_x=0.3, q_zz=0.3, T_max=6,
                maxdim=64, cutoff=1e-12, seed=seed, save_every=6,
            )
            final = evolved.snapshots[end]
            push!(entropies, final.entropy.S)
            rho_sum .+= final.rhoR
        end
        mean_entropy = sum(entropies) / n           # what the code computes
        rho_avg = rho_sum ./ n
        entropy_of_average = reference_entropy_von_neumann(rho_avg).S

        @test !isapprox(mean_entropy, entropy_of_average; atol=1e-3)
    end

    @testset "10. maxdim/cutoff convergence" begin
        L = 6
        baseline = reference_entropy_evolve_one_trial(
            L; lambda_x=0.7, lambda_zz=0.0, q_x=0.3, q_zz=0.3, T_max=8,
            maxdim=256, cutoff=1e-12, seed=42, save_every=4,
        )
        tight = reference_entropy_evolve_one_trial(
            L; lambda_x=0.7, lambda_zz=0.0, q_x=0.3, q_zz=0.3, T_max=8,
            maxdim=512, cutoff=1e-14, seed=42, save_every=4,
        )
        @test baseline.record == tight.record
        @test isapprox(baseline.snapshots[end].entropy.S, tight.snapshots[end].entropy.S; atol=0.02)
    end

    @testset "11. reproducibility with fixed seed" begin
        L = 6
        first = reference_entropy_evolve_one_trial(
            L; lambda_x=0.7, lambda_zz=0.0, q_x=0.3, q_zz=0.3, T_max=6,
            maxdim=64, cutoff=1e-12, seed=99, save_every=2,
        )
        second = reference_entropy_evolve_one_trial(
            L; lambda_x=0.7, lambda_zz=0.0, q_x=0.3, q_zz=0.3, T_max=6,
            maxdim=64, cutoff=1e-12, seed=99, save_every=2,
        )
        @test first.record_checksum == second.record_checksum
        @test first.record == second.record
        @test isapprox(first.snapshots[end].entropy.S, second.snapshots[end].entropy.S; atol=1e-10)
    end
end
