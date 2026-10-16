#!/usr/bin/env julia

"""
Validation / smoke test for the trivial-to-learning dual-EA Binder observable
(`dual_ea_trivial_to_learning_lambda0.7_core.jl`, lambda = 0.7: lambda_x = 0.49, lambda_zz = 0.21).

    julia --project=. dual_ea_trivial_to_learning_lambda0.7_validate.jl

Checks (small L; dense matrices are used ONLY here as the independent reference):
  1. MPS moments vs. brute-force dense calculation over ALL ordered tuples, on genuine
     evolved trajectories of the lambda_x = 0.49, lambda_zz = 0.21 dynamics (both weak
     layers active).
  2. Incremental DP vs. brute-force MPS chain evaluation over all ordered tuples at larger L.
  3. Exact q = 1/2 limit for EVERY trajectory: M2 = 1/N, M4 = (3N^2-2N)/N^4, B = 2/(3N).
  4. Non-mutation of the state and determinism of the dynamics.
  5. Seed derivation identical to `lt_run_fixed_lambda_point`.
"""

using Test
using Random
using LinearAlgebra
using Statistics
using ITensors
using ITensorMPS

include("dual_ea_trivial_to_learning_lambda0.7_core.jl")

const LAMBDA_X = 0.49
const LAMBDA_ZZ = 0.21

# ------------------------------------------------------------
# Dense helpers (reference only)
# ------------------------------------------------------------
const VX = Float64[0 1; 1 0]
const V1 = Matrix{Float64}(I, 2, 2)

kron_sites(ops) = reduce(kron, ops)   # site 1 = most significant factor

function dense_from_doubled_mps(rho::MPS, L::Int)
    sites = siteinds(rho)
    T = rho[1]
    for k in 2:(2L)
        T = T * rho[k]
    end
    A = array(T, sites...)                      # dims ordered (b1,k1,b2,k2,...)
    A = reshape(A, ntuple(_ -> 2, 2L))
    perm = vcat([2i - 1 for i in L:-1:1], [2i for i in L:-1:1])   # (bL..b1, kL..k1)
    return reshape(permutedims(A, perm), 2^L, 2^L)
end

function dense_dual_ea_moments(rho::AbstractMatrix, L::Int)
    N = L - 1
    tr_rho = real(tr(rho))
    xsite(k) = kron_sites([j == k ? VX : V1 for j in 1:L])
    X = Dict(k => xsite(k) for k in 1:L)
    Id = Matrix{Float64}(I, 2^L, 2^L)
    function W(a, b)
        lo, hi = minmax(a, b)
        O = copy(Id)
        for k in lo:(hi - 1)
            O = O * X[k]
        end
        return O
    end
    Ws = Dict((a, b) => W(a, b) for a in 2:L for b in 2:L)
    E(O) = real(tr(rho * O)) / tr_rho

    M2 = sum(E(Ws[(a, b)])^2 for a in 2:L, b in 2:L) / N^2
    M4 = 0.0
    for a in 2:L, b in 2:L, c in 2:L, d in 2:L
        M4 += E(Ws[(a, b)] * Ws[(c, d)])^2
    end
    M4 /= N^4
    return (M2=M2, M4=M4, B=1 - M4 / (3 * M2^2))
end

# Brute-force MPS evaluation over all ordered tuples (no multiplicity shortcut).
function brute_force_mps_moments(rho::MPS, L::Int)
    MI, MX = dual_ea_transfer_matrices(rho, L)
    N = L - 1
    tr_rho = dual_ea_pattern_expectation(MI, MX, falses(L))
    E(eps...) = dual_ea_pattern_expectation(MI, MX, dual_ea_string_pattern(L, eps)) / tr_rho
    M2 = sum(E(a, b)^2 for a in 2:L, b in 2:L) / N^2
    M4 = sum(E(a, b, c, d)^2 for a in 2:L, b in 2:L, c in 2:L, d in 2:L) / N^4
    return (M2=M2, M4=M4, B=1 - M4 / (3 * M2^2))
end

mps_fingerprint(rho::MPS) = [Array(rho[k], inds(rho[k])...) for k in 1:length(rho)]

evolve(L, q, seed; T_max=2L, maxdim=256, cutoff=1e-14) = lt_evolve_one_trial(
    L; lambda_x=LAMBDA_X, lambda_zz=LAMBDA_ZZ, q=q, T_max=T_max,
    maxdim=maxdim, cutoff=cutoff, seed=seed,
)

@testset "dual EA trivial-to-learning observable" begin
    @testset "dense reference: evolved trajectories (λ_x=0.49, λ_zz=0.21)" begin
        for (L, q, seed) in ((4, 0.0, 1), (5, 0.2, 2), (6, 0.3, 3), (5, 0.5, 4))
            ev = evolve(L, q, seed)
            ρ = dense_from_doubled_mps(ev.rho, L)
            @test real(tr(ρ)) ≈ 1.0 atol = 1e-10
            @test maximum(abs.(ρ - ρ')) < 1e-7   # MPS truncation in the dynamics limits exact Hermiticity
            ref = dense_dual_ea_moments(ρ, L)
            res = dual_ea_moments(ev.rho, L)
            @test res.M2 ≈ ref.M2 rtol = 1e-9
            @test res.M4 ≈ ref.M4 rtol = 1e-9
            @test res.B ≈ ref.B rtol = 1e-8
        end
    end

    @testset "DP vs brute-force MPS chains at larger L" begin
        for (L, q, seed) in ((8, 0.25, 5), (10, 0.4, 6))
            ev = evolve(L, q, seed; maxdim=128, cutoff=1e-12)
            ref = brute_force_mps_moments(ev.rho, L)
            res = dual_ea_moments(ev.rho, L)
            @test res.M2 ≈ ref.M2 rtol = 1e-10
            @test res.M4 ≈ ref.M4 rtol = 1e-10
        end
    end

    @testset "exact q = 1/2 limit (every trajectory)" begin
        for (L, seed) in ((4, 11), (6, 12), (8, 13))
            N = L - 1
            ev = evolve(L, 0.5, seed; maxdim=128, cutoff=1e-12)
            res = dual_ea_moments(ev.rho, L)
            @test res.M2 ≈ 1 / N atol = 1e-10
            @test res.M4 ≈ (3N^2 - 2N) / N^4 atol = 1e-10
            @test res.B ≈ 2 / (3N) atol = 1e-9
        end
    end

    @testset "no mutation, deterministic dynamics" begin
        L = 6
        ev = evolve(L, 0.3, 77; maxdim=64, cutoff=1e-12)
        before = mps_fingerprint(ev.rho)
        links_before = [dim(linkind(ev.rho, k)) for k in 1:(2L - 1)]

        r1 = dual_ea_moments(ev.rho, L)
        r2 = dual_ea_moments(ev.rho, L)
        @test r1 == r2
        @test mps_fingerprint(ev.rho) == before
        @test [dim(linkind(ev.rho, k)) for k in 1:(2L - 1)] == links_before

        ev2 = evolve(L, 0.3, 77; maxdim=64, cutoff=1e-12)
        @test mps_fingerprint(ev2.rho) == before
    end

    @testset "driver reuses the learning-to-trivial seed derivation" begin
        seed = 1234
        res = dual_ea_tl_run_point(
            4, 0.3; lambda_x=LAMBDA_X, lambda_zz=LAMBDA_ZZ, ntrials=3, T_max=8,
            maxdim=64, cutoff=1e-12, seed=seed, nboot=50,
        )
        rng = MersenneTwister(seed)
        @test res.trial_seeds == [Int(rand(rng, UInt32)) for _ in 1:3]
        @test res.n_completed == 3
        @test res.B_mean_of_trials ≈ mean(res.Bs)
        @test res.B_ratio_of_means ≈ 1 - mean(res.M4s) / (3 * mean(res.M2s)^2)
    end
end

println("\nCost model (exact 4-point sum, no approximation):")
for (L, chi) in ((16, 64), (32, 128), (40, 256))
    c = dual_ea_cost_estimate(L, chi)
    println(
        "  L=$L chi=$chi: N=$(c.n_endpoints), sorted quads=$(c.n_sorted_quads) ",
        "(ordered $(c.n_ordered_quads)), ~$(round(c.approx_flops / 1e9, digits=2)) GFLOP, ",
        "~$(round(c.approx_memory_bytes / 2^20, digits=1)) MiB",
    )
end
