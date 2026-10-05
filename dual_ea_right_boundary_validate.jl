#!/usr/bin/env julia

"""
Validation / smoke test for the right-boundary dual-EA Binder observable
(`dual_ea_right_boundary_core.jl`).

    julia --project=. dual_ea_right_boundary_validate.jl

Checks (all small L, dense matrices are used ONLY here as the independent reference):
  1. MPS moments vs. brute-force dense calculation over ALL ordered tuples with the
     operators W_ab W_cd multiplied out explicitly (no sorting/multiplicity shortcut):
       real symmetric mixed states with off-diagonals, complex Hermitian states, an
       unnormalised state, and genuine evolved right-boundary trajectories.
  2. Incremental DP vs. brute-force MPS chain evaluation over all ordered tuples at larger L.
  3. Analytic limits (endpoints 2:L): |+x> product state, maximally mixed state.
  4. Non-mutation of the state, no global-RNG consumption, determinism of the dynamics.
  5. Seed derivation identical to the Renyi-2 right-boundary driver.
"""

using Test
using Random
using LinearAlgebra
using Statistics
using ITensors
using ITensorMPS

include("dual_ea_right_boundary_core.jl")

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

function doubled_mps_from_dense(rho::AbstractMatrix, L::Int, sites)
    A = reshape(Array(rho), ntuple(_ -> 2, 2L))                    # dims (bL..b1,kL..k1)
    # target order (b1,k1,b2,k2,...): b_i is dim (L-i+1), k_i is dim (L+(L-i+1))
    perm = Int[]
    for i in 1:L
        push!(perm, L - i + 1)
        push!(perm, L + (L - i + 1))
    end
    A = permutedims(A, perm)
    return MPS(ITensor(A, sites...), sites)
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

function random_density(L::Int; complex_valued::Bool=false, rank::Int=5, rng=MersenneTwister(11))
    dim = 2^L
    ρ = zeros(complex_valued ? ComplexF64 : Float64, dim, dim)
    w = rand(rng, rank)
    w ./= sum(w)
    for r in 1:rank
        ψ = complex_valued ? randn(rng, ComplexF64, dim) : randn(rng, dim)
        ψ ./= norm(ψ)
        ρ .+= w[r] .* (ψ * ψ')
    end
    return ρ
end

product_plus_x(sites, L) = begin
    ts = ITensor[]
    for i in 1:L
        p = MPS(0.5 * [1.0, 1.0, 1.0, 1.0], [sites[2i - 1], sites[2i]])
        push!(ts, p[1], p[2])
    end
    MPS(ts)
end

maximally_mixed(sites, L) = begin
    ts = ITensor[]
    for i in 1:L
        p = MPS(0.5 * [1.0, 0.0, 0.0, 1.0], [sites[2i - 1], sites[2i]])
        push!(ts, p[1], p[2])
    end
    MPS(ts)
end

mps_fingerprint(rho::MPS) = [Array(rho[k], inds(rho[k])...) for k in 1:length(rho)]

@testset "dual EA right-boundary observable" begin
    @testset "multiplicity bookkeeping" begin
        for N in 1:7
            total = 0
            for u in 1:N, v in u:N, w in v:N, x in w:N
                total += dual_ea_multiplicity4(u, v, w, x)
            end
            @test total == N^4
        end
    end

    @testset "dense reference: random mixed states (off-diagonals, repeats)" begin
        for (L, cplx) in ((2, false), (2, true), (3, false), (3, true), (4, false), (4, true), (5, false), (5, true), (6, false))
            sites = siteinds("Qubit", 2L)
            ρ = random_density(L; complex_valued=cplx, rng=MersenneTwister(100L + cplx))
            @test maximum(abs.(ρ - ρ')) < 1e-14
            offdiag = maximum(abs.(ρ - Diagonal(diag(ρ))))
            @test offdiag > 1e-3
            mps = doubled_mps_from_dense(ρ, L, sites)
            @test isapprox(dense_from_doubled_mps(mps, L), ρ; atol=1e-12)   # layout round trip
            ref = dense_dual_ea_moments(ρ, L)
            res = dual_ea_moments(mps, L)
            @test res.M2 ≈ ref.M2 rtol = 1e-10
            @test res.M4 ≈ ref.M4 rtol = 1e-10
            @test res.B ≈ ref.B rtol = 1e-9
            @test res.trace ≈ 1.0 atol = 1e-12

            # unnormalised copy gives identical moments (trace normalisation only)
            scaled = dual_ea_moments(mps * 3.7, L)
            @test scaled.M2 ≈ res.M2 rtol = 1e-12
            @test scaled.M4 ≈ res.M4 rtol = 1e-12
        end
    end

    @testset "dense reference: evolved right-boundary trajectories" begin
        for (L, q, seed) in ((4, 0.0, 1), (5, 0.2, 2), (6, 0.35, 3), (5, 0.5, 4))
            ev = rb_evolve_right_edge_one_trial(
                L; lambda_x=0.7, q=q, T_max=2L, maxdim=256, cutoff=1e-14, seed=seed,
            )
            ρ = dense_from_doubled_mps(ev.rho, L)
            @test real(tr(ρ)) ≈ 1.0 atol = 1e-10
            ref = dense_dual_ea_moments(ρ, L)
            res = dual_ea_moments(ev.rho, L)
            @test res.M2 ≈ ref.M2 rtol = 1e-9
            @test res.M4 ≈ ref.M4 rtol = 1e-9
            @test res.B ≈ ref.B rtol = 1e-8
        end
    end

    @testset "DP vs brute-force MPS chains at larger L" begin
        for (L, q, seed) in ((8, 0.25, 5), (10, 0.4, 6))
            ev = rb_evolve_right_edge_one_trial(
                L; lambda_x=0.7, q=q, T_max=2L, maxdim=128, cutoff=1e-12, seed=seed,
            )
            ref = brute_force_mps_moments(ev.rho, L)
            res = dual_ea_moments(ev.rho, L)
            @test res.M2 ≈ ref.M2 rtol = 1e-10
            @test res.M4 ≈ ref.M4 rtol = 1e-10
        end
    end

    @testset "analytic limits (endpoints 2:L)" begin
        for L in (2, 3, 4, 7, 12, 30)
            sites = siteinds("Qubit", 2L)
            N = L - 1

            plus = dual_ea_moments(product_plus_x(sites, L), L)
            @test plus.M2 ≈ 1.0 atol = 1e-12
            @test plus.M4 ≈ 1.0 atol = 1e-12
            @test plus.B ≈ 2 / 3 atol = 1e-12

            mixed = dual_ea_moments(maximally_mixed(sites, L), L)
            @test mixed.M2 ≈ 1 / N atol = 1e-12
            @test mixed.M4 ≈ (3N^2 - 2N) / N^4 atol = 1e-12
            @test mixed.B ≈ 2 / (3N) atol = 1e-12
        end
    end

    @testset "no mutation, no RNG consumption, deterministic dynamics" begin
        L = 6
        ev = rb_evolve_right_edge_one_trial(
            L; lambda_x=0.7, q=0.3, T_max=2L, maxdim=64, cutoff=1e-12, seed=77,
        )
        before = mps_fingerprint(ev.rho)
        links_before = [dim(linkind(ev.rho, k)) for k in 1:(2L - 1)]

        Random.seed!(2024)
        expected_next = rand()
        Random.seed!(2024)
        r1 = dual_ea_moments(ev.rho, L)
        observed_next = rand()
        r2 = dual_ea_moments(ev.rho, L)

        @test observed_next == expected_next           # global RNG untouched
        @test r1 == r2                                 # deterministic / repeatable
        @test mps_fingerprint(ev.rho) == before        # tensors bitwise unchanged
        @test [dim(linkind(ev.rho, k)) for k in 1:(2L - 1)] == links_before

        ev2 = rb_evolve_right_edge_one_trial(
            L; lambda_x=0.7, q=0.3, T_max=2L, maxdim=64, cutoff=1e-12, seed=77,
        )
        @test mps_fingerprint(ev2.rho) == before       # measuring did not alter later dynamics
    end

    @testset "driver reuses the Renyi-2 seed derivation" begin
        seed = 1234
        res = dual_ea_run_right_edge_point(
            4, 0.3; lambda_x=0.7, lambda_zz=0.0, ntrials=3, T_max=8,
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
for (L, chi) in ((16, 64), (32, 128), (56, 256), (56, 512))
    c = dual_ea_cost_estimate(L, chi)
    println(
        "  L=$L chi=$chi: N=$(c.n_endpoints), sorted quads=$(c.n_sorted_quads) ",
        "(ordered $(c.n_ordered_quads)), ~$(round(c.approx_flops / 1e9, digits=2)) GFLOP, ",
        "~$(round(c.approx_memory_bytes / 2^20, digits=1)) MiB",
    )
end
