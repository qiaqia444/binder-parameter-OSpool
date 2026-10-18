#!/usr/bin/env julia

"""
Validation / smoke test for the GHZ-initial-state right-boundary dual-EA Binder observable
(`dual_ea_right_boundary_ghz_core.jl`).

    julia --project=. dual_ea_right_boundary_ghz_validate.jl

Checks (small L; dense matrices are used ONLY here as the independent reference):
  1. The doubled-MPS initial state equals the dense |GHZ><GHZ| (trace 1, exact entries).
  2. At t = 0 every non-empty dual X-string has zero expectation, so M2 = 1/N,
     M4 = (3N^2-2N)/N^4, B = 2/(3N).
  3. MPS moments vs. brute-force dense calculation over ALL ordered tuples on genuine evolved
     GHZ-start trajectories (operators W_ab W_cd multiplied out explicitly).
  4. Incremental DP vs. brute-force MPS chain evaluation over all ordered tuples at larger L.
  5. Exact q = 1/2 limit for EVERY trajectory: M2 = 1/N, M4 = (3N^2-2N)/N^4, B = 2/(3N).
  6. q = 0 and T = 4L: the weak X measurements collapse the GHZ state onto an X-eigenstate
     product state, so M2 = M4 = 1 and B = 2/3.
  7. The GHZ evolution differs from the all-up evolution only through the initial state
     (checked by evolving the all-up state with the same code path and comparing to
     `rb_evolve_right_edge_one_trial`), plus non-mutation and determinism.
  8. Seed derivation identical to `rb_run_right_edge_point` / the Renyi-2 GHZ driver.
"""

using Test
using Random
using LinearAlgebra
using Statistics
using ITensors
using ITensorMPS

include("dual_ea_right_boundary_ghz_core.jl")

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

evolve_ghz(L, q, seed; T_max=2L, maxdim=256, cutoff=1e-14) = dual_ea_ghz_evolve_right_edge_one_trial(
    L; lambda_x=0.7, q=q, T_max=T_max, maxdim=maxdim, cutoff=cutoff, seed=seed,
)

mixed_limits(N) = (M2=1 / N, M4=(3N^2 - 2N) / N^4, B=2 / (3N))

@testset "dual EA right-boundary GHZ-start observable" begin
    @testset "initial state is |GHZ><GHZ|" begin
        for L in (2, 3, 4, 5)
            sites = siteinds("Qubit", 2L)
            ρ = dense_from_doubled_mps(dual_ea_ghz_initial_state(sites), L)
            ψ = zeros(2^L); ψ[1] = 1 / sqrt(2); ψ[end] = 1 / sqrt(2)
            @test isapprox(ρ, ψ * ψ'; atol=1e-12)
            @test real(tr(ρ)) ≈ 1.0 atol = 1e-12
        end
    end

    @testset "t = 0: all non-empty dual strings vanish (M2 = 1/N)" begin
        for L in (3, 4, 6, 10, 20)
            N = L - 1
            sites = siteinds("Qubit", 2L)
            res = dual_ea_moments(dual_ea_ghz_initial_state(sites), L)
            ex = mixed_limits(N)
            @test res.M2 ≈ ex.M2 atol = 1e-12
            @test res.M4 ≈ ex.M4 atol = 1e-12
            @test res.B ≈ ex.B atol = 1e-12
        end
    end

    @testset "dense reference: evolved GHZ-start trajectories" begin
        for (L, q, seed) in ((4, 0.0, 1), (5, 0.2, 2), (6, 0.35, 3), (5, 0.5, 4))
            ev = evolve_ghz(L, q, seed)
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
            ev = evolve_ghz(L, q, seed; maxdim=128, cutoff=1e-12)
            ref = brute_force_mps_moments(ev.rho, L)
            res = dual_ea_moments(ev.rho, L)
            @test res.M2 ≈ ref.M2 rtol = 1e-10
            @test res.M4 ≈ ref.M4 rtol = 1e-10
        end
    end

    @testset "exact q = 1/2 limit (every trajectory)" begin
        for (L, seed) in ((4, 11), (6, 12), (8, 13))
            ev = evolve_ghz(L, 0.5, seed; maxdim=128, cutoff=1e-12)
            res = dual_ea_moments(ev.rho, L)
            ex = mixed_limits(L - 1)
            @test res.M2 ≈ ex.M2 atol = 1e-10
            @test res.M4 ≈ ex.M4 atol = 1e-10
            @test res.B ≈ ex.B atol = 1e-9
        end
    end

    @testset "q = 0, T = 4L: collapse onto an X-eigenstate product state" begin
        for (L, seed) in ((4, 21), (6, 22), (8, 23))
            ev = evolve_ghz(L, 0.0, seed; T_max=4L, maxdim=128, cutoff=1e-14)
            res = dual_ea_moments(ev.rho, L)
            @test res.M2 ≈ 1.0 atol = 1e-6
            @test res.M4 ≈ 1.0 atol = 1e-6
            @test res.B ≈ 2 / 3 atol = 1e-6
        end
    end

    @testset "same dynamics as the all-up evolution except the initial state" begin
        # Run the GHZ code path with the all-up initial state swapped in: it must reproduce
        # rb_evolve_right_edge_one_trial exactly, so GHZ vs all-up differ ONLY by the initial state.
        L, q, seed = 5, 0.3, 31
        ev_up = rb_evolve_right_edge_one_trial(
            L; lambda_x=0.7, q=q, T_max=2L, maxdim=256, cutoff=1e-14, seed=seed,
        )
        ev_swapped = dual_ea_ghz_evolve_right_edge_one_trial(
            L; lambda_x=0.7, q=q, T_max=2L, maxdim=256, cutoff=1e-14, seed=seed,
            initial_state=rb_initial_repetition_code_state,
        )
        @test dense_from_doubled_mps(ev_swapped.rho, L) ≈ dense_from_doubled_mps(ev_up.rho, L) atol = 1e-12

        # The GHZ start is NOT equivalent to all-up for this observable: products of weak-measurement
        # factors generate the full X^(x)L string, whose GHZ expectation is 1, so the GHZ parity
        # structure feeds into the Born statistics. (Differences are trajectory dependent.)
        dM2 = maximum(1:6) do s
            up = rb_evolve_right_edge_one_trial(8; lambda_x=0.7, q=0.3, T_max=32, maxdim=256, cutoff=1e-12, seed=s)
            gh = evolve_ghz(8, 0.3, s; T_max=32, maxdim=256, cutoff=1e-12)
            abs(dual_ea_moments(up.rho, 8).M2 - dual_ea_moments(gh.rho, 8).M2)
        end
        @test dM2 > 1e-4
    end

    @testset "no mutation, no RNG consumption, deterministic dynamics" begin
        L = 6
        ev = evolve_ghz(L, 0.3, 77; maxdim=64, cutoff=1e-12)
        before = mps_fingerprint(ev.rho)
        links_before = [dim(linkind(ev.rho, k)) for k in 1:(2L - 1)]

        Random.seed!(2024)
        expected_next = rand()
        Random.seed!(2024)
        r1 = dual_ea_moments(ev.rho, L)
        observed_next = rand()
        r2 = dual_ea_moments(ev.rho, L)

        @test observed_next == expected_next           # the observable draws no random numbers
        @test r1 == r2
        @test mps_fingerprint(ev.rho) == before
        @test [dim(linkind(ev.rho, k)) for k in 1:(2L - 1)] == links_before

        ev2 = evolve_ghz(L, 0.3, 77; maxdim=64, cutoff=1e-12)
        @test mps_fingerprint(ev2.rho) == before
    end

    @testset "driver reuses the Renyi-2 seed derivation" begin
        seed = 1234
        res = dual_ea_ghz_run_right_edge_point(
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
for (L, chi) in ((16, 64), (32, 128), (40, 256))
    c = dual_ea_cost_estimate(L, chi)
    println(
        "  L=$L chi=$chi: N=$(c.n_endpoints), sorted quads=$(c.n_sorted_quads) ",
        "(ordered $(c.n_ordered_quads)), ~$(round(c.approx_flops / 1e9, digits=2)) GFLOP, ",
        "~$(round(c.approx_memory_bytes / 2^20, digits=1)) MiB",
    )
end
