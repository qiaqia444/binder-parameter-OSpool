#!/usr/bin/env julia

"""
Validation / smoke test for the Renyi-2 right-boundary Binder measured after the last
X dephasing (`renyi2_right_boundary_xnoise_core.jl`).

    julia --project=. renyi2_right_boundary_xnoise_validate.jl

Checks:
  1. Applying the skipped final ZZ-dephasing layer to the X-noise state reproduces the production
     final state (`rb_evolve_right_edge_one_trial`) bit-for-bit, for the same seed.
  2. Same trial seeds / same point statistics seed derivation as `rb_run_right_edge_point`.
  3. At q = 0 the final ZZ layer is the identity: both observation times coincide.
  4. The production correlator on the X-noise state agrees with an independent dense calculation.
  5. The after-ZZ moments from the production runner are recovered from the X-noise runner's
     trajectories (same trajectories, observed at two times); the observation times differ for q > 0.
  6. No mutation of the state / no consumption of the global RNG by the observable.
  7. Edge case T_max = 1.
"""

using Test
using Random
using LinearAlgebra
using Statistics
using ITensors
using ITensorMPS

include("renyi2_right_boundary_xnoise_core.jl")

mps_fingerprint(rho::MPS) = [Array(rho[k], inds(rho[k])...) for k in 1:length(rho)]

function dense_from_doubled_mps(rho::MPS, L::Int)
    sites = siteinds(rho)
    T = rho[1]
    for k in 2:(2L)
        T = T * rho[k]
    end
    A = reshape(array(T, sites...), ntuple(_ -> 2, 2L))
    perm = vcat([2i - 1 for i in L:-1:1], [2i for i in L:-1:1])
    return reshape(permutedims(A, perm), 2^L, 2^L)
end

# <<rho|Q^n|rho>> / P with Q = sum_i Z_bra,i Z_ket,i, from the dense matrix elements (independent of the MPO)
function dense_renyi2_moments(rho_dense::AbstractMatrix, L::Int)
    P = M2 = M4 = 0.0
    for b in 0:(2^L - 1), k in 0:(2^L - 1)
        w = abs2(rho_dense[b + 1, k + 1])
        Qv = L - 2 * count_ones(xor(b, k))
        P += w
        M2 += w * Qv^2
        M4 += w * Qv^4
    end
    M2 /= (L^2 * P)
    M4 /= (L^4 * P)
    return (M2=M2, M4=M4, B2=1 - M4 / (3 * M2^2), purity=P)
end

@testset "Renyi-2 right boundary, after last X dephasing" begin
    @testset "final ZZ layer applied to the X-noise state == production state" begin
        for (L, q, seed) in ((3, 0.0, 1), (4, 0.2, 2), (5, 0.425, 3), (6, 0.5, 4), (6, 0.3, 5))
            kw = (lambda_x=0.7, q=q, T_max=3L, maxdim=64, cutoff=1e-12, seed=seed)
            xn = rb_xn_evolve_right_edge_one_trial(L; kw...)
            prod = rb_evolve_right_edge_one_trial(L; kw...)

            zz = rb_build_zz_dephasing_gates(xn.sites, L, q)
            rho_after = rb_apply_channel_layer(xn.rho, zz, xn.trace_bra; maxdim=64, cutoff=1e-12)
            @test mps_fingerprint(rho_after) == mps_fingerprint(prod.rho)      # bitwise identical tensors
            @test rb_doubled_trace(xn.rho, xn.trace_bra) ≈ 1.0 atol = 1e-12
            @test size(xn.outcomes) == (3L, L) && xn.seed == seed
        end
    end

    @testset "q = 0: both observation times coincide; q > 0: they differ" begin
        L = 5
        for (q, same) in ((0.0, true), (0.3, false))
            kw = (lambda_x=0.7, q=q, T_max=2L, maxdim=64, cutoff=1e-12, seed=11)
            xn = rb_xn_evolve_right_edge_one_trial(L; kw...)
            prod = rb_evolve_right_edge_one_trial(L; kw...)
            ox = rb_renyi2_binder_one_trajectory(xn.rho, xn.sites, L; maxdim=256, cutoff=1e-14)
            op = rb_renyi2_binder_one_trajectory(prod.rho, prod.sites, L; maxdim=256, cutoff=1e-14)
            @test (isapprox(ox.B2, op.B2; rtol=1e-10) && isapprox(ox.purity, op.purity; rtol=1e-10)) == same
        end
    end

    @testset "production correlator vs dense reference (X-noise state)" begin
        for (L, q, seed) in ((3, 0.2, 21), (4, 0.425, 22), (5, 0.3, 23), (6, 0.5, 24))
            xn = rb_xn_evolve_right_edge_one_trial(L; lambda_x=0.7, q=q, T_max=2L, maxdim=64, cutoff=1e-12, seed=seed)
            ρ = dense_from_doubled_mps(xn.rho, L)
            @test real(tr(ρ)) ≈ 1.0 atol = 1e-10
            ref = dense_renyi2_moments(ρ, L)
            res = rb_renyi2_binder_one_trajectory(xn.rho, xn.sites, L; maxdim=256, cutoff=1e-14)
            @test res.purity ≈ ref.purity rtol = 1e-10
            @test res.M2 ≈ ref.M2 rtol = 1e-9
            @test res.M4 ≈ ref.M4 rtol = 1e-9
            @test res.B2 ≈ ref.B2 rtol = 1e-8
        end
    end

    @testset "runner: same seed derivation as rb_run_right_edge_point; paired trajectories" begin
        L, q, seed, n = 4, 0.3, 1234, 3
        kw = (lambda_x=0.7, lambda_zz=0.0, ntrials=n, T_max=2L, maxdim=64, cutoff=1e-12,
              obs_maxdim=256, obs_cutoff=1e-14, seed=seed, nboot=50)
        xn = rb_xn_run_right_edge_point(L, q; kw...)
        rng = MersenneTwister(seed)
        @test xn.trial_seeds == [Int(rand(rng, UInt32)) for _ in 1:n]
        @test xn.n_valid == n && xn.n_invalid == 0
        @test xn.B ≈ mean(xn.B2_per_trajectory)
        @test xn.B2_ratio_of_mean_moments ≈ 1 - mean(xn.M4_per_trajectory) / (3 * mean(xn.M2_per_trajectory)^2)

        # the production runner analyses the same trajectories, just observed after the last ZZ layer
        prod = rb_run_right_edge_point(L, q; kw...)
        after = map(xn.trial_seeds) do s
            ev = rb_evolve_right_edge_one_trial(L; lambda_x=0.7, q=q, T_max=2L, maxdim=64, cutoff=1e-12, seed=s)
            rb_renyi2_binder_one_trajectory(ev.rho, ev.sites, L; maxdim=256, cutoff=1e-14).B2
        end
        @test mean(after) ≈ prod.B rtol = 1e-8
        @test abs(mean(after) - xn.B) > 1e-3                      # the two observation times really differ at q > 0
    end

    @testset "no mutation / no global-RNG consumption" begin
        L = 6
        xn = rb_xn_evolve_right_edge_one_trial(L; lambda_x=0.7, q=0.3, T_max=2L, maxdim=64, cutoff=1e-12, seed=77)
        before = mps_fingerprint(xn.rho)
        Random.seed!(2024); expected_next = rand()
        Random.seed!(2024)
        o1 = rb_renyi2_binder_one_trajectory(xn.rho, xn.sites, L; maxdim=256, cutoff=1e-14)
        @test rand() == expected_next
        @test mps_fingerprint(xn.rho) == before
        @test rb_renyi2_binder_one_trajectory(xn.rho, xn.sites, L; maxdim=256, cutoff=1e-14) == o1
        xn2 = rb_xn_evolve_right_edge_one_trial(L; lambda_x=0.7, q=0.3, T_max=2L, maxdim=64, cutoff=1e-12, seed=77)
        @test mps_fingerprint(xn2.rho) == before
    end

    @testset "T_max = 1" begin
        xn = rb_xn_evolve_right_edge_one_trial(4; lambda_x=0.7, q=0.3, T_max=1, maxdim=64, cutoff=1e-12, seed=9)
        @test size(xn.outcomes) == (1, 4)
        @test rb_doubled_trace(xn.rho, xn.trace_bra) ≈ 1.0 atol = 1e-12
        @test_throws AssertionError rb_xn_evolve_right_edge_one_trial(4; lambda_x=0.7, q=0.3, T_max=0, maxdim=64, cutoff=1e-12, seed=9)
    end
end
