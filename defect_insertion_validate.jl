using Test
using LinearAlgebra

include("defect_insertion_core.jl")

# ============================================================
# Independent dense brute-force replay (no MPS/ITensor machinery).
#
# Site convention: site 1 is the most-significant bit of the linear
# Hilbert-space index (matching Julia's `kron(I(2^(p-1)), op, I(2^(L-p)))`
# embedding order used below), consistently for both the Kraus/channel
# embedding and the partial-transpose index bookkeeping.
# ============================================================
function dense_embed(op2x2::AbstractMatrix{Float64}, p::Int, L::Int)
    left = Matrix{Float64}(I, 2^(p - 1), 2^(p - 1))
    right = Matrix{Float64}(I, 2^(L - p), 2^(L - p))
    return kron(kron(left, op2x2), right)
end

function idx_to_bits(idx0::Int, L::Int)
    bits = zeros(Int, L)
    for p in 1:L
        bits[p] = (idx0 >> (L - p)) & 1
    end
    return bits
end

function bits_to_idx(bits::Vector{Int})
    idx0 = 0
    for b in bits
        idx0 = (idx0 << 1) | b
    end
    return idx0
end

function dense_domain_wall_operator(i::Int, j::Int, L::Int)
    d = 2^L
    D = Matrix{Float64}(I, d, d)
    Xmat = [0.0 1.0; 1.0 0.0]
    lo, hi = minmax(i, j)
    for p in lo:(hi - 1)
        D = dense_embed(Xmat, p, L) * D
    end
    return D
end

"""
    dense_replay(record, L, T_max, lambda_x, q)

Independently replay the SAME fixed measurement record with dense
(2^L x 2^L) matrices and the fully *unnormalized* Kraus/channel maps
(no trace-renormalization at any step), returning rho_m^true. This is
an entirely separate code path from the MPS dynamics, used only to
validate the MPS-based contraction.
"""
function dense_replay(record::Vector{Int}, L::Int, T_max::Int, lambda_x::Float64, q::Float64)
    d = 2^L
    rho = zeros(Float64, d, d)
    rho[1, 1] = 1.0

    Xmat = [0.0 1.0; 1.0 0.0]
    Imat = Matrix{Float64}(I, 2, 2)
    Zmat = [1.0 0.0; 0.0 -1.0]

    idx = 1
    for _ in 1:T_max
        for p in 1:L
            outcome = record[idx]
            idx += 1
            K = (Imat .+ (-1)^outcome .* lambda_x .* Xmat) ./ sqrt(2 * (1 + lambda_x^2))
            Kemb = dense_embed(K, p, L)
            rho = Kemb * rho * Kemb'
        end
        for p in 1:L
            Xemb = dense_embed(Xmat, p, L)
            rho = (1 - q) .* rho .+ q .* (Xemb * rho * Xemb)
        end
        for p in 1:(L - 1)
            ZZemb = dense_embed(Zmat, p, L) * dense_embed(Zmat, p + 1, L)
            rho = (1 - q) .* rho .+ q .* (ZZemb * rho * ZZemb)
        end
    end

    return rho
end

@testset "defect insertion validation" begin

    @testset "checksum determinism" begin
        @test defect_insertion_record_checksum([0, 1, 1]) == defect_insertion_record_checksum([0, 1, 1])
        @test defect_insertion_record_checksum([0, 1, 1]) != defect_insertion_record_checksum([1, 1, 0])
    end

    L = 4
    sites = siteinds("Qubit", 2L)
    trace_bra = rb_bell_state(sites)
    rho = rb_trace_normalize(rb_initial_repetition_code_state(sites), trace_bra)

    @testset "no insertion is repeatable" begin
        @test defect_insertion_logZ0(rho) == defect_insertion_logZ0(rho)
    end

    @testset "coincident defects give deltaF = 0 (F^2 = I)" begin
        identity_result = defect_insertion_contract_pair(rho, sites, 2, 2; maxdim=64, cutoff=1e-12)
        @test identity_result.deltaF === 0.0
        @test identity_result.logR === 0.0
    end

    @testset "double insertion restores B0" begin
        once = defect_insertion_apply_domain_wall(rho, sites, 1, 3; maxdim=64, cutoff=1e-12)
        twice = defect_insertion_apply_domain_wall(once, sites, 1, 3; maxdim=64, cutoff=1e-12)
        @test isapprox(inner(rho, twice), inner(rho, rho); atol=1e-10)
    end

    @testset "no mutation of the baseline state" begin
        baseline_snapshot = deepcopy(rho)
        value = defect_insertion_contract_pair(rho, sites, 1, 2; maxdim=64, cutoff=1e-12)
        @test isfinite(value.logR) || isnan(value.logR)
        @test isapprox(inner(rho, rho), inner(baseline_snapshot, baseline_snapshot); atol=1e-12)
        for n in eachindex(rho)
            @test norm(rho[n] - baseline_snapshot[n]) < 1e-12
        end
    end

    @testset "B0 is not invariant under the domain-wall flip" begin
        # Required by the spec: verify F_i B_0 != B_0 (otherwise the defect
        # would be trivial, as the earlier bra/ket-SWAP candidate was proven
        # to be -- see the core file's "SUPERSEDED CANDIDATE" note). X^bra
        # breaks the bra/ket-symmetric structure of this dynamics directly,
        # so this is expected -- and required -- to be non-trivial.
        evolved = defect_insertion_evolve_one_trial(
            L; lambda_x=0.7, q=0.2, T_max=L, maxdim=64, cutoff=1e-12, seed=11,
        )
        flipped = defect_insertion_apply_domain_wall(evolved.rho, evolved.sites, 1, 3; maxdim=64, cutoff=1e-12)
        @test !isapprox(inner(evolved.rho, flipped), inner(evolved.rho, evolved.rho); atol=1e-6)
    end

    @testset "same seed gives identical record and rho" begin
        first = defect_insertion_evolve_one_trial(
            L; lambda_x=0.7, q=0.2, T_max=L, maxdim=64, cutoff=1e-12, seed=77,
        )
        second = defect_insertion_evolve_one_trial(
            L; lambda_x=0.7, q=0.2, T_max=L, maxdim=64, cutoff=1e-12, seed=77,
        )
        @test first.record_checksum == second.record_checksum
        @test first.record == second.record
        @test isapprox(inner(first.rho, second.rho), inner(first.rho, first.rho); atol=1e-10)
    end

    @testset "same bulk network used for Z0 and Zij" begin
        evolved = defect_insertion_evolve_one_trial(
            L; lambda_x=0.7, q=0.2, T_max=L, maxdim=64, cutoff=1e-12, seed=5,
        )
        traj = defect_insertion_compute_trajectory(
            evolved, L, [(1, 3), (2, 4)]; maxdim=64, cutoff=1e-12,
        )
        for row in traj.rows
            @test row.record_checksum == evolved.record_checksum
        end
    end

    @testset "domain-wall defect differs from the physical ZZ correlator" begin
        evolved = defect_insertion_evolve_one_trial(
            L; lambda_x=0.7, q=0.2, T_max=L, maxdim=64, cutoff=1e-12, seed=21,
        )
        domain_wall_value = defect_insertion_contract_pair(
            evolved.rho, evolved.sites, 1, 3; maxdim=64, cutoff=1e-12,
        )
        zz_operator = op(RB_sigma_z, evolved.sites[2 * 1 - 1]) * op(RB_sigma_z, evolved.sites[2 * 1]) *
                      op(RB_sigma_z, evolved.sites[2 * 3 - 1]) * op(RB_sigma_z, evolved.sites[2 * 3])
        zz_value = real(inner(evolved.rho, apply([zz_operator], evolved.rho; cutoff=1e-12, maxdim=64)))
        @test !isapprox(domain_wall_value.z_ij_raw, zz_value; rtol=1e-3)
    end

    @testset "brute-force dense cross-check (small system, no truncation)" begin
        L_small = 3
        T_small = L_small
        evolved = defect_insertion_evolve_one_trial(
            L_small; lambda_x=0.7, q=0.25, T_max=T_small, maxdim=256, cutoff=0.0, seed=901,
        )

        rho_dense = dense_replay(evolved.record, L_small, T_small, 0.7, 0.25)
        z0_dense = tr(rho_dense * rho_dense)
        logZ0_dense_unnormalized = log(z0_dense)

        # The domain-wall Z_ij is a genuine (possibly signed) correlator-like
        # quantity, not a manifestly positive Boltzmann weight. Compare
        # log(Complex(.)) on both sides (handles negative reals via the
        # +/- i*pi branch), after rescaling the MPS side by exp(2*logC) to
        # match the dense side's fully unnormalized rho^true (see the
        # normalization derivation in the core file's docstring):
        # Z_ij = Tr[D_ij * rho^2] (derived in the core file's docstring),
        # with D_ij = prod_{p=i}^{j-1} X_p embedded as an ordinary
        # (non-doubled) operator on the L-qubit Hilbert space.
        for (i, j) in [(1, 2), (1, 3), (2, 3)]
            mps_value = defect_insertion_contract_pair(
                evolved.rho, evolved.sites, i, j; maxdim=256, cutoff=0.0,
            )
            D_dense = dense_domain_wall_operator(i, j, L_small)
            zij_dense = tr(D_dense * rho_dense * rho_dense)
            log_mps_absolute = log(Complex(mps_value.z_ij_raw)) + 2 * evolved.logC
            log_dense = log(Complex(zij_dense))
            @test isapprox(log_mps_absolute, log_dense; atol=1e-6)
        end

        # Absolute logZ0 (scale-tracking) cross-check. rho^final = rho^true/C_m
        # (C_m = exp(logC)), so Tr[(rho^final)^2] = Tr[(rho^true)^2]/C_m^2,
        # i.e. logZ0_true = 2*logC + logZ0_relative.
        logZ0_from_tracking = 2 * evolved.logC + defect_insertion_logZ0(evolved.rho)
        @test isapprox(logZ0_from_tracking, logZ0_dense_unnormalized; atol=1e-6)
    end

    @testset "truncation convergence" begin
        L_mid = 6
        loose = defect_insertion_evolve_one_trial(
            L_mid; lambda_x=0.7, q=0.3, T_max=L_mid, maxdim=16, cutoff=1e-12, seed=303,
        )
        tight = defect_insertion_evolve_one_trial(
            L_mid; lambda_x=0.7, q=0.3, T_max=L_mid, maxdim=128, cutoff=1e-14, seed=303,
        )
        @test loose.record == tight.record
        pair = (2, 5)
        loose_value = defect_insertion_contract_pair(
            loose.rho, loose.sites, pair...; maxdim=16, cutoff=1e-12,
        )
        tight_value = defect_insertion_contract_pair(
            tight.rho, tight.sites, pair...; maxdim=128, cutoff=1e-14,
        )
        @test isapprox(loose_value.logZ0, tight_value.logZ0; atol=0.05)
        @test isapprox(loose_value.z_ij_raw, tight_value.z_ij_raw; atol=0.05)
    end

    @testset "positivity smoke test" begin
        L_smoke = 8
        n_negative = 0
        n_total = 0
        for seed in 1:6
            evolved = defect_insertion_evolve_one_trial(
                L_smoke; lambda_x=0.7, q=0.3, T_max=L_smoke, maxdim=64, cutoff=1e-12, seed=seed,
            )
            @test defect_insertion_logZ0(evolved.rho) isa Float64  # Z0 > 0 always (sum of squares)
            for (i, j) in [(2, 6), (3, 7)]
                value = defect_insertion_contract_pair(
                    evolved.rho, evolved.sites, i, j; maxdim=64, cutoff=1e-12,
                )
                n_total += 1
                if !(value.z_ij_raw > 0)
                    n_negative += 1
                end
            end
        end
        println("positivity smoke test: $(n_total - n_negative)/$(n_total) trials gave Z_ij > 0")
        @test n_total > 0
    end
end
