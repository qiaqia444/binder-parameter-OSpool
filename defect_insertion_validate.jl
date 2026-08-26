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

    @testset "1. D_ii = I so R_ii = 1 exactly" begin
        identity_result = defect_insertion_contract_pair(rho, sites, 2, 2; maxdim=64, cutoff=1e-12)
        @test identity_result.R === 1.0
        @test identity_result.absR === 1.0
        @test identity_result.signR === 1.0
        @test identity_result.string_log_magnitude === 0.0
        @test identity_result.is_numerical_zero === false
        @test identity_result.r == 0
    end

    @testset "2. D_ij^2 = I" begin
        once = defect_insertion_apply_domain_wall(rho, sites, 1, 3; maxdim=64, cutoff=1e-12)
        twice = defect_insertion_apply_domain_wall(once, sites, 1, 3; maxdim=64, cutoff=1e-12)
        @test isapprox(inner(rho, twice), inner(rho, rho); atol=1e-10)
    end

    @testset "3. applied segment is exactly p = i,...,j-1" begin
        @test collect(defect_insertion_domain_wall_segment(3, 7)) == collect(3:6)
        @test collect(defect_insertion_domain_wall_segment(7, 3)) == collect(3:6)  # order-independent
        @test collect(defect_insertion_domain_wall_segment(5, 5)) == Int[]         # empty for i == j
        @test collect(defect_insertion_domain_wall_segment(4, 5)) == [4]           # r = 1: single site
    end

    @testset "5. |R_ij| <= 1 (Cauchy-Schwarz, D_ij unitary+Hermitian)" begin
        for seed in 1:5
            evolved = defect_insertion_evolve_one_trial(
                L; lambda_x=0.7, q=0.3, T_max=L, maxdim=64, cutoff=1e-12, seed=seed,
            )
            for (i, j) in [(1, 2), (1, 3), (1, 4), (2, 4)]
                value = defect_insertion_contract_pair(
                    evolved.rho, evolved.sites, i, j; maxdim=64, cutoff=1e-12,
                )
                @test value.absR <= 1.0 + 1e-8
            end
        end
    end

    @testset "6. no mutation of the baseline state" begin
        baseline_snapshot = deepcopy(rho)
        value = defect_insertion_contract_pair(rho, sites, 1, 2; maxdim=64, cutoff=1e-12)
        @test isfinite(value.R)
        @test isapprox(inner(rho, rho), inner(baseline_snapshot, baseline_snapshot); atol=1e-12)
        for n in eachindex(rho)
            @test norm(rho[n] - baseline_snapshot[n]) < 1e-12
        end
    end

    @testset "B0 is not invariant under the domain-wall flip" begin
        # X^bra breaks the bra/ket-symmetric structure of this dynamics
        # directly, unlike the earlier bra/ket-SWAP candidate (see core
        # file's "SUPERSEDED CANDIDATE" note), so this is expected to be
        # non-trivial.
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

    @testset "same bulk network used for z0 and zij" begin
        evolved = defect_insertion_evolve_one_trial(
            L; lambda_x=0.7, q=0.2, T_max=L, maxdim=64, cutoff=1e-12, seed=5,
        )
        for (i, j) in [(1, 3), (2, 4)]
            value = defect_insertion_contract_pair(
                evolved.rho, evolved.sites, i, j; maxdim=64, cutoff=1e-12,
            )
            @test isfinite(value.z0) && value.z0 > 0
        end
    end

    @testset "7. imaginary-to-total contraction ratio below tolerance" begin
        evolved = defect_insertion_evolve_one_trial(
            L; lambda_x=0.7, q=0.2, T_max=L, maxdim=64, cutoff=1e-12, seed=21,
        )
        for (i, j) in [(1, 3), (2, 4)]
            value = defect_insertion_contract_pair(
                evolved.rho, evolved.sites, i, j; maxdim=64, cutoff=1e-12,
            )
            @test value.contraction_error < 1e-8
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
        @test !isapprox(domain_wall_value.zij, zz_value; rtol=1e-3)
    end

    @testset "4. brute-force dense cross-check (small system, no truncation)" begin
        L_small = 3
        T_small = L_small
        evolved = defect_insertion_evolve_one_trial(
            L_small; lambda_x=0.7, q=0.25, T_max=T_small, maxdim=256, cutoff=0.0, seed=901,
        )

        rho_dense = dense_replay(evolved.record, L_small, T_small, 0.7, 0.25)
        z0_dense = tr(rho_dense * rho_dense)

        # zij = Tr[D_ij * rho^2] (derived in the core file's docstring), with
        # D_ij = prod_{p=i}^{j-1} X_p embedded as an ordinary (non-doubled)
        # operator on the L-qubit Hilbert space. Compare the SIGNED,
        # absolute-scale zij/z0 directly (accounting for the MPS-side
        # trace-renormalization factor C_m via evolved.logC).
        for (i, j) in [(1, 2), (1, 3), (2, 3)]
            mps_value = defect_insertion_contract_pair(
                evolved.rho, evolved.sites, i, j; maxdim=256, cutoff=0.0,
            )
            D_dense = dense_domain_wall_operator(i, j, L_small)
            zij_dense = tr(D_dense * rho_dense * rho_dense)

            C2 = exp(2 * evolved.logC)
            @test isapprox(mps_value.z0 * C2, z0_dense; atol=1e-6, rtol=1e-6)
            @test isapprox(mps_value.zij * C2, zij_dense; atol=1e-6, rtol=1e-6)

            # R is scale-invariant: matches directly without any C_m correction.
            R_dense = zij_dense / z0_dense
            @test isapprox(mps_value.R, R_dense; atol=1e-8)
        end
    end

    @testset "11. truncation / cutoff convergence" begin
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
        @test isapprox(loose_value.z0, tight_value.z0; atol=0.05)
        @test isapprox(loose_value.R, tight_value.R; atol=0.05)
    end

    @testset "8./9. bookkeeping: counts sum to total, nothing silently discarded" begin
        L_smoke = 8
        n_total = 0
        n_positive = 0
        n_negative = 0
        n_zero = 0
        collected = NamedTuple[]

        for seed in 1:8
            evolved = defect_insertion_evolve_one_trial(
                L_smoke; lambda_x=0.7, q=0.3, T_max=L_smoke, maxdim=64, cutoff=1e-12, seed=seed,
            )
            for (i, j) in [(2, 6), (3, 7)]
                value = defect_insertion_contract_pair(
                    evolved.rho, evolved.sites, i, j; maxdim=64, cutoff=1e-12,
                )
                push!(collected, value)  # every value retained -- R < 0 included
                n_total += 1
                if value.is_numerical_zero
                    n_zero += 1
                elseif value.R > 0
                    n_positive += 1
                else
                    n_negative += 1
                end
            end
        end

        @test n_positive + n_negative + n_zero == n_total
        @test length(collected) == n_total  # nothing filtered out of `collected`
        @test any(v -> v.R < 0, collected) || n_negative == 0  # negative R, if any, is present in `collected`
        println(
            "positivity bookkeeping: n_total=$n_total n_positive=$n_positive ",
            "n_negative=$n_negative n_zero=$n_zero",
        )
    end

    @testset "10. translation-equivalent pairs both compute" begin
        L_test = 16
        evolved = defect_insertion_evolve_one_trial(
            L_test; lambda_x=0.7, q=0.0, T_max=L_test, maxdim=64, cutoff=1e-12, seed=1,
        )
        for r in sort(collect(keys(DISORDER_STRING_PAIRS[L_test])))
            values = [
                defect_insertion_contract_pair(
                    evolved.rho, evolved.sites, i, j; maxdim=64, cutoff=1e-12,
                )
                for (i, j) in DISORDER_STRING_PAIRS[L_test][r]
            ]
            @test all(v -> isfinite(v.R), values)
            @test length(values) == 2
        end
    end

    @testset "q = 0.5 investigation: does R vanish within numerical accuracy?" begin
        # The full-dephasing boundary q = 0.5 previously produced an
        # anomalously large string_log_magnitude (~50). Investigate whether
        # this is a genuine effect or numerical noise around R ~ 0, by
        # comparing MPS (loose vs. tight truncation) against an independent
        # dense brute-force calculation at small L.
        L_q5 = 4
        loose = defect_insertion_evolve_one_trial(
            L_q5; lambda_x=0.7, q=0.5, T_max=L_q5, maxdim=32, cutoff=1e-10, seed=42,
        )
        tight = defect_insertion_evolve_one_trial(
            L_q5; lambda_x=0.7, q=0.5, T_max=L_q5, maxdim=256, cutoff=0.0, seed=42,
        )
        @test loose.record == tight.record

        rho_dense = dense_replay(tight.record, L_q5, L_q5, 0.7, 0.5)
        z0_dense = tr(rho_dense * rho_dense)

        for (i, j) in [(1, 2), (1, 3), (2, 4)]
            loose_value = defect_insertion_contract_pair(
                loose.rho, loose.sites, i, j; maxdim=32, cutoff=1e-10,
            )
            tight_value = defect_insertion_contract_pair(
                tight.rho, tight.sites, i, j; maxdim=256, cutoff=0.0,
            )
            D_dense = dense_domain_wall_operator(i, j, L_q5)
            zij_dense = tr(D_dense * rho_dense * rho_dense)
            R_dense = zij_dense / z0_dense

            println(
                "q=0.5 (i=$i,j=$j): loose R=$(loose_value.R) tight R=$(tight_value.R) ",
                "dense R=$R_dense  (loose zij=$(loose_value.zij), tight zij=$(tight_value.zij), ",
                "dense zij=$zij_dense, dense z0=$z0_dense)",
            )
            if abs(R_dense) < 1e-6
                println("  -> string overlap vanishes within numerical accuracy at (i=$i,j=$j)")
            end
            @test isapprox(tight_value.R, R_dense; atol=1e-6)
        end
    end
end
