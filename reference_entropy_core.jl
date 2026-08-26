"""
Reference-entropy diagnostic core.

Reuses the EXISTING general-lambda doubled-MPS Born-sampling dynamics
from `lambda_scan_susceptibilities_itensorcorrelators.jl` (functions
prefixed `ls_`) UNCHANGED: weak-X measurement -> X-dephasing -> weak-ZZ
measurement -> ZZ-dephasing, with

    lambda_x  = delta * lambda,
    lambda_zz = delta * (1 - lambda),
    q_x = q_zz = q,

already documented there (delta = 0.7 for the right-boundary slice;
lambda = 1 gives lambda_x = 0.7, lambda_zz = 0, the primary q-scan
requested here). No Kraus operator, channel gate, sampling probability,
normalization rule, or layer order is modified in this file.

The ONLY new physics is the INITIAL STATE. Instead of the trivial
product state |0...0><0...0| on L system qubits, this file adds one
reference qubit R (NEVER measured, dephased, or gated) prepared in a
GHZ-Bell state with the logical |0...0>/|1...1> branches of the L
system qubits:

    |Psi_0>_{QR} = (|0...0>_Q |0>_R + |1...1>_Q |1>_R) / sqrt(2).

The doubled bra/ket MPS representation is extended by ONE extra
site-pair (bra_R at site 2L+1, ket_R at site 2L+2), appended AFTER the L
system pairs. Because every `ls_*` dynamics function takes `L` as an
explicit argument and only ever loops over positions 1..L (see
`ls_apply_weak_x_layer`, `ls_apply_weak_zz_layer`,
`ls_build_x_dephasing_gates`, `ls_build_zz_dephasing_gates`), calling
them with `L` = system size (not L+1) on the (L+1)-pair `sites` array
structurally guarantees that no channel or gate ever touches the
reference site -- checked explicitly in `reference_entropy_validate.jl`.

Vectorization convention (matching the rest of this repo): site 2p-1 =
bra_p (row/operator-left index), site 2p = ket_p (column/operator-right
index), for p = 1..L (system), p = L+1 (reference).

|Psi_0><Psi_0| = (1/2) * sum_{b,b' in {0,1}} |b>^{L+1}_bra <b'|^{L+1}_ket,
a sum of 4 UNENTANGLED-across-bra/ket product terms (each individually a
computational basis state across all 2(L+1) sites), built via ordinary
MPS addition -- no bond-dimension-4 hand construction is needed.

Reference entropy S_R(m) is computed per Born-sampled trajectory (record
m) by an EXACT partial trace over the L system qubits of the (already
trace-normalized) doubled-MPS rho, leaving the reference bra/ket legs
open, giving the reduced 2x2 matrix rho_R(m) directly (not via any MPS
bond-entropy shortcut, since rho_R can be mixed). S_R is then computed
per trajectory and averaged over trajectories (sample mean over
Born-sampled records) -- NEVER by first averaging density matrices and
then computing entropy of the average (see module docstring in
`reference_entropy_validate.jl`'s averaging test for the numerical
demonstration that these differ).
"""

using Random
using SHA
using LinearAlgebra
using ITensors
using ITensorMPS

include("lambda_scan_susceptibilities_itensorcorrelators.jl")

const REFERENCE_ENTROPY_EIGENVALUE_TOL = 1e-9

function reference_entropy_record_checksum(record::AbstractVector{<:Integer})
    return bytes2hex(sha256(join(record, ",")))
end

# ============================================================
# Explicit product-state MPS builder (bond dim 1; used only to build
# the 4 basis-state terms summed into the GHZ-Bell initial state).
# ============================================================
function reference_entropy_product_state(sites, bits::Vector{Int})
    tensors = ITensor[]
    for (n, s) in enumerate(sites)
        t = ITensor(s)
        t[s => bits[n] + 1] = 1.0
        push!(tensors, t)
    end
    return MPS(tensors)
end

"""
    reference_entropy_initial_state(sites, L)

|Psi_0><Psi_0| for the GHZ-Bell state entangling all L system qubits
with the reference qubit R (site pair L+1), vectorized as a doubled
bra/ket MPS over the given `sites` (length 2*(L+1)).
"""
function reference_entropy_initial_state(sites, L::Int)
    N = length(sites)
    @assert N == 2 * (L + 1)

    function pattern(bra_bit::Int, ket_bit::Int)
        bits = Vector{Int}(undef, N)
        for p in 1:(L + 1)
            bits[2p - 1] = bra_bit
            bits[2p] = ket_bit
        end
        return bits
    end

    term00 = reference_entropy_product_state(sites, pattern(0, 0))
    term01 = reference_entropy_product_state(sites, pattern(0, 1))
    term10 = reference_entropy_product_state(sites, pattern(1, 0))
    term11 = reference_entropy_product_state(sites, pattern(1, 1))

    return 0.5 * (term00 + term01 + term10 + term11)
end

# ============================================================
# Exact partial trace over the L system qubits, leaving the reference
# qubit's bra/ket legs open.
# ============================================================
function reference_entropy_bell_tensor(bra_idx, ket_idx)
    t = ITensor(bra_idx, ket_idx)
    t[bra_idx => 1, ket_idx => 1] = 1.0
    t[bra_idx => 2, ket_idx => 2] = 1.0
    return t
end

"""
    reference_entropy_reduced_rho_R(rho, sites, L)

Exact partial trace of the doubled-MPS `rho` over the L system qubits,
returning the 2x2 reduced reference density matrix as a plain Julia
matrix (rows/cols ordered by `sites[2L+1]` (bra_R), `sites[2L+2]`
(ket_R)). Since `rho` is already trace-normalized over ALL 2(L+1) sites
(Tr_{Q,R}[rho] = 1), Tr_R[rho_R] = Tr_{Q,R}[rho] = 1: no additional
normalization is needed here.
"""
function reference_entropy_reduced_rho_R(rho::MPS, sites, L::Int)
    env = reference_entropy_bell_tensor(sites[1], sites[2]) * rho[1] * rho[2]
    for p in 2:L
        bell = reference_entropy_bell_tensor(sites[2p - 1], sites[2p])
        env = env * (bell * rho[2p - 1] * rho[2p])
    end
    reduced = env * rho[2L + 1] * rho[2L + 2]
    bra_R, ket_R = sites[2L + 1], sites[2L + 2]
    return Array(reduced, bra_R, ket_R)
end

"""
    reference_entropy_von_neumann(M; eig_tol)

Hermitize `M`, verify Tr = 1, clip eigenvalues in [-eig_tol, 0) to zero
(fail loudly for anything more negative), and return the von Neumann
entropy in bits plus the diagnostic quantities required for the raw CSV
row (eigenvalues, trace/Hermiticity error, minimum eigenvalue).
"""
function reference_entropy_von_neumann(M::AbstractMatrix; eig_tol::Float64=REFERENCE_ENTROPY_EIGENVALUE_TOL)
    trace_val = real(tr(M))
    trace_error = abs(trace_val - 1.0)
    hermiticity_error = maximum(abs.(M .- M'))

    sym = Hermitian((M .+ M') ./ 2)
    raw_eigs = eigvals(sym)

    clipped = Float64[]
    for e in raw_eigs
        if e < 0
            @assert e >= -eig_tol "Reference density matrix eigenvalue too negative: $e (tolerance $eig_tol)"
            push!(clipped, 0.0)
        else
            push!(clipped, e)
        end
    end

    total = sum(clipped)
    @assert total > 0 "All reduced-reference eigenvalues are zero/negative."
    normalized = clipped ./ total

    S = 0.0
    for eta in normalized
        if eta > 0
            S -= eta * log2(eta)
        end
    end

    return (
        S=S, eig1=normalized[1], eig2=normalized[2],
        trace_error=trace_error, hermiticity_error=hermiticity_error,
        minimum_eigenvalue=minimum(clipped),
    )
end

"""
    reference_entropy_evolve_one_trial(L; lambda_x, lambda_zz, q_x, q_zz,
                                        T_max, maxdim, cutoff, seed, save_every)

Evolve one Born-sampled trajectory of the (Q,R) system using the
UNMODIFIED `ls_*` dynamics (weak-X -> X-dephasing -> weak-ZZ ->
ZZ-dephasing, called with `L` = system size so the reference site is
never touched), starting from the GHZ-Bell (Q,R) initial state.
Snapshots (time, reduced rho_R, entropy diagnostics, max bond
dimension) are recorded at t=0 and every `save_every` steps thereafter
(and always at t=T_max), preserving the FULL time series per trajectory
for later trajectory-level bootstrap resampling.
"""
function reference_entropy_evolve_one_trial(
    L::Int;
    lambda_x::Float64,
    lambda_zz::Float64,
    q_x::Float64,
    q_zz::Float64,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    seed::Int,
    save_every::Int=1,
)
    rng = MersenneTwister(seed)
    sites = siteinds("Qubit", 2 * (L + 1))
    trace_bra = ls_trace_state(sites)
    rho = ls_trace_normalize(reference_entropy_initial_state(sites, L), trace_bra)

    x_dephasing_gates = ls_build_x_dephasing_gates(sites, L, q_x)
    zz_dephasing_gates = ls_build_zz_dephasing_gates(sites, L, q_zz)

    record = Int[]
    max_trace_error = abs(ls_trace(rho, trace_bra) - 1.0)
    max_linkdim = maxlinkdim(rho)

    snapshots = NamedTuple[]
    function record_snapshot!(t::Int)
        M = reference_entropy_reduced_rho_R(rho, sites, L)
        entropy = reference_entropy_von_neumann(M)
        push!(snapshots, (time=t, rhoR=M, entropy=entropy))
    end

    record_snapshot!(0)

    for t in 1:T_max
        rho, x_outcomes = ls_apply_weak_x_layer(
            rho, sites, L, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )
        append!(record, x_outcomes)

        rho = ls_apply_channel_layer(rho, x_dephasing_gates, trace_bra; maxdim=maxdim, cutoff=cutoff)

        rho, zz_outcomes = ls_apply_weak_zz_layer(
            rho, sites, L, lambda_zz, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )
        append!(record, zz_outcomes)

        rho = ls_apply_channel_layer(rho, zz_dephasing_gates, trace_bra; maxdim=maxdim, cutoff=cutoff)

        max_trace_error = max(max_trace_error, abs(ls_trace(rho, trace_bra) - 1.0))
        max_linkdim = max(max_linkdim, maxlinkdim(rho))

        if t % save_every == 0 || t == T_max
            record_snapshot!(t)
        end
    end

    return (
        rho=rho, sites=sites,
        snapshots=snapshots,
        record=record,
        record_checksum=reference_entropy_record_checksum(record),
        max_trace_error=max_trace_error,
        max_interphysical_linkdim=max_linkdim,
    )
end
