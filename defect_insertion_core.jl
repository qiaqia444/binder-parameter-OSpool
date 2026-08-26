"""
Right-boundary defect-pair partition-function diagnostic.

This file reuses the *existing*, unmodified right-boundary dynamics from
`renyi2_right_boundary_susceptibility_core.jl` (doubled bra/ket MPS
representation of a density matrix rho_m; site 2p-1 = bra_p, site 2p =
ket_p for physical position p=1..L) to generate one measurement record m
and its final state rho_m at T_max = L. No bulk tensor, Kraus operator,
or channel gate is redefined here; only per-step *bookkeeping* (recording
outcomes and the log-scale removed by trace-normalization) is added on
top of the existing per-site/per-layer primitives.

======================================================================
Derivation: what Z_m^(0) and Z_m^(ij) actually are in this formalism
======================================================================

n=2 replica partition function (no defect)
-------------------------------------------
The doubled bra/ket MPS already vectorizes rho_m: site 2p-1 carries the
"bra" (row) index and site 2p carries the "ket" (column) index of
rho_m's matrix element rho_m(a,b) at physical position p. The ordinary
ITensor inner product of rho_m with itself is

    <<rho_m|rho_m>> = sum_{a,b} conj(rho_m(a,b)) rho_m(a,b)
                    = Tr[rho_m^dagger rho_m] = Tr[rho_m^2]

using Hermiticity of a physical density matrix. Tr[rho_m^2] *is* the n=2
Renyi partition function with the ordinary ("untwisted") boundary
condition:

    Z_m^(0) := Tr[rho_m^2] = <<rho_m|rho_m>>.

This equals `rb_hilbert_schmidt_norm(rho_m)` -- but here that equality is
derived explicitly (not assumed) and is checked against an independent
dense brute-force calculation in `defect_insertion_validate.jl`
(validation items 6-7).

Right-boundary domain-wall defect (replica-overlap sign flip)
----------------------------------------------------------------
The Renyi-2 replica-overlap variable at site p is

    mu_p = Z_p^bra Z_p^ket    (the same operator whose sum Q = sum_p mu_p
                               already defines the existing Binder/
                               susceptibility observables in this repo).

Inserting a defect pair (i,j) means creating a domain wall in mu_p: mu_p
is REVERSED (mu_p -> -mu_p) at every site along the segment [i, j-1],
with the two endpoints i and j left as the two domain-wall kinks (mu_i,
mu_j themselves are not separately pinned; they are simply where the
sign change begins/ends). A single-site operator that reverses mu_p
without touching any other site is the physical X operator acting on
the BRA leg alone:

    X_p^bra (Z_p^bra Z_p^ket) X_p^bra = -Z_p^bra Z_p^ket,

since X anticommutes with Z on the bra leg and commutes with Z on the
untouched ket leg. The domain-wall (disorder-string) operator for the
pair (i,j) is therefore

    D_{i,j} = prod_{p=i}^{j-1} X_p^bra,          D_{i,j}^2 = Identity,

(product over the HALF-OPEN segment [i, j-1], matching the requested
separation r = |i-j| = L/2 and reducing to the identity when i = j, an
empty product). This is the standard Kramers-Wannier-type disorder
string construction: a local order-parameter mu_p flip supported on a
finite segment, whose two endpoints are the only places the operator's
support begins/ends. Z_m^{(ij)} is this string sandwiched around the
SAME final state used for Z_m^{(0)}:

    Z_m^{(0)}   := <<rho_m | rho_m>>          = Tr[rho_m^2],
    Z_m^{(ij)}  := <<rho_m | D_{i,j} | rho_m>>.

This is the local, involutive, *boundary-only* modification requested:
- it acts only on the already-evolved final bra legs at i,...,j-1;
- it never touches a bulk tensor, Kraus operator, or channel gate;
- D_{i,j}^2 = I (a plain product of commuting single-site involutions),
  so i = j gives back the ordinary boundary condition exactly
  (Z_m^{(ii)} = Z_m^{(0)}, Delta F = 0);
- B_0 is *not* invariant under D_{i,j} for i != j: X^bra breaks the
  bra/ket-symmetric structure of this dynamics directly (unlike a
  bra/ket SWAP, which was found to be an exact symmetry of every state
  in this model -- see "PROVEN LIMITATION" note retained below for the
  earlier, ruled-out candidate and why it failed).

Why this is *not* the physical correlator
    <rho | Z_(2i-1) Z_(2i) Z_(2j-1) Z_(2j) | rho>
------------------------------------------------
Z_i^bra Z_i^ket is *diagonal* in the {|00>,|01>,|10>,|11>} bra/ket basis
of site i (eigenvalues +1,-1,-1,+1): it only reweights branches, never
mixes them, and it commutes with itself trivially (it is its own
square-root of identity in a trivial sense: applying it twice is a
no-op that changes nothing at all, since it is diagonal and squares to
the identity on EVERY branch simultaneously). X_p^bra is *off-diagonal*
(it exchanges |0*> <-> |1*> on the bra leg only) and, crucially, acts on
ONLY the bra leg -- not both legs symmetrically -- which is exactly what
lets it reverse the SIGN of mu_p instead of leaving it invariant. This
is numerically checked directly in `defect_insertion_validate.jl`.

Normalization / preserved scale factors
------------------------------------------
`renyi2_right_boundary_susceptibility_core.jl` trace-renormalizes rho
after every weak-measurement site and every channel layer. Because every
step of the dynamics (Kraus maps, dephasing channels) is *linear*, doing
this is exactly equivalent to computing the fully unnormalized
rho_m^true once and dividing by a single overall positive scalar C_m
(the accumulated product of Born/trace weights): rho_m^final = rho_m^true / C_m.
Since

    logZ_m^{(ij)} - logZ_m^{(0)}
      = log(Tr[rho_m^final (rho_m^final)^{Tij}]) - log(Tr[(rho_m^final)^2])
      = log(Tr[rho_m^true (rho_m^true)^{Tij}] / C_m^2)
        - log(Tr[(rho_m^true)^2] / C_m^2)

the common factor C_m^2 cancels *exactly* in logR = logZ_ij - logZ_0
(the requested primary observable), independent of how/when the
dynamics renormalizes. This file nonetheless *also* tracks the
accumulated log-scale `logC` explicitly (Born probability of every
sampled weak-X outcome, plus the trace removed by every channel-layer
renormalization) so that the reported `logZ0`/`logZij` are absolute, and
so that the cancellation claim itself is checked numerically
(`"scale_invariance"` validation test): logC + log(Tr[(rho^final)^2])
must equal the fully unnormalized log(Tr[(rho^true)^2]) from an
independent dense replay of the same fixed record.

======================================================================
SUPERSEDED CANDIDATE (kept for the record): a bra/ket SWAP at just the
two isolated points i,j is an exact symmetry of every rho_m here
======================================================================
An earlier version of this file used F_p = SWAP(bra_p, ket_p) applied
only at the two isolated points i and j (not the segment between them),
i.e. treating the defect as a partial transpose of rho_m at two points.
That construction was numerically found -- and then proven analytically
-- to be an EXACT symmetry of every rho_m produced by this dynamics, for
every i, j, L, q, and seed tested (Z_m^{(ij)} = Z_m^{(0)} to machine
precision, always), because every operator used anywhere in this
dynamics (the weak-X Kraus operator, the X-dephasing mixture components
{I, X}, and the ZZ-dephasing mixture components {I, Z tensor Z}) is real
and symmetric (O = O^T), and a real-symmetric-only dynamics starting
from a real product state preserves invariance under partial transpose
at every site, for every step, as an exact algebraic identity. Two
*isolated point* flips can therefore never create a nontrivial domain
wall in this representation.

This is exactly the contingency the request itself anticipated ("verify
that B_0 is not invariant under the proposed local flip... in that
case, identify the correct fixed/twisted boundary tensor instead").
Following a subsequent clarification that the defect is a domain-wall
*segment* flip (not two isolated point flips), this file now uses
`D_{i,j} = prod_{p=i}^{j-1} X_p^bra` instead (see above), which is *not*
a symmetry of rho_m (validated in `defect_insertion_validate.jl`).
"""

using Random
using SHA
using LinearAlgebra
using ITensors
using ITensorMPS

include("renyi2_right_boundary_susceptibility_core.jl")

"""
    DISORDER_STRING_PAIRS[L][r]

For each system size L, a Dict mapping separation r to two translated
endpoint pairs (i,j) with |i-j| = r, used to measure the r-dependence
(and, for r = L/2, comparability with the earlier single-length pairs).
The two pairs at each r are translated copies within the SAME
trajectory (not independent disorder samples): average them within a
trajectory before computing error bars across trajectories.
"""
const DISORDER_STRING_PAIRS = Dict(
    16 => Dict(
        2 => [(4, 6),  (5, 7)],
        4 => [(4, 8),  (5, 9)],
        6 => [(4, 10), (5, 11)],
        8 => [(4, 12), (5, 13)],
    ),
    24 => Dict(
        3  => [(6, 9),  (7, 10)],
        6  => [(6, 12), (7, 13)],
        9  => [(6, 15), (7, 16)],
        12 => [(6, 18), (7, 19)],
    ),
    32 => Dict(
        4  => [(8, 12), (9, 13)],
        8  => [(8, 16), (9, 17)],
        12 => [(8, 20), (9, 21)],
        16 => [(8, 24), (9, 25)],
    ),
)

function defect_insertion_record_checksum(record::AbstractVector{<:Integer})
    return bytes2hex(sha256(join(record, ",")))
end

# ============================================================
# Instrumented layers: same physics as
# `rb_apply_weak_x_layer` / `rb_apply_channel_layer`, reusing the exact
# same per-site/per-gate primitives, but additionally returning the
# log-scale removed by trace-renormalization at each step.
# ============================================================
function defect_insertion_weak_x_layer_tracked(
    rho::MPS,
    sites,
    L::Int,
    lambda_x::Float64,
    trace_bra::MPS,
    rng::AbstractRNG;
    maxdim::Int,
    cutoff::Float64,
)
    outcomes = Vector{Int}(undef, L)
    logscale = 0.0

    for i in 1:L
        rho, outcome, probabilities = rb_apply_weak_x_measurement_site(
            rho, sites, i, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )
        outcomes[i] = outcome
        logscale += log(probabilities[outcome + 1])
    end

    return rho, outcomes, logscale
end

function defect_insertion_channel_layer_tracked(
    rho::MPS, gates::Vector{ITensor}, trace_bra::MPS; maxdim::Int, cutoff::Float64,
)
    for gate in gates
        rho = apply([gate], rho; cutoff=cutoff, maxdim=maxdim)
    end

    trace = rb_doubled_trace(rho, trace_bra)
    @assert isfinite(trace) && trace > 0 "Non-positive channel-layer trace: $trace"

    return rho / trace, log(trace)
end

# ============================================================
# One Born-sampled trajectory: record m and final rho_m at T_max = L.
# ============================================================
function defect_insertion_evolve_one_trial(
    L::Int;
    lambda_x::Float64,
    q::Float64,
    T_max::Int,
    maxdim::Int,
    cutoff::Float64,
    seed::Int,
)
    @assert T_max == L "The defect diagnostic requires T_max = L."

    rng = MersenneTwister(seed)
    sites = siteinds("Qubit", 2L)
    trace_bra = rb_bell_state(sites)
    rho = rb_trace_normalize(rb_initial_repetition_code_state(sites), trace_bra)

    x_gates = rb_build_x_dephasing_gates(sites, L, q)
    zz_gates = rb_build_zz_dephasing_gates(sites, L, q)

    record = Int[]
    logC = 0.0
    max_trace_error = abs(rb_doubled_trace(rho, trace_bra) - 1.0)
    max_linkdim = rb_max_interphysical_linkdim(rho, L)

    for _ in 1:T_max
        rho, outcomes, logscale_x = defect_insertion_weak_x_layer_tracked(
            rho, sites, L, lambda_x, trace_bra, rng; maxdim=maxdim, cutoff=cutoff,
        )
        append!(record, outcomes)
        logC += logscale_x

        rho, logscale_xdeph = defect_insertion_channel_layer_tracked(
            rho, x_gates, trace_bra; maxdim=maxdim, cutoff=cutoff,
        )
        logC += logscale_xdeph

        rho, logscale_zzdeph = defect_insertion_channel_layer_tracked(
            rho, zz_gates, trace_bra; maxdim=maxdim, cutoff=cutoff,
        )
        logC += logscale_zzdeph

        max_trace_error = max(
            max_trace_error, abs(rb_doubled_trace(rho, trace_bra) - 1.0),
        )
        max_linkdim = max(max_linkdim, rb_max_interphysical_linkdim(rho, L))
    end

    return (
        rho=rho,
        sites=sites,
        trace_bra=trace_bra,
        record=record,
        record_checksum=defect_insertion_record_checksum(record),
        logC=logC,
        max_trace_error=max_trace_error,
        max_interphysical_linkdim=max_linkdim,
    )
end

# ============================================================
# Right-boundary domain-wall defect: D_{i,j} = prod_{p=i}^{j-1} X_p^bra.
# ============================================================
function defect_insertion_domain_wall_segment(i::Int, j::Int)
    lo, hi = minmax(i, j)
    return lo:(hi - 1)
 end

function defect_insertion_apply_domain_wall(
    rho::MPS, sites, i::Int, j::Int; maxdim::Int, cutoff::Float64,
)
    for p in defect_insertion_domain_wall_segment(i, j)
        rho = apply([op(RB_sigma_x, sites[2p - 1])], rho; cutoff=cutoff, maxdim=maxdim)
    end
    return rho
end

# ============================================================
# R_{ij} = zij/z0 for the same, unmodified rho_m -- a signed Renyi-2
# disorder-string overlap, not a partition-function ratio (see module
# docstring). No trajectory is ever discarded and R is never clipped or
# replaced by |R|.
# ============================================================

"""
Absolute threshold on |R| below which R is classified as numerically
zero (`is_numerical_zero = true`, `string_log_magnitude = Inf`). This is
a fixed, documented convention, not an arbitrary positive floor: R is
never altered, only classified.
"""
const DEFECT_INSERTION_ZERO_ATOL = 1e-12

"""
    defect_insertion_contract_pair(rho, sites, i, j; maxdim, cutoff)

Compute z0 = Tr[rho^2], zij = <<rho | D_ij | rho>> (via the disorder
string D_{i,j} = prod_{p=i}^{j-1} X_p^bra), and the SIGNED overlap
R = zij/z0, for the SAME `rho` used for both contractions. `rho` is
never mutated (`apply` returns a new MPS). If i == j, D_ii = Identity
exactly (empty segment), so R_ii = 1 exactly, with no contraction
performed. R is never discarded, clipped, or replaced by |R|; a
numerically-zero R is flagged via `is_numerical_zero` and reported with
`string_log_magnitude = Inf` rather than NaN.
"""
function defect_insertion_contract_pair(
    rho::MPS, sites, i::Int, j::Int; maxdim::Int, cutoff::Float64,
)
    z0_complex = inner(rho, rho)
    contraction_error_z0 = abs(imag(z0_complex)) / max(abs(z0_complex), eps())
    z0 = real(z0_complex)
    @assert isfinite(z0) && z0 > 0 "Non-positive/non-finite z0 = Tr[rho^2]: $z0"

    if i == j
        return (
            i=i, j=j, r=0, z0=z0, zij=z0, R=1.0, absR=1.0, signR=1.0,
            string_log_magnitude=0.0, is_numerical_zero=false,
            contraction_error=contraction_error_z0,
        )
    end

    flipped = defect_insertion_apply_domain_wall(rho, sites, i, j; maxdim=maxdim, cutoff=cutoff)
    zij_complex = inner(rho, flipped)
    contraction_error = max(
        contraction_error_z0,
        abs(imag(zij_complex)) / max(abs(zij_complex), eps()),
    )
    zij = real(zij_complex)

    R = zij / z0
    absR = abs(R)
    is_zero = absR <= DEFECT_INSERTION_ZERO_ATOL
    signR = is_zero ? 0.0 : sign(R)
    string_log_magnitude = is_zero ? Inf : -log(absR)

    return (
        i=i, j=j, r=abs(i - j), z0=z0, zij=zij, R=R, absR=absR, signR=signR,
        string_log_magnitude=string_log_magnitude, is_numerical_zero=is_zero,
        contraction_error=contraction_error,
    )
end
