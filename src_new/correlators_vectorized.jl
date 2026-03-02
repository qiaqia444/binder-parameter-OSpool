"""
ITensorCorrelators - Compute expectation values for vectorized density matrices.

This module computes <1|O|ρ> where |ρ> is a vectorized density matrix MPS
(each site has dim=d², encoding a d×d matrix).
"""

module ITensorCorrelators

using ITensors, ITensorMPS

export correlator

# Manual parity function (counts inversions in a permutation)
function parity(perm)
    n = length(perm)
    inversions = 0
    for i in 1:(n-1)
        for j in (i+1):n
            if perm[i] > perm[j]
                inversions += 1
            end
        end
    end
    return inversions % 2
end

# ------------------------------------------------------------
# Build local bra tensor for <1| : contracts a fused d^2 site index
# with vec(I) (i=j entries = 1, others = 0)
# ------------------------------------------------------------
function onebra_site_tensor(s::Index)
  dim_s = dim(s)
  d = Int(round(sqrt(dim_s)))
  if d * d != dim_s
    error("Site index dim = $dim_s is not a perfect square; expected d^2 for superket |rho>.")
  end
  b = ITensor(dag(s))
  # vec(I): α = i + (i-1)*d (for i=j)
  for i in 1:d
    α = i + (i - 1) * d
    b[dag(s) => α] = 1.0
  end
  return b
end

function make_onebra(psi::MPS)
  s = siteinds(psi)
  return [onebra_site_tensor(s[n]) for n in 1:length(psi)]
end

# ------------------------------------------------------------
# Public API
# ------------------------------------------------------------
function correlator(
  psi,
  cor_ops,
  op_sites,
)
  indices = similar(op_sites)
  op_sorted = similar(op_sites)

  for (i, op_site) in enumerate(op_sites)
    op_sorted[i] = tuple(sort([op_site...])...)
    indices[i] = tuple(sortperm([op_site...])...)
  end

  i = 1
  ind_sites = unique(indices)[i]
  op_sites_ord = op_sites
  cor_ops_ord = cor_ops

  C = correlator_recursive_compact(psi, cor_ops_ord, op_sites_ord; indices=ind_sites)
  return C
end

# ------------------------------------------------------------
# Recursive correlator: computes <1| O |psi> (psi = |rho>)
# ------------------------------------------------------------
function correlator_recursive_compact(
  psi,
  ops,
  sites;
  indices=nothing,
)
  if indices === nothing
    indices = collect(1:length(ops))
  end

  # Special case: no operators, just compute Tr(ρ) = <1|ρ>
  if length(ops) == 0 || length(sites) == 0 || (length(sites) > 0 && length(sites[1]) == 0)
    onebra = make_onebra(psi)
    L = ITensor(1.0)
    for n in 1:length(psi)
      L *= psi[n] * onebra[n]
    end
    C = Dict{Tuple{Vararg{Int64}},ComplexF64}()
    C[()] = scalar(L)
    return C
  end

  sites = sort(sites)
  N = length(sites[1])

  C = Dict{Tuple{Vararg{Int64}},ComplexF64}()

  # Pre-build <1| tensors (product bra)
  onebra = make_onebra(psi)

  # First index, repeats, perms
  inds_ord = [sort([sites[idx]...])[1] for idx in 1:length(sites)]
  repeats  = [count(==(sort([sites[idx]...])[1]), sites[idx]) - 1 for idx in 1:length(sites)]
  perms    = [sortperm([sites[idx]...])[1:(repeats[idx] + 1)] for idx in 1:length(sites)]

  op_inds = unique([(inds_ord[idx], repeats[idx], perms[idx]) for idx in 1:length(sites)])
  op_inds = sort(op_inds, by = x -> x[1])

  s = siteinds(psi)

  L = ITensor(1.0)
  counter = 1
  jw = 0
  element = zeros(Int64, N)

  add_operator_fermi(
    op_inds, sites, L, counter, element, N, ops, s, psi, onebra, C, indices, jw
  )
  return C
end

# ------------------------------------------------------------
# Fermionic/JW version
# ------------------------------------------------------------
function add_operator_fermi(
  op_inds,
  sites_ind_prev,
  L_prev,
  counter,
  element,
  N,
  ops,
  s,
  psi,
  onebra,
  C,
  indices,
  jw,
)
  for (a, op_ind) in enumerate(op_inds)
    repeat   = op_ind[2]
    perm_ind = op_ind[3]
    site     = op_ind[1]

    element[counter:(counter + repeat)] .= site

    # Initialize L_prev = <1| psi(1..site-1) with JW if needed
    if counter == 1
      L_prev = ITensor(1.0)
      if site > 1
        for str in 1:(site - 1)
          if jw % 2 != 0
            L_prev = L_prev * apply(op("F", s[str]), psi[str]) * onebra[str]
          else
            L_prev = L_prev * psi[str] * onebra[str]
          end
        end
      end
    end

    L = L_prev
    op_psi = psi[site]

    # if JW odd, apply F before site ops
    if jw % 2 != 0
      op_psi = apply(op("F", s[site]), op_psi)
    end

    # apply ops on this site (handle repeats)
    for i in 0:repeat
      op_psi = apply(op(ops[perm_ind[counter + repeat - i]], s[site]), op_psi)
    end

    # update JW count
    jw_next = jw
    for i in 0:repeat
      if has_fermion_string(ops[perm_ind[counter + repeat - i]], s[site])
        jw_next += 1
        op_psi = apply(op("F", s[site]), op_psi)
      end
    end

    # contract this site with <1|
    L = L * op_psi * onebra[site]

    # base case: last operator inserted
    if counter + repeat == N
      # contract tail to the right
      if site < length(psi)
        for str in (site + 1):length(psi)
          if jw_next % 2 != 0
            L = L * apply(op("F", s[str]), psi[str]) * onebra[str]
          else
            L = L * psi[str] * onebra[str]
          end
        end
      end

      # reorder elements by perm_ind
      perm_elem = element[sortperm(perm_ind)]

      # fermionic sign from reordering
      ferm_sites = Int.(perm_elem[findall(x -> has_fermion_string(x, s[site]), ops)])
      par = 1 - 2 * parity(sortperm(ferm_sites))

      C[tuple(perm_elem...)] = par * scalar(L)

    else
      # consistent site tuples
      sites_ind = sites_ind_prev[findall(
        x -> sort([x...])[counter + repeat] == site, sites_ind_prev
      )]

      # remove immediate repeats already handled
      deleteat!(
        sites_ind, findall(x -> sort([x...])[counter + repeat + 1] == site, sites_ind)
      )

      # next repeats/perms
      repeat_next = [
        count(==(sort([sites_ind[idx]...])[counter + repeat + 1]), sites_ind[idx]) - 1
        for idx in 1:length(sites_ind)
      ]

      inds_ord = [
        sort([sites_ind[idx]...])[counter + repeat + 1] for idx in 1:length(sites_ind)
      ]

      perms = [
        sortperm([sites_ind[idx]...])[1:(counter + repeat + repeat_next[idx] + 1)]
        for idx in 1:length(sites_ind)
      ]

      op_inds_next = unique([
        (inds_ord[idx], repeat_next[idx], perms[idx]) for idx in 1:length(sites_ind)
      ])
      op_inds_next = sort(op_inds_next, by = x -> x[1])

      # ensure same perm prefix
      op_inds_next = op_inds_next[findall(
        x -> x[3][1:length(perm_ind)] == perm_ind, op_inds_next
      )]

      next_site = op_inds_next[1][1]

      # contract intermediate sites to next operator site
      for str in (site + 1):(next_site - 1)
        if jw_next % 2 != 0
          L = L * apply(op("F", s[str]), psi[str]) * onebra[str]
        else
          L = L * psi[str] * onebra[str]
        end
      end

      add_operator_fermi(
        op_inds_next,
        sites_ind,
        L,
        counter + repeat + 1,
        element,
        N,
        ops,
        s,
        psi,
        onebra,
        C,
        indices,
        jw_next,
      )
    end

    # advance cached L_prev to next possible site in this loop
    if (a < length(op_inds))
      next_site = op_inds[a + 1][1]
      for str in site:(next_site - 1)
        if jw % 2 != 0
          L_prev = L_prev * apply(op("F", s[str]), psi[str]) * onebra[str]
        else
          L_prev = L_prev * psi[str] * onebra[str]
        end
      end
    end
  end
end

end # module
