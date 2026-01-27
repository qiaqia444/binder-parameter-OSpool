# Pure States vs Density Matrices: When to Use Each

## The Critical Distinction

### Pure State Trajectories (Currently in dynamics.jl)
**When to use:** Weak measurements ONLY (no dephasing)

**Physics:** 
- Each trajectory remains pure: |ψ⟩
- Stochastically sample which Kraus operator to apply
- Average over many trajectories gives ensemble density matrix

**Implementation:**
```julia
if rand() < P_x
    ψ → X|ψ⟩  # Apply X with probability P_x
else
    ψ → |ψ⟩   # Do nothing
end
```

**Correct for:** Quantum trajectories, unraveling, pure measurement protocols

---

### Density Matrix Evolution (NEW: channels.jl)
**When to use:** Dephasing or decoherence present

**Physics:**
- State is genuinely mixed: ρ (diagonal or general)
- Apply ALL Kraus operators with proper weights
- Deterministic evolution of the density matrix

**Implementation:**
```julia
ρ → (1-P_x)ρ + P_x·X·ρ·X†
```

**Correct for:** Environmental decoherence, dephasing, dissipation

---

## The Physics of Dephasing

### X Dephasing Channel
The quantum channel is:
```
ℰ_X(ρ) = (1-P_x)ρ + P_x·X·ρ·X†
```

For **diagonal states** where X²=I:
```
ℰ_X(ρ) = (1-P_x)ρ + P_x·X·ρ·X = [(1-P_x)I + P_x·X]·ρ
```

This is NOT a single pure state - it's a **classical mixture** of two outcomes!

### ZZ Dephasing Channel
```
ℰ_ZZ(ρ) = (1-P_zz)ρ + P_zz·(Z⊗Z)·ρ·(Z⊗Z)†
```

For diagonal states:
```
ℰ_ZZ(ρ) = [(1-P_zz)I + P_zz·(Z⊗Z)]·ρ
```

---

## Why This Matters for Your Simulations

### Scenario 1: Weak Measurements ONLY
✅ **Use:** PureStateMPS with quantum trajectories (current dynamics.jl)
```julia
function evolve_one_trial(L::Int; lambda_x, lambda_zz, kwargs...)
    ψ, sites = create_up_state_mps(L)  # Pure state
    # Apply weak measurements (stochastic)
    for t in 1:T_max
        for i in 1:L
            ψ = sample_and_apply(ψ, KX0[i], KX1[i], [i])
        end
    end
    return ψ
end
```

### Scenario 2: Weak Measurements + Dephasing
❌ **WRONG:** PureStateMPS with stochastic dephasing (current approach)
```julia
# This is WRONG! Each trajectory stays pure but should be mixed
if rand() < P_x
    ψ = apply X gate  # Stochastic sampling
end
```

✅ **CORRECT:** DiagonalStateMPS with deterministic channels
```julia
function evolve_one_trajectory_density_matrix(L::Int; lambda_x, lambda_zz, P_x, P_zz, kwargs...)
    # Start with diagonal state (classical mixture)
    state = zero_state(DiagonalStateMPS, L)
    
    for t in 1:T_max
        # Weak measurements (still sample outcomes)
        for i in 1:L
            state, _, _ = measure_with_outcome(state, Sx, lambda_x, i)
        end
        
        # Dephasing (deterministic, apply both Kraus operators)
        for i in 1:L
            state = apply_x_dephasing_channel(state, i, P_x)  # ρ → (1-P)ρ + P·X·ρ
        end
        
        # More measurements...
        for i in 1:(L-1)
            state, _, _ = measure_with_outcome(state, SzSz, lambda_zz, i)
        end
        
        # More dephasing
        for i in 1:(L-1)
            state = apply_zz_dephasing_channel(state, i, i+1, P_zz)
        end
    end
    
    return state
end
```

---

## Updated Binder Parameter Calculation

### For Pure State Trajectories (No Dephasing)
Keep current approach - it's correct:
```julia
function ea_binder_mc(L; lambda_x, lambda_zz, ntrials, kwargs...)
    for t in 1:ntrials
        ψ, sites = evolve_one_trial(L; lambda_x, lambda_zz)
        M2sq, M4sq = compute_correlators(ψ, sites)
        # ... accumulate statistics
    end
end
```

### For Density Matrix Evolution (With Dephasing)
Use DiagonalStateMPS:
```julia
function ea_binder_mc_dephasing_correct(L; lambda_x, lambda_zz, P_x, P_zz, ntrials, kwargs...)
    for t in 1:ntrials
        # Evolve density matrix (not pure state)
        state = evolve_density_matrix(L; lambda_x, lambda_zz, P_x, P_zz)
        
        # Extract MPS from wrapper
        ρ = get_mps(state)
        sites = siteinds(ρ)
        
        # Compute correlators on density matrix
        M2sq, M4sq = compute_correlators(ρ, sites)
        # ... accumulate statistics
    end
end
```

---

## Comparison: What Changes

| Aspect | Pure State (Old) | Density Matrix (New) |
|--------|-----------------|---------------------|
| **State Type** | MPS (pure) | DiagonalStateMPS |
| **Initial State** | `create_up_state_mps(L)` | `zero_state(DiagonalStateMPS, L)` |
| **Weak Measurements** | Sample outcome | Sample outcome (same!) |
| **Dephasing** | `if rand() < P: apply X` | `ρ → (1-P)ρ + P·X·ρ` |
| **Trajectory Nature** | Pure throughout | Mixed throughout |
| **Physical Meaning** | One realization | Ensemble-averaged |

---

## Key Insight

**The fundamental question:**
- Do you want to simulate **individual quantum trajectories** (each pure, then average)?
- Or simulate the **ensemble-averaged density matrix** directly?

**For weak measurements only:** Both give same result, trajectories are more efficient

**For weak measurements + dephasing:** 
- ❌ Pure state trajectories don't correctly represent dephasing
- ✅ Density matrix evolution is the correct approach

---

## Recommendation

1. **Keep current `dynamics.jl`** for pure measurement protocols
2. **Use new `channels.jl`** when dephasing is present
3. **Create separate evolution functions:**
   - `evolve_pure_trajectories()` - for measurement-only
   - `evolve_density_matrix()` - for measurement + dephasing

This matches the src_1 architecture where DiagonalStateMPS and MixedStateMPS were designed specifically for dephasing scenarios.
