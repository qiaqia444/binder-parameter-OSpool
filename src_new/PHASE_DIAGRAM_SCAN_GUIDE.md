# Phase Diagram Scan Guide

## Understanding the Phase Diagram

```
     P (dephasing)
     ^
  1  |________________________
     |        AREA            |
     |       REGIME          |
     |                        |
  Pc |----LEFT BOUNDARY------|---
     |                        |
     |    VOLUME-LAW          |
     |     ENTANGLED          |
  0  |________________________|_____> λ (measurement)
     0         λc             1
     
     BOTTOM BOUNDARY
```

## Three Different Scans

### 1. Bottom Boundary (P = 0, vary λ)
**Physics:** Pure measurement-induced phase transition

**Correct Method:** ✅ Pure state quantum trajectories
```julia
# Use current ea_binder_mc (NO dephasing)
result = ea_binder_mc(L; lambda_x=λ, lambda_zz=λ, ntrials=1000)
```

**Why it works:** No dephasing → states remain pure → trajectory unraveling is exact

---

### 2. Left Boundary (λ = fixed, vary P)  
**Physics:** Dephasing-induced transition with measurement background

**WRONG Method:** ❌ Pure state quantum trajectories (current approach)
```julia
# WRONG - uses stochastic dephasing on pure states!
result = ea_binder_mc_dephasing(L; lambda_x=0.5, lambda_zz=0.5, P_x=P, P_zz=P, ntrials=1000)
```

**Correct Method:** ✅ Density matrix evolution
```julia
# CORRECT - uses DiagonalStateMPS with deterministic dephasing
include("src_new/dynamics_density_matrix.jl")
result = ea_binder_density_matrix(L; lambda_x=0.5, lambda_zz=0.5, P_x=P, P_zz=P, ntrials=1000)
```

**Why density matrices are required:**
- Dephasing creates classical mixtures: ρ → (1-P)ρ + P·X·ρ·X
- Cannot be represented as single pure state
- Trajectory approach samples which Kraus → WRONG physics!

---

### 3. Interior Points (vary both λ and P)
**Correct Method:** ✅ Density matrix evolution
```julia
result = ea_binder_density_matrix(L; lambda_x=λ, lambda_zz=λ, P_x=P, P_zz=P, ntrials=1000)
```

---

## What You're Doing Now

Based on your description "left boundary where dephasing occurs but X measurement is at some strength":

**Your Goal:** Scan P from 0 to 1 at fixed λ (e.g., λ = 0.5)

**Current Code:** `run_dephasing_p05.jl` with `ea_binder_mc_dephasing`
- Sets P_x = 0.5, P_zz = 0.5 (FIXED)
- Varies lambda_x (scanning λ axis)
- Uses pure state trajectories

**This is actually the WRONG scan direction!** You're scanning λ at fixed P, not P at fixed λ.

---

## How to Fix Your Scans

### For Left Boundary Scan (Fixed λ, Vary P)

Create new file: `run_left_boundary_scan.jl`

```julia
using JSON
using Random

include("src_new/dynamics_density_matrix.jl")

function main()
    # Fixed measurement strength (on left boundary)
    lambda_x = 0.5
    lambda_zz = 0.5
    
    # Scan dephasing from 0 to 1
    P_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    L = 12
    ntrials = 1000
    
    results = []
    for P in P_values
        println("\nRunning P = $P (λ = $lambda_x)")
        
        # Use density matrix evolution (CORRECT for dephasing)
        result = ea_binder_density_matrix(
            L;
            lambda_x = lambda_x,
            lambda_zz = lambda_zz,
            P_x = P,
            P_zz = P,
            ntrials = ntrials,
            maxdim = 256,
            seed = 42
        )
        
        push!(results, Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "P_x" => P,
            "P_zz" => P,
            "B" => result.B,
            "B_mean" => result.B_mean_of_trials,
            "B_std" => result.B_std_of_trials,
            "S2" => result.S2_bar,
            "S4" => result.S4_bar
        ))
    end
    
    # Save results
    output_file = "left_boundary_L$(L)_lambda$(lambda_x).json"
    open(output_file, "w") do f
        JSON.print(f, results, 4)
    end
    
    println("\n✓ Results saved to: $output_file")
end

main()
```

### For Bottom Boundary Scan (P = 0, Vary λ)

Your existing code is fine:
```julia
# Use ea_binder_mc (no dephasing version)
result = ea_binder_mc(L; lambda_x=λ, lambda_zz=λ, ntrials=1000)
```

---

## Summary

| Scan Type | P value | λ value | Correct Method |
|-----------|---------|---------|----------------|
| **Bottom** | 0 (fixed) | vary | `ea_binder_mc` (pure states) |
| **Left** | vary | fixed | `ea_binder_density_matrix` (density matrix) |
| **Interior** | vary | vary | `ea_binder_density_matrix` (density matrix) |

**Critical:** Any scan with P > 0 MUST use density matrix evolution, not pure state trajectories!
