# Left Boundary Physics: What to Scan

## The Phase Diagram

```
     P (dephasing)
     ^
  1  |________________________
     |                       |
     |   AREA REGIME         |
     |   (Classical)         |
     |                       |
  Pc |----LEFT BOUNDARY------|---
     |                       |
     |    VOLUME-LAW          |
     |   (Entangled)          |
  0  |________________________|_____> λ (measurement)
     0         λc             1
```

## Different Left Boundary Scenarios

### **Option 1: Pure Dephasing (No Measurements)**
```bash
julia run_left_boundary_scan.jl --L 12 --lambda_x 0.0 --lambda_zz 0.0 --P_zz_mode match
```

**Physics:**
- No measurements at all
- Pure environmental decoherence
- Studies: P = 0 (coherent) → P = 1 (fully dephased)

**Phase transition:** From quantum coherent to classical mixture

---

### **Option 2: X Measurements vs X Dephasing (RECOMMENDED)**
```bash
julia run_left_boundary_scan.jl --L 12 --lambda_x 0.3 --lambda_zz 0.0 --P_zz_mode zero
```

**Physics:**
- X measurements at strength λ_x = 0.3 (trying to create entanglement)
- No ZZ measurements (λ_zz = 0)
- X dephasing varies from P_x = 0 to 1
- No ZZ dephasing (P_zz = 0)

**Competition:**
- X measurements → create quantum correlations
- X dephasing → destroy quantum coherence

**This is the CLEANEST left boundary!** Pure competition between measurement-induced entanglement and dephasing-induced decoherence in the same basis.

---

### **Option 3: Both Measurements, Both Dephasing (Symmetric)**
```bash
julia run_left_boundary_scan.jl --L 12 --lambda_x 0.5 --lambda_zz 0.5 --P_zz_mode match
```

**Physics:**
- Both X and ZZ measurements
- Both X and ZZ dephasing (matched: P_zz = P_x)
- Symmetric protocol

**More complex:** Multiple competing effects

---

### **Option 4: X Measurements, Asymmetric Dephasing**
```bash
julia run_left_boundary_scan.jl --L 12 --lambda_x 0.3 --lambda_zz 0.0 --P_zz_mode fixed --P_zz_value 0.1
```

**Physics:**
- X measurements only
- X dephasing varies (P_x)
- Fixed weak ZZ dephasing (P_zz = 0.1)

**Studies:** How cross-basis dephasing affects transitions

---

## My Recommendation

For **cleanest physics** and **easiest interpretation**, use **Option 2**:

```bash
# Small test (L=8, quick)
julia run_left_boundary_scan.jl \
    --L 8 \
    --lambda_x 0.3 \
    --lambda_zz 0.0 \
    --P_steps 11 \
    --P_zz_mode zero \
    --ntrials 200

# Full simulation (L=12, production)
julia run_left_boundary_scan.jl \
    --L 12 \
    --lambda_x 0.3 \
    --lambda_zz 0.0 \
    --P_steps 21 \
    --P_zz_mode zero \
    --ntrials 1000
```

## Why λ_x = 0.3?

- **Too small (λ < 0.2):** Weak measurement effects, hard to see phase transition
- **Sweet spot (λ ≈ 0.3-0.4):** Strong enough to create correlations, but still competes with dephasing
- **Too large (λ > 0.5):** Measurements dominate, less interesting competition

You can scan different λ_x values (0.2, 0.3, 0.4, 0.5) to see how the critical P_c changes!

## Expected Physics

At fixed λ_x:
- **P_x ≈ 0:** Volume-law entanglement (measurements win)
- **P_x ≈ P_c:** Phase transition
- **P_x ≈ 1:** Area-law (dephasing wins, classical state)

The critical P_c should **increase** with λ_x (stronger measurements resist dephasing better).
