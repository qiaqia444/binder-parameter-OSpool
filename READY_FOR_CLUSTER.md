# Ready for Cluster Deployment

## ✅ Cleanup Complete

All incorrect dephasing code has been removed:
- ❌ Old simulation scripts (run_dephasing*.jl) 
- ❌ Old plotting scripts (plot_dephasing*.jl)
- ❌ Old test files (test_dephasing.jl)
- ❌ Old analysis scripts (analyze_dephasing*.jl)
- ❌ Old result directories and archives
- ❌ Old cluster submission scripts
- ❌ Old parameter files

## ✓ Correct Implementation

### Core Physics (src_new/)
- `types.jl` - State type wrappers (PureStateMPS, DiagonalStateMPS, MixedStateMPS)
- `channels.jl` - Proper dephasing channels: ρ → (1-P)ρ + P·M·ρ·M
- `dynamics_density_matrix.jl` - Correct density matrix evolution

### Simulation Scripts
- `run_left_boundary_scan.jl` - Full-featured scan with ArgParse
- `run_left_boundary_simple.jl` - Simple version (edit parameters at top)

### Testing
- `test_density_matrix.jl` - ✅ All tests pass

### Documentation  
- `DENSITY_MATRIX_VS_PURE_STATE.md` - Theory explanation
- `LEFT_BOUNDARY_PHYSICS.md` - Boundary scan physics

## 🎯 Left Boundary Physics

**Goal**: Compare X measurements (λ_x) vs X dephasing (P_x)

**Parameters**:
- Fixed: λ_zz = 0.0, P_zz = 0.0 (no ZZ effects)
- Scan: P_x from 0 to 1
- Measurement: λ_x = 0.3 (fixed)

**Why this is clean**:
- Pure competition: X measurements vs X dephasing
- λ_zz=0 means no ZZ measurements to complicate things
- P_zz=0 means only X dephasing channel active
- Clear phase boundary expected

## 🚀 Cluster Deployment Plan

1. **Test locally** (already done ✓):
   ```julia
   julia test_density_matrix.jl  # All tests pass
   ```

2. **Push to GitHub**:
   ```bash
   git add src_new/ run_left_boundary_*.jl test_density_matrix.jl *.md
   git commit -m "Correct density matrix implementation for dephasing"
   git push
   ```

3. **On cluster**:
   ```bash
   git pull
   julia run_left_boundary_scan.jl --L 12 --lambda_x 0.3 --lambda_zz 0.0 --P_zz_mode zero --trials 100
   ```

## 📊 Typical Command Examples

### Small test (L=8, quick):
```bash
julia run_left_boundary_simple.jl
```

### Production run (L=12, many trials):
```bash
julia run_left_boundary_scan.jl --L 12 --lambda_x 0.3 --lambda_zz 0.0 --P_zz_mode zero --trials 100 --steps 50
```

### Large system (L=16):
```bash  
julia run_left_boundary_scan.jl --L 16 --lambda_x 0.3 --lambda_zz 0.0 --P_zz_mode zero --trials 200 --steps 100 --maxdim 64
```

## ⚠️ Important Notes

- **Never run L > 10 on Mac** - use cluster for large systems
- Always use `--P_zz_mode zero` for left boundary scan
- Results saved to `left_boundary_results/` directory
- Each run creates timestamped output file

## 📈 Expected Results

- P_x = 0: Measurement-induced entanglement (low B)
- P_x = 1: Area-law phase (high B) 
- Phase transition somewhere in between
- Should see clear crossing in B vs P_x plots
