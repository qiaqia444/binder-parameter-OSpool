# Learning-to-Trivial Lambda0.7 - Cluster Deployment Guide

## Files Ready ✓

All files for the learning_to_trivial_lambda0.7 pipeline have been created and pushed to GitHub:

### New Core Module
- **`src_new/renyi2_dynamics.jl`** - Independent Rényi-2 Binder dynamics pipeline
  - Standalone evolution function with fast T_max=2*L
  - MPO-based moment calculations
  - Ensemble Binder with full diagnostics
  - No dependency on old EA Binder code

### Main Calculation
- **`run_learning_to_trivial_lambda0.7_scan.jl`** - Single-point computation
  - Usage: `julia run_learning_to_trivial_lambda0.7_scan.jl L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix`
  - Example: `julia run_learning_to_trivial_lambda0.7_scan.jl 8 0.49 0.21 0.3 0.3 100 12345 1 output`

### Validation
- **`test_learning_to_trivial_lambda0.7.jl`** - Quick test (L=8, P_x=0.3, 50 trials)
  - ✓ PASSED: B=0.1228, purity=0.02, valid=100%
  - Runtime: ~2 minutes

### HTCondor Batch Infrastructure
- **`jobs/jobs_learning_to_trivial_lambda0.7.submit`** - HTCondor job descriptor
- **`jobs/run_learning_to_trivial_lambda0.7.sh`** - Job wrapper script
- **`jobs/make_params_learning_to_trivial_lambda0.7.jl`** - Parameter generator
- **`jobs/params_learning_to_trivial_lambda0.7.txt`** - 2,200 job parameters

### Support Scripts
- **`submit_learning_to_trivial_lambda0.7.sh`** - Batch submission wrapper
- **`collect_learning_to_trivial_lambda0.7_results.sh`** - Results collection
- **`analyze_learning_to_trivial_lambda0.7.jl`** - Post-processing analysis

## Configuration
- **Parameters**: λ_x = 0.49 (fixed), λ_zz = 0.21 (fixed)
- **Scan**: P_x ∈ [0, 0.5] in 11 steps
- **System Sizes**: L ∈ [8, 10, 12, 14, 16]
- **Total Jobs**: 5 L × 11 P × 40 samples = 2,200 jobs
- **Trials per Job**: 100
- **Time Limit**: 24 hours per job (conservative)
- **Resources**: 4 CPUs, 8 GB RAM per job

## GitHub Status

✓ Repository: https://github.com/qiaqia444/binder-parameter-OSpool.git
✓ Latest Commit: 5c241be - "Add learning_to_trivial_lambda0.7 complete pipeline with new renyi2_dynamics"

## Cluster Deployment Instructions

### 1. SSH to Cluster
```bash
ssh <cluster>
```

### 2. Clone/Pull Latest Code
**First time:**
```bash
git clone https://github.com/qiaqia444/binder-parameter-OSpool.git
cd binder-parameter-OSpool
```

**Already cloned:**
```bash
cd binder-parameter-OSpool
git pull origin main
```

### 3. Generate Parameters (if not already done)
```bash
cd jobs
julia make_params_learning_to_trivial_lambda0.7.jl
cd ..
```

This creates `jobs/params_learning_to_trivial_lambda0.7.txt` with 2,200 job specifications.

### 4. Submit Batch Jobs
```bash
./submit_learning_to_trivial_lambda0.7.sh
```

Or manually:
```bash
cd jobs
condor_submit jobs_learning_to_trivial_lambda0.7.submit
cd ..
```

### 5. Monitor Progress
```bash
condor_q                    # All jobs
condor_q -name <schedd>     # Specific scheduler
condor_q | grep lambda0.7   # Filter to lambda0.7 jobs
```

### 6. Collect Results (After All Jobs Complete)
```bash
./collect_learning_to_trivial_lambda0.7_results.sh
```

Creates:
- `learning_to_trivial_lambda0.7_results_YYYYMMDD_HHMM/` - directory with all results
- `learning_to_trivial_lambda0.7_results_YYYYMMDD_HHMM.tar.gz` - compressed archive

### 7. Analyze Results
```bash
julia analyze_learning_to_trivial_lambda0.7.jl
```

Outputs statistics by (L, P_x) pair.

## Quick Test on Cluster (Before Full Batch)
```bash
# Test single configuration
julia run_learning_to_trivial_lambda0.7_scan.jl 8 0.49 0.21 0.1 0.1 10 42 1 test_output

# Should complete in ~25 seconds with 10 trials
# Check output/test_output.json
```

## Physics Summary

**Measurement-Induced Phase Transition**
- λ_x = 0.49: weak X measurements
- λ_zz = 0.21: weak ZZ measurements  
- P_x, P_zz: dephasing probabilities

**Observable**: Rényi-2 Binder B = 1 - M₄/(3M₂²)
- Pure state: B → 2/3
- Maximally mixed: B → 2/(3L) ≈ 0.083 for L=8

**Expected Behavior**
- Low P: Learning phase → B close to 2/3
- High P: Trivial/mixed phase → B close to 2/(3L)
- Intermediate P: Transition region

## File Locations

```
binder-parameter-OSpool/
├── src_new/
│   └── renyi2_dynamics.jl          # Core module
├── run_learning_to_trivial_lambda0.7_scan.jl
├── test_learning_to_trivial_lambda0.7.jl
├── submit_learning_to_trivial_lambda0.7.sh
├── collect_learning_to_trivial_lambda0.7_results.sh
├── analyze_learning_to_trivial_lambda0.7.jl
└── jobs/
    ├── run_learning_to_trivial_lambda0.7.sh
    ├── make_params_learning_to_trivial_lambda0.7.jl
    ├── jobs_learning_to_trivial_lambda0.7.submit
    └── params_learning_to_trivial_lambda0.7.txt
```

## Support

For issues or modifications, refer to the inline documentation in each script.
Key parameters can be adjusted in `jobs/make_params_learning_to_trivial_lambda0.7.jl` before generating the parameter file.
