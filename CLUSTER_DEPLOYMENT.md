# Left Boundary Scan - Cluster Deployment

## ✅ Files Ready for Cluster

### Parameter Generation
- `jobs/make_params_left_boundary.jl` - Generates parameter file ✓
- `jobs/params_left_boundary.txt` - **550 jobs** generated ✓

### Cluster Submission
- `jobs/run_left_boundary.sh` - Job execution script
- `jobs/jobs_left_boundary.submit` - HTCondor submit file
- `submit_left_boundary.sh` - Main submission script

### Results Collection
- `collect_left_boundary_results.sh` - Organizes and archives results

## 📊 Simulation Parameters

**Fixed values** (left boundary physics):
- λ_x = 0.3 (X measurement strength)
- λ_zz = 0.0 (no ZZ measurements)
- P_zz = 0.0 (no ZZ dephasing)

**Scan range**:
- P_x: [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
- 11 dephasing values

**System sizes**:
- L = [8, 10, 12, 14, 16]
- 5 system sizes for finite-size scaling

**Statistics**:
- 400 trials per job
- 10 samples per configuration (for error bars)
- 550 total jobs

## 🚀 Deployment Instructions

### 1. On your Mac (before cluster)
```bash
# Verify tests pass
julia test_density_matrix.jl

# Push to GitHub
git add src_new/ run_left_boundary_scan.jl jobs/ submit_left_boundary.sh collect_left_boundary_results.sh
git commit -m "Left boundary scan: dephasing-induced phase transition"
git push
```

### 2. On the cluster
```bash
# Pull latest code
git pull

# Generate parameters (already done, but in case you need to regenerate)
cd jobs
julia make_params_left_boundary.jl
cd ..

# Submit jobs
./submit_left_boundary.sh

# Monitor jobs
condor_q

# Check logs
tail -f jobs/logs/left_boundary_condor.log
```

### 3. After jobs complete
```bash
# Collect results
./collect_left_boundary_results.sh

# This creates:
# - left_boundary_results_<timestamp>/
# - left_boundary_results_<timestamp>.tar.gz
```

## 📈 Expected Runtime

- L=8: ~2-5 min per job
- L=10: ~5-10 min per job
- L=12: ~10-20 min per job
- L=14: ~20-40 min per job
- L=16: ~40-80 min per job

**Total**: ~550 jobs × ~15 min average = **~137 CPU-hours**

With cluster parallelization: **~2-4 hours wall time**

## 📁 Output Format

Each job produces JSON file:
```json
{
  "parameters": {
    "L": 12,
    "lambda_x": 0.3,
    "lambda_zz": 0.0,
    "P_x": 0.25,
    "P_zz": 0.0,
    "trials": 400
  },
  "results": {
    "B": 0.5234,
    "B_err": 0.0145,
    "M2_squared": 0.6789,
    "M4_squared": 0.4567
  }
}
```

## 🎯 Next Steps After Data Collection

1. **Analysis script**: Parse all JSON files, compute averages
2. **Plotting script**: B vs P_x for each L, identify crossing
3. **Finite-size scaling**: Extract critical P_x^c from crossing point

Ready to run on cluster! 🚀
