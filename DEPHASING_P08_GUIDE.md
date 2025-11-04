# Dephasing P=0.8 Job Submission Guide

## Files Created

### Main Simulation Files
1. **run_dephasing_p08.jl** - Julia script for P=0.8 dephasing simulations
   - Sets P_x = P_zz = 0.8
   - Otherwise identical to P=0.2 version

### HTCondor Job Files (in `jobs/` directory)
2. **run_dephasing_p08.sh** - Wrapper script for HTCondor
   - Creates output directory
   - Installs Julia packages
   - Runs the simulation

3. **dephasing_p08_jobs.submit** - HTCondor submit file
   - Configured for OSPool with your container
   - 4 CPUs, 8 GB RAM per job

4. **params_dephasing_p08.txt** - Parameter file with 153 jobs
   - 3 system sizes: L = 8, 12, 16
   - 17 lambda values: 0.1-0.9 (dense sampling around 0.46-0.54)
   - 3 samples per parameter combination
   - Seeds: 5001-5153

5. **collect_dephasing_p08_results.sh** - Results collection script
   - Run after jobs complete to gather all JSON files

## Job Statistics
- Total jobs: **153**
- Seeds: 5001 to 5153
- Expected output files: 153 JSON files
- Output prefix format: `dephasing_p08_L{L}_lam{lambda}_P0.8_s{sample}`

## Submission Commands (on cluster)

```bash
# Navigate to jobs directory
cd jobs/

# Submit the jobs
condor_submit dephasing_p08_jobs.submit

# Check job status
condor_q

# After jobs complete, collect results
./collect_dephasing_p08_results.sh

# Create tarball for transfer (adjust the date)
tar czf results_dephasing_p08_2024-11-03.tar.gz results_dephasing_p08/

# Transfer to Mac using magic wormhole
wormhole send results_dephasing_p08_2024-11-03.tar.gz
```

## Expected Directory Structure After Jobs Complete

```
jobs/
├── dephasing_p08_L8_lam0.10_P0.8_s1/
│   └── output/
│       └── dephasing_p08_L8_lam0.10_P0.8_s1.json
├── dephasing_p08_L8_lam0.10_P0.8_s2/
│   └── output/
│       └── dephasing_p08_L8_lam0.10_P0.8_s2.json
├── ... (151 more job directories)
└── results_dephasing_p08/
    ├── dephasing_p08_L8_lam0.10_P0.8_s1.json
    ├── dephasing_p08_L8_lam0.10_P0.8_s2.json
    └── ... (153 total result files)
```

## Comparison with Previous Runs

| Parameter | P=0.01 | P=0.2 | P=0.8 (NEW) |
|-----------|--------|-------|-------------|
| Number of jobs | 150 | 153 | 153 |
| Seeds | 1001-1150 | 3001-3153 | 5001-5153 |
| System sizes | L=8,12,16 | L=8,12,16 | L=8,12,16 |
| Lambda values | 17 values | 17 values | 17 values |
| Samples | 3 per point | 3 per point | 3 per point |

## Notes
- The P=0.8 jobs use the same parameter scan as P=0.2
- Output format is consistent with previous runs
- All necessary setup (mkdir, Pkg.instantiate) is included in the wrapper script
- Results can be analyzed together with P=0.01 and P=0.2 data
