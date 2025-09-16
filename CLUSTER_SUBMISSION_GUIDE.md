# Edwards-Anderson Binder Parameter Cluster Submission Guide

## Quick Start

1. **Connect to submit machine:**
   ```bash
   ssh qia.wang@ap40.uw.osg-htc.org
   cd /home/qia.wang/binder-parameter-OSpool/jobs
   ```

2. **Submit test jobs (RECOMMENDED FIRST):**
   ```bash
   condor_submit jobs_adaptive_test.submit
   ```

3. **Monitor jobs:**
   ```bash
   condor_q -nobatch          # Check status
   watch -n 30 'condor_q'     # Auto-refresh
   ```

4. **After test success, submit production:**
   ```bash
   condor_submit jobs_adaptive_production.submit
   ```

## Job Details

### Test Jobs (9 jobs)
- **File**: `jobs_adaptive_test.submit`
- **Parameters**: `params_adaptive_test.txt`
- **Config**: L=8,12,16 × λ=0.3,0.5,0.7 × 1 sample
- **Trials**: 200 per job
- **Resources**: 2 CPU, 2GB RAM
- **Purpose**: Verify L=16 works with fixed BinderSim.jl

### Production Jobs (153 jobs)  
- **File**: `jobs_adaptive_production.submit`
- **Parameters**: `params_adaptive.txt`
- **Config**: L=8,12,16 × 17 λ values × 3 samples
- **Trials**: 2000 per job
- **Resources**: 4 CPU, 6GB RAM
- **Purpose**: Full simulation study

## Monitoring Commands

```bash
# Job status
condor_q -nobatch                    # All jobs
condor_q -run                        # Running jobs only
condor_q -hold                       # Held jobs only

# Job analysis
condor_q -analyze                    # Why jobs aren't running
condor_q -better-analyze             # Detailed analysis

# Logs and output
ls -la logs/                         # Check log files
tail -f logs/[cluster].[process].out # Follow specific job
ls -la output/                       # Check results

# Cancel jobs if needed
condor_rm $(whoami)                  # Cancel all your jobs
condor_rm [cluster]                  # Cancel specific cluster
```

## What's Different Now

✅ **Fixed BinderSim.jl**: Now works for L>12 using the same approach as your notebook
✅ **Robust correlators**: Uses direct correlator() function instead of complex chunking  
✅ **Central region**: Analyzes central 60% of sites for efficiency
✅ **Adaptive parameters**: Automatically adjusts maxdim/cutoff based on system size

## Expected Runtime

- **Test jobs**: ~5-10 minutes each (200 trials)
- **Production jobs**: ~30-60 minutes each (2000 trials)
- **Total test time**: ~1 hour for all 9 jobs
- **Total production time**: ~5-10 hours for all 153 jobs

## Results Structure

```
output/
├── adaptive_test_L8_lam0.3_s1.json     # Test results
├── L8_lam0.1_s1.json                   # Production results
├── L8_lam0.1_s2.json
├── ...
└── L16_lam0.9_s3.json
```

Each JSON file contains:
- Binder parameter value
- Statistical information  
- Runtime metadata
- System parameters
