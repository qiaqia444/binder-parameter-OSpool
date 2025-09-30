# Bipartite Entropy Analysis Workflow

This document describes the complete workflow for running bipartite entropy calculations on the cluster and analyzing results on your Mac.

## Overview

The simulation calculates bipartite entanglement entropy using weak measurements with complementary strengths:
- **System sizes**: L = [8, 12, 16]  
- **Parameter relationship**: λₓ = λ, λ_zz = 1 - λ
- **Lambda range**: 0.0 to 1.0 in 0.1 steps
- **Trials**: 500 per job for statistical accuracy
- **Total jobs**: 33 (11 lambda values × 3 system sizes)

## Step-by-Step Workflow

### 1. Submit Jobs to Cluster

```bash
# On cluster (OSPool)
cd binder-parameter-OSpool/
condor_submit jobs/bipartite_jobs.submit
```

Monitor progress:
```bash
condor_q  # Check job status
```

### 2. Collect Results (On Cluster)

After jobs complete:
```bash
# On cluster
./collect_bipartite_results.sh
```

This creates:
- `bipartite_results_YYYYMMDD_HHMM/` - Directory with all results
- `bipartite_results_YYYYMMDD_HHMM.tar.gz` - Compressed archive

### 3. Transfer Data to Mac

Using magic-wormhole:
```bash
# On cluster
wormhole send bipartite_results_YYYYMMDD_HHMM.tar.gz

# On Mac
wormhole receive
```

### 4. Analyze Results (On Mac)

```bash
# On Mac
./analyze_bipartite_workflow.sh bipartite_results_YYYYMMDD_HHMM.tar.gz
```

This automatically:
- Extracts the archive
- Loads and processes all JSON results
- Creates plots similar to your reference
- Generates analysis reports

## Generated Files

### Results Data
- `bipartite_L*_lambda*.json` - Individual job results
- `collection_summary.txt` - Summary of collected data

### Analysis Output
- `bipartite_entropy_analysis_YYYYMMDD_HHMM.png` - Main entropy vs λ plot
- `entropy_scaling_YYYYMMDD_HHMM.png` - Entropy scaling with system size  
- `bipartite_analysis_report_YYYYMMDD_HHMM.txt` - Detailed analysis report

## Expected Results

The plot should show:
- **Peak entropy** around λ = 0.5-0.6 (balanced measurement regime)
- **System size dependence** with larger L showing higher peak entropy
- **Sharp transition** at λ = 0.8-1.0 where entropy drops rapidly
- **Area law behavior** at λ = 0 and λ = 1 (measurement-dominated regimes)

## File Structure

```
binder-parameter-OSpool/
├── src/BipartiteEntropy.jl          # Core simulation module
├── run_bipartite.jl                 # Main run script  
├── jobs/
│   ├── bipartite_jobs.submit        # HTCondor submission file
│   ├── run_bipartite.sh             # Job runner script
│   ├── params_bipartite.txt         # Parameter file (33 jobs)
│   └── make_bipartite_params.jl     # Parameter generator
├── collect_bipartite_results.sh     # Results collection script
├── analyze_bipartite.jl             # Analysis script for Mac
└── analyze_bipartite_workflow.sh    # Complete workflow script
```

## Troubleshooting

### If jobs fail:
1. Check log files: `logs/bipartite_*.err`
2. Look for `bipartite_*_FAILED.json` files
3. Check memory/disk requirements in submit file

### If analysis fails:
1. Verify Julia packages: `Plots`, `JSON`, `Statistics`, `DataFrames`, `Glob`
2. Check that result files contain valid data
3. Ensure file paths are correct

### If plots look wrong:
1. Check if all jobs completed successfully
2. Verify parameter ranges in results
3. Look at analysis report for data summary

## Manual Analysis

If you need to customize analysis:

```julia
# Load results manually
using JSON, Plots
results = []
for file in glob("bipartite_L*.json", "results_directory")
    push!(results, JSON.parsefile(file))
end

# Extract data for plotting
# ... custom analysis code ...
```

## Contact

For questions about the simulation setup or analysis workflow, refer to the main project documentation or check the code comments in the Julia modules.