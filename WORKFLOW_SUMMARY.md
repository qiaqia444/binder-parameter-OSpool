# Forced Measurement Binder Parameter Calculation - Complete Workflow

## Job Setup Summary
- **Total Jobs**: 153 (51 per system size L=8,12,16)
- **Parameter Space**: 17 λ values from 0.1 to 0.9, with 3 samples each
- **Container**: `/ospool/ap40/data/qia.wang/container.sif` (Julia 1.11 + ITensors)
- **Output Format**: Individual JSON files with complete metadata

## File Naming Scheme
Output files follow the pattern: `forced_L{L}_lam{lambda}_s{sample}.json`

Examples:
- `forced_L8_lam0.1_s1.json` - L=8, λ=0.1, sample 1
- `forced_L12_lam0.46_s3.json` - L=12, λ=0.46, sample 3  
- `forced_L16_lam0.9_s2.json` - L=16, λ=0.9, sample 2

## Job Submission Commands

```bash
# Submit forced measurement jobs
cd /path/to/binder-parameter-OSpool/jobs
condor_submit jobs_forced.submit

# Monitor job progress  
condor_q [your_username]
watch condor_q [your_username]

# Check job status
condor_q -analyze [job_id]  # If jobs are held
```

## Results Collection (On OSPool)

After all jobs complete, run the collection script:
```bash
cd /path/to/binder-parameter-OSpool
./collect_results.sh
```

This will:
1. Organize results by system size (L=8,12,16)
2. Check for failed jobs
3. Create summary report
4. Generate compressed archive for transfer

## Data Transfer to Mac

Use Magic Wormhole for secure transfer:
```bash
# On OSPool (after running collect_results.sh)
wormhole send forced_binder_results_YYYYMMDD_HHMM.tar.gz

# On Mac
wormhole receive [code-from-wormhole]
tar -xzf forced_binder_results_*.tar.gz
```

## Analysis on Mac

Run the analysis script in the extracted directory:
```bash
cd path/to/extracted/results
julia analyze_results.jl
```

This generates:
- `binder_vs_lambda_forced.png` - Main physics plot
- `data_coverage_forced.png` - Data quality visualization  
- `processed_binder_results_forced.csv` - Statistical summary

## Expected Results Structure

```
results/
├── forced_measurements/
│   ├── L8/          # 51 JSON files  
│   ├── L12/         # 51 JSON files
│   └── L16/         # 51 JSON files
├── analysis/
│   ├── collection_summary.txt
│   └── *_FAILED.json (if any)
└── [analysis outputs after running analyze_results.jl]
```

## Each JSON Output Contains
- **Parameters**: L, lambda_x, lambda_zz, lambda, seed, sample
- **Results**: binder_parameter, S2_bar, S4_bar, statistics
- **Metadata**: maxdim, cutoff, measurement_type, success status
- **Computational**: ntrials_completed, execution details

## Physics Interpretation
- **Binder Parameter B = 1 - S4/(3*S2²)** from forced +1 measurements
- **Critical Point**: Look for B ≈ 0.6 crossing between different L values
- **Finite-Size Scaling**: B(L,λ) behavior near criticality
- **Comparison**: Results complement standard quantum trajectory simulations

## Troubleshooting

If jobs fail:
1. Check logs in `logs/forced_*.err` files
2. Verify container accessibility: `singularity exec /ospool/ap40/data/qia.wang/container.sif julia --version`
3. Check parameter format in `params_adaptive.txt`
4. Monitor resource usage in HTCondor logs

## Quality Checks

Expected completion:
- **153 total jobs** (51 per L value)
- **17 λ values** from 0.1 to 0.9  
- **3 samples** per (L,λ) combination
- **Coverage**: Dense sampling near critical region (λ ≈ 0.5)

The results will provide high-precision Binder parameter data for identifying the measurement-induced phase transition in the Edwards-Anderson model.