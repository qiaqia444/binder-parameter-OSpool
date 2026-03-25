# Memory-to-Trivial Transition Scan

## Physics Configuration

This simulation suite explores the **memory-to-trivial phase transition** with fixed measurement strengths while scanning dephasing probabilities.

### Parameters

- **λ_x (X measurement strength)**: 0.1 (FIXED)
- **λ_zz (ZZ measurement strength)**: 0.7 (FIXED)
- **P_x = P_zz (dephasing probability)**: 0.0 → 0.5 (SCANNED)
- **System sizes**: L = 8, 10, 12, 14, 16 (finite-size scaling)
- **Trials per job**: 100
- **Samples per configuration**: 40
- **Total jobs**: 2,200

### Naming Convention

Output files follow the pattern:
```
memory_to_trivial_L{L}_lx{λ_x:.2f}_lzz{λ_zz:.2f}_P{P:.2f}_s{sample}
```

Example: `memory_to_trivial_L12_lx0.10_lzz0.70_P0.25_s5.json`

## File Structure

### Main Simulation Script
- **[run_memory_to_trivial_scan.jl](run_memory_to_trivial_scan.jl)** — Julia script that runs the actual simulation
  - Uses density matrix evolution (DiagonalStateMPS)
  - Handles parameter arguments from cluster jobs
  - Outputs JSON results

### Parameter Generation
- **[jobs/make_params_memory_to_trivial.jl](jobs/make_params_memory_to_trivial.jl)** — Generates parameter combinations
  - Creates `params_memory_to_trivial.txt` with all job configurations
  - Organized by system size L, then dephasing probability P

### Generated Parameter File
- **[jobs/params_memory_to_trivial.txt](jobs/params_memory_to_trivial.txt)** — List of all 2,200 parameter sets
  - Format: `L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix`

### Cluster Submission
- **[jobs/jobs_memory_to_trivial.submit](jobs/jobs_memory_to_trivial.submit)** — HTCondor submit file
  - Requests 4 CPUs, 8GB memory per job
  - Transfers `src_new/`, `Project.toml`, `Manifest.toml`
  - Uses container: `/ospool/ap40/data/qia.wang/container.sif`

### Job Execution
- **[jobs/run_memory_to_trivial.sh](jobs/run_memory_to_trivial.sh)** — Bash script executed by HTCondor
  - Sets up Julia environment
  - Runs `run_memory_to_trivial_scan.jl` with specified parameters
  - Outputs JSON to `output/` directory

### Submission Commands

**Generate parameters** (if needed to regenerate):
```bash
cd jobs
julia make_params_memory_to_trivial.jl
cd ..
```

**Submit jobs to cluster**:
```bash
./submit_memory_to_trivial.sh
```

Or manually:
```bash
cd jobs
condor_submit jobs_memory_to_trivial.submit
cd ..
```

### Results Collection

**After jobs complete** (on cluster):
```bash
./collect_memory_to_trivial_results.sh
```

This will:
1. Create timestamped results directory with L subdirectories
2. Collect all output JSON files, organized by system size
3. Check for failures and archive them separately
4. Create compressed tar.gz archive
5. Provide magic wormhole transfer instructions

**To transfer to Mac**:
```bash
# On cluster after collection
wormhole send memory_to_trivial_results_YYYYMMDD_HHMM.tar.gz

# On Mac
wormhole receive
tar -xzf memory_to_trivial_results_YYYYMMDD_HHMM.tar.gz
```

## Physics Interpretation

The memory-to-trivial transition explores:
- **Memory phase** (low P): System retains entanglement and quantum coherence
- **Trivial phase** (high P): Dephasing dominates, system becomes classical

With λ_x = 0.1 and λ_zz = 0.7, we examine this transition in a regime with significant ZZ measurements and weak X measurements.

The **Binder parameter** as a function of P reveals:
- Critical behavior near the phase transition
- System size dependence (finite-size scaling)
- Order of the transition (first vs second order)

## Monitoring Jobs

**On cluster** (while jobs are running):
```bash
condor_q          # Check job status
condor_rm <job>   # Cancel specific job if needed
tail -f jobs/logs/memory_to_trivial_*.out  # Monitor output
```

## Expected Output

Each job produces a JSON file with:
- Metadata (parameters, timestamp, system config)
- Binder parameter measurements
- Statistics (mean, std) across Monte Carlo trajectories

Results structure:
```
memory_to_trivial_results_YYYYMMDD_HHMM/
├── L8/
│   ├── memory_to_trivial_L8_lx0.10_lzz0.70_P0.00_s1.json
│   ├── memory_to_trivial_L8_lx0.10_lzz0.70_P0.00_s2.json
│   └── ... (40 files per P value)
├── L10/
├── L12/
├── L14/
├── L16/
└── failed/  (if any jobs failed)
```
