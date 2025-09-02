# Binder Parameter Simulation on OSPool
This project computes the Edwards-Anderson Binder parameter for quantum many-body systems using Monte Carlo simulations with tensor networks (ITensors.jl).

## Repository Organization

### Directory Structure
```
binder-parameter-OSpool/
├── downloaded_results/     # Compressed archives from cluster runs
├── output/                 # Extracted simulation results (JSON files)
├── logs/                   # Execution logs
├── analysis_results/       # Analysis outputs (plots, CSV statistics)
├── jobs/                   # Job submission files and parameters
├── src/                    # Source code
├── containers/             # Container definitions
└── notebooks/              # Jupyter notebooks
```

### Data Management Workflow

**1. Download Results from Cluster:**
- Compressed tar.gz files → `downloaded_results/`

**2. Extract and Organize:**
```bash
# List available archives
./extract_results.sh list

# Extract specific archive
./extract_results.sh extract L8_10_12_complete_results.tar.gz

# Extract all archives
./extract_results.sh extract-all

# Check current status
./extract_results.sh status
```

**3. Analysis:**
```bash
# Analyze JSON results and generate plots
julia analyze_json_results.jl
```

Results are saved to `analysis_results/` directory.

## Simulation Parameters (L=8,10,12 Complete)

- **System sizes (L)**: [8, 10, 12] ✅ **COMPLETED**
- **Lambda range**: 0.1 to 0.9 with fine grid around critical region (22 values total)
  - Coarse grid: λ ∈ {0.1, 0.2, 0.3, 0.35, 0.4, 0.6, 0.65, 0.7, 0.8, 0.9}
  - Fine grid: λ ∈ {0.42, 0.44, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55}
- **Samples per parameter**: 20 independent runs for statistics
- **Monte Carlo trials**: 2000 per job
- **Results**: 1,302 completed simulations using ITensorCorrelators

## Methods Comparison

### ITensorCorrelators (Completed)
- **Status**: ✅ Complete (1,302 results)
- **Method**: Built-in correlation functions from ITensorCorrelators.jl
- **Files**: `run_standard.jl`, `src/BinderSim.jl`

### Manual Correlators (Ready for Deployment)
- **Status**: 🚀 Infrastructure ready (1,320 jobs prepared)
- **Method**: Direct tensor contractions for method validation
- **Files**: `run_manual.jl`, `src/BinderSimManual.jl`, `src/ManualCorrelators.jl`

## Quick Start

### Extract and Analyze Existing Results
```bash
# Extract completed simulation data
./extract_results.sh extract-all

# Analyze and generate plots
julia analyze_json_results.jl
```

### Deploy Manual Correlators (Next Step)
```bash
# Submit manual correlator jobs to cluster
cd jobs/
condor_submit manual.submit
```

## File Structure

- `src/BinderSim.jl` - Main simulation module
- `run.jl` - Entry point for cluster jobs
- `jobs/` - HTCondor job definitions and parameter generation
- `containers/` - Singularity container definition
- `notebooks/` - Original Jupyter notebook with development code

## Dependencies

- ITensors.jl, ITensorMPS.jl, ITensorCorrelators.jl
