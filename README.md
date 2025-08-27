# Binder Parameter Simulation on OSPool
This project computes the Edwards-Anderson Binder parameter for quantum many-body systems using Monte Carlo simulations with tensor networks (ITensors.jl).

## Simulation Parameters

- **System sizes (L)**: [12, 16, 20, 24, 28] (5 values)
- **Lambda range**: 0.1 to 0.9 with fine grid around λ=0.5 (17 values total)
  - Coarse grid: λ ∈ {0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9} (8 values)
  - Fine grid: λ ∈ {0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54} (9 values)
- **Samples per parameter**: 3 independent runs for statistics
- **Time evolution**: T_max = 2L
- **Monte Carlo trials**: 1000 per job

## Quick Start

### On Cluster (OSPool)
```bash
# Test locally
./run_workflow.sh test

# Generate parameter files
./run_workflow.sh params

# Submit to cluster
./run_workflow.sh submit

# Monitor jobs
./run_workflow.sh status
```

### On Local Mac (After Jobs Complete)
```bash
# Download all results and create plots
./create_plots.sh

# Or step by step:
./download_results.sh          # Download results from cluster
julia analyze_and_plot.jl      # Create plots like original notebook
```

This will generate:
- `binder_vs_lambda_combined.png` - Main plot with all L values
- `binder_vs_lambda_L*.png` - Individual plots for each system size
- `binder_results_summary.csv` - Data table for further analysis

## File Structure

- `src/BinderSim.jl` - Main simulation module
- `run.jl` - Entry point for cluster jobs
- `jobs/` - HTCondor job definitions and parameter generation
- `containers/` - Singularity container definition
- `notebooks/` - Original Jupyter notebook with development code

## Dependencies

- ITensors.jl, ITensorMPS.jl, ITensorCorrelators.jl
