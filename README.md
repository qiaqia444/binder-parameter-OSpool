# Binder Parameter Simulation on OSPool

This project computes the Edwards-Anderson Binder parameter for quantum many-body systems using Monte Carlo simulations with tensor networks (ITensors.jl).

## Simulation Parameters

- **System sizes (L)**: [12, 16, 20, 24, 28]
- **Lambda range**: 0.1 to 0.9 with fine grid (0.02 step) around λ=0.5
- **Time evolution**: T_max = 2L
- **Monte Carlo trials**: 1000 per job
- **Total jobs**: 255 (5 sizes × 17 λ values × 3 samples)
- **Total computational cost**: 255,000 Monte Carlo trials

## Physics

The simulation studies weak measurements on a 1D spin chain:
- Weak X measurements with strength λₓ = λ
- Weak ZZ measurements with strength λ_zz = 1-λ
- Computes Edwards-Anderson Binder parameter: B = 1 - S₄/(3S₂²)

## Quick Start

```bash
# Test locally
./run_workflow.sh test

# Generate parameter files
./run_workflow.sh params

# Submit to cluster (on OSPool)
./run_workflow.sh submit

# Monitor jobs
./run_workflow.sh status

# Collect results
./run_workflow.sh collect
```

## File Structure

- `src/BinderSim.jl` - Main simulation module
- `run.jl` - Entry point for cluster jobs
- `jobs/` - HTCondor job definitions and parameter generation
- `containers/` - Singularity container definition
- `notebooks/` - Original Jupyter notebook with development code

## Dependencies

- ITensors.jl, ITensorMPS.jl, ITensorCorrelators.jl
- Julia 1.10+ with required packages (see Project.toml)
