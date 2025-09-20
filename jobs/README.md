# Jobs Directory - MIPT Simulations

This directory contains organized job submission files for MIPT simulations.

## Directory Structure

```
jobs/
├── params/                # Generated parameter files  
├── make_params.jl        # Unified parameter generator
├── jobs_template.submit  # Job submission template
├── run.sh               # Execution script
├── archive/             # Old files (archived)
└── README.md            # This file
```

## Unified Parameter Generator

**Single script handles all cases**: `make_params.jl`

### Usage
```bash
julia make_params.jl [mode] [output_file]
```

### Modes

| Mode | Description | Jobs | Runtime | Use Case |
|------|-------------|------|---------|----------|
| `test` | Quick testing | 18 | ~36 min | Development & validation |
| `critical` | Critical region | 243 | ~8 hours | Finite-size scaling |
| `standard` | Standard coverage | 255 | ~8.5 hours | Phase diagram |
| `production` | Full production | 1,840 | ~61 hours | Publication quality |

### Examples
```bash
# Quick test
julia make_params.jl test

# Critical region analysis  
julia make_params.jl critical params/critical.txt

# Standard phase diagram
julia make_params.jl standard

# Production run
julia make_params.jl production params/prod.txt
```

## Quick Start Workflow

### 1. Generate Parameters
```bash
cd jobs
julia make_params.jl test    # Start with test mode
```

### 2. Setup Job Submission
```bash
cp jobs_template.submit jobs.submit
# Edit jobs.submit to uncomment desired parameter file line
```

### 3. Submit Jobs
```bash
condor_submit jobs.submit
```

## Parameter Details

### Test Mode
- **Systems**: L = [8, 12]
- **Lambda**: [0.3, 0.5, 0.7] 
- **Samples**: 3 per point
- **Focus**: Quick validation

### Critical Mode  
- **Systems**: L = [8, 12, 16]
- **Lambda**: 0.46-0.54 (fine sampling around λ = 0.5)
- **Samples**: 10 per point
- **Focus**: Critical point analysis

### Standard Mode
- **Systems**: L = [8, 12, 16] 
- **Lambda**: 0.1-0.9 with critical region detail
- **Samples**: 5 per point
- **Focus**: Complete phase diagram

### Production Mode
- **Systems**: L = [8, 12, 16, 20]
- **Lambda**: Comprehensive coverage
- **Samples**: 20 per point  
- **Focus**: High-statistics results

## File Outputs

- **Parameter files**: `params/*.txt`
- **Job results**: `../output/L{L}_lam{λ}_s{s}.json`
- **Job logs**: `logs/job_*.{log,out,error}`

## Tips

1. **Always test first**: Start with `test` mode
2. **Monitor jobs**: Check `condor_q` and log files
3. **Critical region**: Use for crossing point analysis
4. **Production**: Only for final high-quality results

## Dependencies

- Julia with ITensors, ITensorMPS packages
- HTCondor cluster access
- Sufficient storage in `../output/`
