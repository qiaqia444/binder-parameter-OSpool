# Corrected MIPT Jobs - L = [8, 12, 16] with Fixed Weak Measurements

## What was fixed:
1. **Weak measurement operators** now exactly match the sparse matrix implementation
2. **System sizes** focused on L = [8, 12, 16] for crossing point analysis  
3. **Parameter mapping**: λ_x = λ, λ_zz = 1-λ (unchanged)

## Job details:
- **Total jobs**: 153
- **System sizes**: L = 8, 12, 16
- **Lambda values**: 17 values (coarse + fine around λ = 0.5)
- **Samples per λ**: 3 independent runs
- **Trials per sample**: 2000

## Key fix:
The weak measurement operators now use the correct form:
```julia
# X weak measurement: (I + (-1)^outcome * λ * X) / √(2(1 + λ²))
# ZZ weak measurement: (I⊗I + (-1)^outcome * λ * ZZ) / √(2(1 + λ²))
```
This exactly matches your sparse matrix implementation.

## Expected result:
With this fix, you should see **crossing behavior** at λ = 0.5 where all system sizes L = 8, 12, 16 have similar Binder parameter values, rather than the monotonic decrease we saw before.

## Submission:
```bash
./submit_corrected_jobs.sh
```

## Monitoring:
```bash
condor_q
watch -n 30 'condor_q | tail -10'
```
