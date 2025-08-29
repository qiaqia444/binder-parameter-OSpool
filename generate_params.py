#!/usr/bin/env python3

# Generate all 255 parameter combinations for production run
# Based on README specifications

# System sizes from README
Ls = [12, 16, 20, 24, 28]

# Lambda values from README
coarse_lambdas = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
fine_lambdas = [0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54]
lambdas = sorted(list(set(coarse_lambdas + fine_lambdas)))

# Parameters
samples_per_lambda = 3
maxdim = 256
cutoff = 1e-12
ntrials = 1000
chunk4 = 50000
base_seed = 1234
outdir = "output"

lines = []
job_id = 0

for L in Ls:
    for lam in lambdas:
        for sample in range(1, samples_per_lambda + 1):
            job_id += 1
            lambda_x = round(lam, 6)
            lambda_zz = round(1.0 - lam, 6)
            seed = base_seed + job_id
            out_prefix = f"L{L}_lam{lam:.3f}_s{sample}"
            
            line = f"--L {L} --lambda_x {lambda_x} --lambda_zz {lambda_zz} --lambda {lam} --maxdim {maxdim} --cutoff {cutoff} --ntrials {ntrials} --chunk4 {chunk4} --seed {seed} --sample {sample} --out_prefix {out_prefix} --outdir {outdir}"
            lines.append(line)

# Write to file
with open("jobs/params_production.txt", "w") as f:
    for line in lines:
        f.write(line + "\n")

print(f"Generated {len(lines)} parameter combinations")
print(f"System sizes: {len(Ls)}")
print(f"Lambda values: {len(lambdas)}")
print(f"Samples per lambda: {samples_per_lambda}")
print(f"Total: {len(Ls)} × {len(lambdas)} × {samples_per_lambda} = {len(lines)} jobs")
