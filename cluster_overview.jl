#!/usr/bin/env julia
"""
Quick overview of the simulation setup
"""

println("Binder Parameter Simulation Setup")

# Lambda values
lambda_coarse = 0.1:0.1:0.9
lambda_fine = 0.4:0.02:0.6
all_lambdas = sort(unique(vcat(collect(lambda_coarse), collect(lambda_fine))))

println("Simulation Parameters:")
println("  System sizes (L): [12, 16, 20, 24, 28]")
println("  Lambda values ($(length(all_lambdas)) points): $all_lambdas")
println("  Samples per (L,λ): 3")
println("  Monte Carlo trials per job: 1000")
println("  Time evolution: T_max = 2L")
println()

total_jobs = 5 * length(all_lambdas) * 3
total_trials = total_jobs * 1000

println("💻 Computational Scale:")
println("  Total jobs: $total_jobs")
println("  Total Monte Carlo trials: $(total_trials)")
println("  Memory per job: 8 GB")
println("  Estimated runtime per job: 1-4 hours (depends on L)")
println()

println("🔬 Physics:")
println("  Weak X measurements: λₓ = λ")
println("  Weak ZZ measurements: λ_zz = 1-λ")
println("  Edwards-Anderson Binder parameter: B = 1 - S₄/(3S₂²)")
println()

println("Files Ready for Cluster:")
println("  ✓ src/BinderSim.jl - Simulation module")
println("  ✓ run.jl - Main execution script")
println("  ✓ jobs/jobs.submit - HTCondor submission file")
println("  ✓ jobs/params.txt - Parameter sets ($(total_jobs) jobs)")
println("  ✓ containers/image.def - Singularity container")
println()

println("Next Steps on Cluster:")
println("  1. git clone https://github.com/qiaqia444/binder-parameter-OSpool.git")
println("  2. cd binder-parameter-OSpool")
println("  3. Edit jobs/jobs.submit (set your project name)")
println("  4. ./run_workflow.sh container  # Build container")
println("  5. ./run_workflow.sh submit     # Submit jobs")
println("  6. ./run_workflow.sh status     # Monitor progress")
println("  7. ./run_workflow.sh collect    # Analyze results")
println()

# Estimate total computational time
println("Computational Estimate:")
println("  Conservative: ~$(round(total_trials * 0.001 / 3600, digits=1)) CPU-hours")
println("  Optimistic: ~$(round(total_trials * 0.0005 / 3600, digits=1)) CPU-hours")
println("  With $(total_jobs) parallel jobs: 1-4 hours wall-clock time")
