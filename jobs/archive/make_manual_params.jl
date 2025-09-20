#!/usr/bin/env julia

using Random

# Generate parameter files for L = [8, 10, 12] study using MANUAL correlators
# This creates the same parameter set but for manual correlator testing

function generate_manual_params()
    # Same system sizes as standard run
    Ls = [8, 10, 12]
    
    # Same lambda values for direct comparison
    lambdas = [
        0.1, 0.2, 0.3, 0.35, 0.4,           # Pre-critical region
        0.42, 0.44, 0.46, 0.47, 0.48, 0.49, # Fine sampling around critical point  
        0.50, 0.51, 0.52, 0.53, 0.54,       # Critical region
        0.55, 0.6, 0.65, 0.7, 0.8, 0.9      # Post-critical region
    ]
    
    # Reduced trials for manual correlators (computationally expensive)
    ntrials = 500  # Reduced from 2000 for feasibility
    seeds = [1000 + i for i in 1:20]  # Same 20 seeds for comparison
    
    return Ls, lambdas, ntrials, seeds
end

function format_manual_line(L, lambda_x, lambda_zz, lambda, ntrials, seed, sample, out_prefix)
    return "$(L) $(lambda_x) $(lambda_zz) $(lambda) $(ntrials) $(seed) $(sample) $(out_prefix)"
end

function main()
    Ls, lambdas, ntrials, seeds = generate_manual_params()
    
    # Generate manual correlator parameter set
    println("Generating parameter files for MANUAL correlators L = $Ls")
    println("Lambda values: $lambdas")
    println("Seeds per (L,λ): $(length(seeds))")
    println("Trials per job: $ntrials (reduced for manual correlators)")
    
    total_jobs = length(Ls) * length(lambdas) * length(seeds)
    println("Total jobs: $total_jobs")
    
    # Create parameter file for manual correlators
    open("jobs/params_manual_L8_10_12.txt", "w") do f
        job_count = 0
        for L in Ls
            for lambda in lambdas
                lambda_x = lambda
                lambda_zz = 1.0 - lambda
                
                for (seed_idx, seed) in enumerate(seeds)
                    job_count += 1
                    out_prefix = "manual_L$(L)_lam$(lambda)_s$(seed_idx)"
                    line = format_manual_line(L, lambda_x, lambda_zz, lambda, ntrials, seed, seed_idx, out_prefix)
                    println(f, line)
                end
            end
        end
        println("Generated $job_count parameter lines")
    end
    
    # Create small test parameter file
    test_lambdas = [0.3, 0.5, 0.7]  # Just 3 lambda values for testing
    open("params_manual_test_L8_10_12.txt", "w") do f
        job_count = 0
        test_seeds = seeds[1:3]  # Just 3 seeds for testing
        
        for L in Ls
            for lambda in test_lambdas
                lambda_x = lambda
                lambda_zz = 1.0 - lambda
                
                for (seed_idx, seed) in enumerate(test_seeds)
                    job_count += 1
                    out_prefix = "manual_test_L$(L)_lam$(lambda)_s$(seed_idx)"
                    line = format_manual_line(L, lambda_x, lambda_zz, lambda, 10, seed, seed_idx, out_prefix)  # Very reduced for test
                    println(f, line)
                end
            end
        end
        println("Generated $job_count test parameter lines")
    end
    
    println("\nManual correlator parameter files created:")
    println("  params_manual_L8_10_12.txt - Manual correlators run ($total_jobs jobs, $ntrials trials each)")
    println("  params_manual_test_L8_10_12.txt - Test run ($(length(Ls) * length(test_lambdas) * 3) jobs)")
    
    println("\nComparison with ITensorCorrelators:")
    println("  Standard (ITensorCorrelators): 2000 trials per job")
    println("  Manual (tensor contractions): $ntrials trials per job")
    println("  This allows direct method comparison on same parameter space")
    
    println("\nRecommended cluster deployment strategy:")
    println("  1. Test with params_manual_test_L8_10_12.txt first")
    println("  2. If successful, deploy params_manual_L8_10_12.txt")
    println("  3. Compare results with standard ITensorCorrelators data")
    println("  4. Analyze method differences and computational efficiency")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
