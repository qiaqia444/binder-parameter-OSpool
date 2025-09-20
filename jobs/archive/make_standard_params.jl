#!/usr/bin/env julia

using Random

# Generate parameter files for L = [8, 10, 12] study on OSPool cluster
# This focuses on the reliable ITensorCorrelators approach

function generate_standard_params()
    # System sizes that work reliably with ITensorCorrelators
    Ls = [8, 10, 12]
    
    # Lambda values for phase transition study - fine grid around critical point
    lambdas = [
        0.1, 0.2, 0.3, 0.35, 0.4,           # Pre-critical region
        0.42, 0.44, 0.46, 0.47, 0.48, 0.49, # Fine sampling around critical point  
        0.50, 0.51, 0.52, 0.53, 0.54,       # Critical region
        0.55, 0.6, 0.65, 0.7, 0.8, 0.9      # Post-critical region
    ]
    
    # Standard simulation parameters
    ntrials = 2000
    seeds = [1000 + i for i in 1:20]  # 20 different seeds per (L, λ) combination
    
    return Ls, lambdas, ntrials, seeds
end

function format_standard_line(L, lambda_x, lambda_zz, lambda, ntrials, seed, sample, out_prefix)
    return "$(L) $(lambda_x) $(lambda_zz) $(lambda) $(ntrials) $(seed) $(sample) $(out_prefix)"
end

function main()
    Ls, lambdas, ntrials, seeds = generate_standard_params()
    
    # Generate comprehensive parameter set
    println("Generating parameter files for L = $Ls")
    println("Lambda values: $lambdas")
    println("Seeds per (L,λ): $(length(seeds))")
    
    total_jobs = length(Ls) * length(lambdas) * length(seeds)
    println("Total jobs: $total_jobs")
    
    # Create parameter file
    open("params_standard_L8_10_12.txt", "w") do f
        job_count = 0
        for L in Ls
            for lambda in lambdas
                lambda_x = lambda
                lambda_zz = 1.0 - lambda
                
                for (seed_idx, seed) in enumerate(seeds)
                    job_count += 1
                    out_prefix = "std_L$(L)_lam$(lambda)_s$(seed_idx)"
                    line = format_standard_line(L, lambda_x, lambda_zz, lambda, ntrials, seed, seed_idx, out_prefix)
                    println(f, line)
                end
            end
        end
        println("Generated $job_count parameter lines")
    end
    
    # Create small test parameter file
    test_lambdas = [0.3, 0.5, 0.7]  # Just 3 lambda values for testing
    open("params_test_L8_10_12.txt", "w") do f
        job_count = 0
        test_seeds = seeds[1:3]  # Just 3 seeds for testing
        
        for L in Ls
            for lambda in test_lambdas
                lambda_x = lambda
                lambda_zz = 1.0 - lambda
                
                for (seed_idx, seed) in enumerate(test_seeds)
                    job_count += 1
                    out_prefix = "test_L$(L)_lam$(lambda)_s$(seed_idx)"
                    line = format_standard_line(L, lambda_x, lambda_zz, lambda, 50, seed, seed_idx, out_prefix)  # Reduced ntrials for test
                    println(f, line)
                end
            end
        end
        println("Generated $job_count test parameter lines")
    end
    
    println("\nParameter files created:")
    println("  params_standard_L8_10_12.txt - Full production run ($total_jobs jobs)")
    println("  params_test_L8_10_12.txt - Test run ($(length(Ls) * length(test_lambdas) * 3) jobs)")
    
    println("\nRecommended cluster deployment strategy:")
    println("  1. Test with params_test_L8_10_12.txt first")
    println("  2. If successful, deploy params_standard_L8_10_12.txt")
    println("  3. Analyze L=8,10,12 scaling behavior")
    println("  4. Use results to optimize approach for larger L")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
