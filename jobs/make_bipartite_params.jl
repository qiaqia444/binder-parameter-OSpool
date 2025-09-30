using JSON
using Printf

function generate_bipartite_parameters()
    """
    Generate parameter sets for bipartite entropy calculations
    System sizes: L = [8, 12, 16]
    Parameter relationship: lambda_x = lambda, lambda_zz = 1 - lambda
    """
    
    # System sizes to explore
    L_values = [8, 12, 16]
    
    # Lambda parameter from 0 to 1 with complementary relationship
    # lambda_x = lambda, lambda_zz = 1 - lambda
    lambda_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    parameters = []
    job_id = 1
    
    println("Generating bipartite entropy parameter sets...")
    
    for L in L_values
        for lambda in lambda_values
            lambda_x = round(lambda, digits=1)
            lambda_zz = round(1.0 - lambda, digits=1)
            
            param_set = Dict(
                "job_id" => job_id,
                "L" => L,
                "lambda" => lambda,
                "lambda_x" => lambda_x,
                "lambda_zz" => lambda_zz,
                "T_max" => 2 * L,  # Evolution time scales with system size
                "ntrials" => 500   # Fixed number of trials per job
            )
            
            push!(parameters, param_set)
            job_id += 1
        end
    end
    
    println("Generated $(length(parameters)) parameter sets")
    println("System sizes: $L_values")
    println("Lambda values: $lambda_values") 
    println("Parameter relationship: λ_x = λ, λ_zz = 1 - λ")
    
    # Save parameter sets as JSON
    output_file = "jobs/bipartite_params.json"
    mkpath("jobs")
    
    open(output_file, "w") do f
        JSON.print(f, parameters, 2)
    end
    
    println("Parameters saved to: $output_file")
    
    # Create HTCondor-compatible parameter file
    params_file = "jobs/params_bipartite.txt"
    open(params_file, "w") do f
        for param in parameters
            # Format: L lambda_x lambda_zz lambda ntrials seed sample out_prefix
            lambda_str = @sprintf("%.1f", param["lambda"])
            out_prefix = "bipartite_L$(param["L"])_lambda$(lambda_str)"
            println(f, "$(param["L"]) $(param["lambda_x"]) $(param["lambda_zz"]) $(param["lambda"]) $(param["ntrials"]) $(param["job_id"] + 1000) $(param["job_id"]) $(out_prefix)")
        end
    end
    
    println("HTCondor parameters saved to: $params_file")
    
    # Also create a summary file
    summary = Dict(
        "total_jobs" => length(parameters),
        "L_values" => L_values,
        "lambda_values" => lambda_values,
        "parameter_relationship" => "lambda_x = lambda, lambda_zz = 1 - lambda",
        "trials_per_job" => 500,
        "description" => "Bipartite entanglement entropy with complementary weak measurement strengths"
    )
    
    summary_file = "jobs/bipartite_summary.json"
    open(summary_file, "w") do f
        JSON.print(f, summary, 2)
    end
    
    println("Summary saved to: $summary_file")
    
    return parameters
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_bipartite_parameters()
end