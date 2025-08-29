#!/usr/bin/env julia

# Script to generate HTCondor submit file with command line arguments
# This avoids all JSON parsing issues

function generate_submit_file()
    # Parameter ranges
    L_values = [12]
    lambda_x_values = [0.1]
    lambda_zz_values = [0.9]
    lambda_values = [0.1]
    maxdim_values = [256]
    cutoff_values = [1e-12]
    ntrials_values = [1000]
    chunk4_values = [50000]
    seed_base = 1235
    samples = 1:3
    
    submit_content = """universe                = vanilla
executable              = jobs/run.sh
log                     = logs/\$(Cluster).log
output                  = logs/\$(Cluster).\$(Process).out
error                   = logs/\$(Cluster).\$(Process).err
request_cpus            = 4
request_memory          = 4 GB
request_disk            = 1 GB
+JobDurationCategory    = "Medium"
#+ProjectName            = "qia.wang"
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = src/, run.jl, Project.toml, Manifest.toml
requirements            = (HAS_SINGULARITY =?= TRUE)

"""
    
    job_count = 0
    for L in L_values
        for lambda_x in lambda_x_values
            for lambda_zz in lambda_zz_values
                for lambda in lambda_values
                    for maxdim in maxdim_values
                        for cutoff in cutoff_values
                            for ntrials in ntrials_values
                                for chunk4 in chunk4_values
                                    for sample in samples
                                        seed = seed_base + sample - 1
                                        out_prefix = "L$(L)_lam$(lambda)_s$(sample)"
                                        
                                        args = "--L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --maxdim $maxdim --cutoff $cutoff --ntrials $ntrials --chunk4 $chunk4 --seed $seed --sample $sample --out_prefix $out_prefix --outdir output"
                                        
                                        submit_content *= "arguments = $args\n"
                                        submit_content *= "queue\n\n"
                                        job_count += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    println("Generated $job_count jobs")
    return submit_content
end

# Generate and write the submit file
content = generate_submit_file()
open("jobs/jobs.submit", "w") do io
    write(io, content)
end

println("Generated jobs/jobs.submit with command line arguments")
