#!/bin/bash

# Submit script for critical scaling analysis jobs
# L=24,28,32,36 at lambda=0.5 with enhanced resources

echo "=== Critical Scaling Analysis Job Submission ==="
echo "Submitting jobs for L=24,28,32,36 at lambda=0.5 (critical point)"
echo "Enhanced resources: 72h runtime, 16GB memory, 20GB disk"

# Change to jobs directory
cd "$(dirname "$0")/jobs"

# Check if parameter file exists
if [[ ! -f "params_scaling_critical.txt" ]]; then
    echo "ERROR: params_scaling_critical.txt not found!"
    echo "Please run this script from the project root directory"
    exit 1
fi

# Show parameter file summary
echo ""
echo "Parameter file summary:"
echo "Number of jobs: $(wc -l < params_scaling_critical.txt)"
echo "System sizes: L=24,28,32,36"
echo "Lambda values: 0.5 (critical point)"
echo "Samples per size: 3"
echo "Total jobs: 12"

echo ""
echo "First few parameter combinations:"
head -5 params_scaling_critical.txt

# Create logs directory if it doesn't exist
mkdir -p logs

# Submit the jobs
echo ""
echo "Submitting critical scaling jobs to HTCondor..."
condor_submit scaling_critical_jobs.submit

if [[ $? -eq 0 ]]; then
    echo ""
    echo "✅ Jobs submitted successfully!"
    echo ""
    echo "Monitor job status with:"
    echo "  condor_q"
    echo "  condor_q -nobatch"
    echo "  watch -n 30 condor_q"
    echo ""
    echo "Check job details:"
    echo "  condor_q -l [job_id]"
    echo ""
    echo "Remove jobs if needed:"
    echo "  condor_rm [job_id]  # Remove specific job"
    echo "  condor_rm \$USER     # Remove all your jobs"
    echo ""
    echo "Expected completion time: 24-72 hours for largest systems"
    echo "Results will be in output/ directory"
else
    echo "❌ Job submission failed!"
    exit 1
fi