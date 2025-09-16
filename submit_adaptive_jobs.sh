#!/bin/bash

echo "=== Edwards-Anderson Binder Parameter Cluster Submission ==="
echo "System sizes: L = 8, 12, 16"
echo "Lambda values: 17 total (coarse + fine)"
echo "Trials per job: 2000 (production), 200 (test)"
echo "T_max: 2L (automatic)"
echo ""

# Function to submit jobs and show status
submit_and_monitor() {
    local submit_file=$1
    local job_type=$2
    local param_file=$3
    
    echo "=== Submitting $job_type jobs ==="
    echo "Parameter file: $param_file"
    echo "Job count: $(wc -l < $param_file)"
    echo ""
    
    # Submit jobs
    echo "Submitting jobs..."
    condor_submit $submit_file
    
    echo ""
    echo "Jobs submitted! Use these commands to monitor:"
    echo "  condor_q -nobatch                    # Check job status"
    echo "  condor_q -analyze                    # Check why jobs aren't running"
    echo "  watch -n 30 'condor_q -nobatch'     # Auto-refresh job status"
    echo ""
}

# Check if we're on the submit machine
if ! command -v condor_submit &> /dev/null; then
    echo "ERROR: condor_submit not found. Are you on the submit machine?"
    echo "Please run: ssh qia.wang@ap40.uw.osg-htc.org"
    exit 1
fi

# Check if parameter files exist
if [ ! -f "params_adaptive_test.txt" ]; then
    echo "ERROR: params_adaptive_test.txt not found"
    echo "Please run the parameter generation script first"
    exit 1
fi

if [ ! -f "params_adaptive.txt" ]; then
    echo "ERROR: params_adaptive.txt not found"
    echo "Please run the parameter generation script first"
    exit 1
fi

# Ask user what they want to submit
echo "What would you like to submit?"
echo "1) Test jobs first (9 jobs, 200 trials each) - RECOMMENDED"
echo "2) Production jobs (153 jobs, 2000 trials each)"
echo "3) Both (test first, then production)"
echo ""
read -p "Enter choice (1-3): " choice

case $choice in
    1)
        submit_and_monitor "jobs_adaptive_test.submit" "TEST" "params_adaptive_test.txt"
        ;;
    2)
        echo "WARNING: This will submit 153 production jobs with 2000 trials each."
        read -p "Are you sure? (y/N): " confirm
        if [[ $confirm =~ ^[Yy]$ ]]; then
            submit_and_monitor "jobs_adaptive_production.submit" "PRODUCTION" "params_adaptive.txt"
        else
            echo "Cancelled."
        fi
        ;;
    3)
        submit_and_monitor "jobs_adaptive_test.submit" "TEST" "params_adaptive_test.txt"
        echo ""
        echo "Test jobs submitted. Monitor them first, then run this script again to submit production jobs."
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo ""
echo "=== Additional useful commands ==="
echo "Monitor job progress:"
echo "  condor_q -nobatch"
echo "  condor_q -run"
echo "  condor_q -hold"
echo ""
echo "Check specific job output:"
echo "  ls -la logs/"
echo "  tail -f logs/[cluster].[process].out"
echo ""
echo "When jobs complete, check results:"
echo "  ls -la output/"
echo "  wc -l output/*.json"
echo ""
echo "Cancel all jobs if needed:"
echo "  condor_rm \$(whoami)"
