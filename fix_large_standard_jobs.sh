#!/bin/bash

echo "=========================================="
echo "Fixing Large Standard Jobs Time Limit Issue"
echo "=========================================="

# Check current held jobs
HELD_COUNT=$(condor_q qia.wang -constraint 'JobStatus==5' -format "%d\n" ClusterId | wc -l)
RUNNING_COUNT=$(condor_q qia.wang -constraint 'JobStatus==2' -format "%d\n" ClusterId | wc -l)

echo "Current job status:"
echo "  Held jobs (time limit exceeded): $HELD_COUNT"
echo "  Still running: $RUNNING_COUNT"
echo ""

if [ "$HELD_COUNT" -gt 0 ]; then
    echo "Problem: Large systems (L=32, L=36) exceeded 20-hour time limit"
    echo "Solution: Remove held jobs and resubmit with 48-hour limit + more memory"
    echo ""
    
    echo "Held jobs will be removed and resubmitted with:"
    echo "  ✓ Time limit: 48 hours (was 20 hours)"  
    echo "  ✓ Memory: 12GB (was 8GB)"
    echo "  ✓ Disk: 10GB (was 8GB)"
    echo "  ✓ Optimized threading for large systems"
    echo "  ✓ Reduced chunk sizes for L≥28"
    echo ""
    
    read -p "Proceed with fixing held jobs? (y/N): " -n 1 -r
    echo ""
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Step 1: Removing held jobs..."
        condor_rm -constraint 'JobStatus==5' qia.wang
        
        echo "Step 2: Waiting for cleanup..."
        sleep 5
        
        echo "Step 3: Resubmitting with fixed parameters..."
        if [ -f "jobs/large_standard_jobs.submit" ]; then
            cd jobs
            condor_submit large_standard_jobs.submit
            cd ..
            echo ""
            echo "✅ Jobs resubmitted successfully!"
            echo ""
            echo "Monitor progress with:"
            echo "  condor_q qia.wang"
            echo ""
            echo "Expected runtimes with new limits:"
            echo "  L=20: 2-4 hours"
            echo "  L=24: 4-8 hours" 
            echo "  L=28: 8-16 hours"
            echo "  L=32: 16-32 hours"
            echo "  L=36: 24-48 hours"
        else
            echo "❌ Error: jobs/large_standard_jobs.submit not found!"
            echo "Run 'git pull origin main' first to get updated files"
        fi
    else
        echo "Operation cancelled. Held jobs remain."
        echo ""
        echo "Manual cleanup commands:"
        echo "  condor_rm -constraint 'JobStatus==5' qia.wang  # Remove held jobs"
        echo "  cd jobs && condor_submit large_standard_jobs.submit   # Resubmit"
    fi
else
    echo "✅ No held jobs found. System is healthy!"
fi

echo ""
echo "Final status:"
condor_q qia.wang