#!/bin/bash

echo "=========================================="
echo "Emergency Cleanup of Held Jobs"
echo "=========================================="

echo "WARNING: This will remove ALL held jobs for qia.wang"
echo ""

# Show what will be removed
echo "Jobs to be removed:"
condor_q qia.wang -constraint 'JobStatus==5' -format "Will remove: %d.%d\n" ClusterId ProcId

HELD_COUNT=$(condor_q qia.wang -constraint 'JobStatus==5' -format "%d\n" ClusterId | wc -l)
echo ""
echo "Total jobs to remove: $HELD_COUNT"

if [ "$HELD_COUNT" -gt 0 ]; then
    read -p "Proceed with removing $HELD_COUNT held jobs? (y/N): " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing all held jobs..."
        condor_rm -constraint 'JobStatus==5' qia.wang
        
        echo "Waiting for cleanup..."
        sleep 3
        
        echo "Done! Current status:"
        condor_q qia.wang
        
        echo ""
        echo "Next steps:"
        echo "1. git pull origin main  # Get updated job files"
        echo "2. ./submit_large_standard.sh  # Resubmit with better parameters"
    else
        echo "Operation cancelled."
    fi
else
    echo "No held jobs to remove."
fi