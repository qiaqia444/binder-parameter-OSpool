#!/bin/bash

# Submit corrected MIPT jobs for L = [8, 12, 16] with fixed weak measurements
echo "Submitting corrected MIPT jobs for crossing point analysis..."
echo "System sizes: L = 8, 12, 16"
echo "Total jobs: 153"
echo "Fixed weak measurement operators to match sparse matrix implementation"

cd jobs
condor_submit jobs.submit -queue params from params_production.txt

echo "Jobs submitted!"
echo ""
echo "Monitor with:"
echo "  condor_q"
echo "  condor_q -better-analyze"
echo "  watch -n 30 'condor_q | tail -10'"
echo ""
echo "When complete, transfer results with:"
echo "  scp -r output/ your_local_machine:"
