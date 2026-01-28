#!/bin/bash

# Submit both P=0.2 dephasing and RBIM jobs
# Run on OSPool submit node

echo "=========================================="
echo "Submitting Multiple Job Sets"
echo "=========================================="

# Create logs directory
mkdir -p logs

# Submit P=0.2 dephasing (680 jobs)
echo ""
echo "1. Submitting P=0.2 dephasing (680 jobs)..."
condor_submit dephasing_p02_jobs.submit
echo "   Status: P=0.2 jobs submitted"

# Submit RBIM (1680 jobs)
echo ""
echo "2. Submitting RBIM (1680 jobs)..."
condor_submit rbim_jobs.submit
echo "   Status: RBIM jobs submitted"

echo ""
echo "=========================================="
echo "Total jobs submitted: 2360"
echo "  - P=0.2 dephasing: 680 jobs"
echo "  - RBIM (λ_x=0): 1680 jobs"
echo "=========================================="
echo ""
echo "Monitor with:"
echo "  condor_q"
echo ""
echo "Check user summary:"
echo "  condor_q -submitter"
echo ""
