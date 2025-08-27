#!/bin/bash
# Download Binder simulation results from OSPool cluster

CLUSTER_USER="qia.wang"
CLUSTER_HOST="ap40.uw.osg-htc.org"
LOCAL_DIR="./binder-simulation-results"

echo "Downloading Binder parameter simulation results..."

# Create local directory
mkdir -p "$LOCAL_DIR"

# Download result files
echo "Downloading output files..."
scp -r ${CLUSTER_USER}@${CLUSTER_HOST}:~/binder-parameter-OSpool/output/ "$LOCAL_DIR/"

# Download analysis files
echo "Downloading analysis files..."
scp ${CLUSTER_USER}@${CLUSTER_HOST}:~/binder-parameter-OSpool/results_summary.csv "$LOCAL_DIR/" 2>/dev/null || true
scp ${CLUSTER_USER}@${CLUSTER_HOST}:~/binder-parameter-OSpool/*.png "$LOCAL_DIR/" 2>/dev/null || true

# Download logs
echo "Downloading log files..."
scp -r ${CLUSTER_USER}@${CLUSTER_HOST}:~/binder-parameter-OSpool/logs/ "$LOCAL_DIR/"

echo "Download complete. Results are in: $LOCAL_DIR"
result_count=$(ls $LOCAL_DIR/output/*.json 2>/dev/null | wc -l)
echo "Result files: $result_count/255"
