#!/usr/bin/env bash
set -euo pipefail

# Find Singularity executable
if command -v singularity >/dev/null 2>&1; then
    SINGULARITY_CMD="singularity"
elif command -v apptainer >/dev/null 2>&1; then
    SINGULARITY_CMD="apptainer"
elif [ -f "/usr/bin/singularity" ]; then
    SINGULARITY_CMD="/usr/bin/singularity"
elif [ -f "/usr/local/bin/singularity" ]; then
    SINGULARITY_CMD="/usr/local/bin/singularity"
else
    echo "ERROR: Neither singularity nor apptainer found!"
    echo "Available commands:"
    which -a singularity apptainer || echo "None found"
    echo "PATH: $PATH"
    exit 1
fi

echo "Using container runtime: $SINGULARITY_CMD"

# Use transferred container file (transferred as part of job)
CONTAINER_PATH="./container.sif"

# Check if container file exists
if [ ! -f "$CONTAINER_PATH" ]; then
    echo "ERROR: Container file not found at $CONTAINER_PATH"
    echo "Current directory contents:"
    ls -la
    exit 1
fi

echo "Using container: $CONTAINER_PATH"

# Create project directory structure if needed
mkdir -p output
mkdir -p .julia_local

# Set up environment for writable Julia depot
export JULIA_DEPOT_PATH="$(pwd)/.julia_local:/tmp/julia"

# Try different Singularity execution methods
echo "Attempting to run container..."

# Method 1: Try standard execution with bind mounts
if $SINGULARITY_CMD exec --bind $(pwd):/work --pwd /work $CONTAINER_PATH julia --version >/dev/null 2>&1; then
    echo "Using standard bind mount mode"
    $SINGULARITY_CMD exec --bind $(pwd):/work --pwd /work $CONTAINER_PATH julia --project=. run.jl "$@"
# Method 2: Try with --fakeroot
elif $SINGULARITY_CMD exec --fakeroot --bind $(pwd):/work --pwd /work $CONTAINER_PATH julia --version >/dev/null 2>&1; then
    echo "Using fakeroot mode"
    $SINGULARITY_CMD exec --fakeroot --bind $(pwd):/work --pwd /work $CONTAINER_PATH julia --project=. run.jl "$@"
# Method 3: Try without bind mounts (copy files)
else
    echo "Trying without bind mounts - copying files to temp directory"
    ORIGINAL_DIR=$(pwd)
    WORK_DIR="/tmp/job_work_$$"
    mkdir -p $WORK_DIR
    
    # Copy essential files only (avoid copying directories that might cause issues)
    cp *.toml *.jl *.sif $WORK_DIR/ 2>/dev/null || true
    cp -r src $WORK_DIR/ 2>/dev/null || true
    
    cd $WORK_DIR
    mkdir -p output
    
    echo "Running in temporary directory: $WORK_DIR"
    $SINGULARITY_CMD exec ./container.sif julia --project=. run.jl "$@"
    
    # Copy results back
    echo "Copying results back to $ORIGINAL_DIR"
    cp -r output/* $ORIGINAL_DIR/output/ 2>/dev/null || true
fi

echo "Container execution completed"
