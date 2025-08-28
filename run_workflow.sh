#!/usr/bin/env bash
set -e

echo "=== Binder Parameter Simulation Workflow ==="
echo

# Check if we're in the right directory
if [ ! -f "Project.toml" ] || [ ! -d "src" ]; then
    echo "Error: Please run this script from the project root directory"
    exit 1
fi

# Set Julia path
JULIA_PATH="julia"
if [ ! -f "$JULIA_PATH" ]; then
    echo "Error: Julia not found at $JULIA_PATH"
    echo "Please update JULIA_PATH in this script to point to your Julia installation"
    exit 1
fi

# Function to show usage
show_usage() {
    echo "Usage: $0 [COMMAND]"
    echo
    echo "Commands:"
    echo "  test        - Test the simulation locally with a small example"
    echo "  params      - Generate parameter files for cluster submission"
    echo "  container   - Build the Singularity container (requires singularity/apptainer)"
    echo "  submit      - Submit jobs to HTCondor cluster"
    echo "  status      - Check job status"
    echo "  collect     - Collect and analyze results after jobs complete"
    echo "  clean       - Clean up generated files"
    echo "  all         - Run test, params, and show submission instructions"
    echo
    echo "Examples:"
    echo "  $0 test                    # Test locally"
    echo "  $0 params                  # Generate parameters with defaults (L=[12,16,20,24,28], λ=0.1-0.9, 1000 trials)"
    echo "  $0 params --Ls 12,16,20 --lmin 0.2 --lmax 0.8 --samples 5 --ntrials 500"
    echo "  $0 submit                  # Submit to cluster"
}

# Parse command
COMMAND=${1:-help}

case $COMMAND in
    "params")
        shift  # Remove 'params' from arguments
        echo "Generating parameter files..."
        mkdir -p jobs logs output
        $JULIA_PATH jobs/make_params.jl "$@"
        
        # Show what was generated
        if [ -f "jobs/params.txt" ]; then
            NUM_JOBS=$(wc -l < jobs/params.txt)
            echo "Generated $NUM_JOBS jobs in jobs/params.txt"
            echo "First few parameter sets:"
            head -3 jobs/params.txt | jq .
        fi
        ;;
        
    "container")
        echo "Building Singularity container..."
        if ! command -v singularity &> /dev/null && ! command -v apptainer &> /dev/null; then
            echo "Error: Neither singularity nor apptainer found"
            echo "Please install one of them to build containers"
            exit 1
        fi
        
        BUILDER="singularity"
        if command -v apptainer &> /dev/null; then
            BUILDER="apptainer"
        fi
        
        mkdir -p containers
        $BUILDER build containers/container.sif containers/image.def
        echo "Container built: containers/container.sif"
        ;;
        
    "submit")
        echo "Submitting jobs to HTCondor..."
        
        # Check if parameter file exists
        if [ ! -f "jobs/params.txt" ]; then
            echo "Error: Parameter file jobs/params.txt not found"
            echo "Run '$0 params' first to generate parameters"
            exit 1
        fi
        
        # Check if container exists
        if [ ! -f "containers/container.sif" ]; then
            echo "Warning: Container file containers/container.sif not found"
            echo "You may need to build it with '$0 container' or ensure it's available on the cluster"
        fi
        
        # Create necessary directories
        mkdir -p logs output
        
        # Submit jobs
        condor_submit jobs/jobs.submit
        ;;
        
    "status")
        echo "Checking job status..."
        condor_q
        ;;
        
    "collect")
        echo "Collecting and analyzing results..."
        $JULIA_PATH analyze_results.jl
        ;;
        
    "clean")
        echo "Cleaning up generated files..."
        rm -rf logs/* output/* jobs/params.txt
        echo "Cleaned logs, output, and parameter files"
        ;;
        
    "all")
        echo "Running complete workflow (parameter generation)..."
        echo
        echo "Generating parameter files..."
        mkdir -p jobs logs output
        $JULIA_PATH jobs/make_params.jl
        echo
        echo "Setup complete! Next steps:"
        echo "1. Build container (if needed): $0 container"
        echo "2. Edit jobs/jobs.submit to set your project name"
        echo "3. Submit jobs: $0 submit"
        echo "4. Monitor jobs: $0 status"
        echo "5. Collect results: $0 collect"
        ;;
        
    "help"|*)
        show_usage
        ;;
esac
