#!/bin/bash
#
# Launch N workers of generate_lyapunov.py simultaneously
#
# Usage:
#   ./scripts/run_lyapunov_workers.sh N [options]
#
# Examples:
#   ./scripts/run_lyapunov_workers.sh 4                    # Launch 4 workers with defaults
#   ./scripts/run_lyapunov_workers.sh 8 --test             # Launch 8 workers in test mode
#   ./scripts/run_lyapunov_workers.sh 4 --stop-time 10.0   # Custom stop time
#

set -e

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 N [additional_args...]"
    echo "  N: Number of workers to launch"
    echo "  additional_args: Any additional arguments to pass to generate_lyapunov.py"
    echo ""
    echo "Examples:"
    echo "  $0 4                     # Launch 4 workers"
    echo "  $0 8 --test              # Launch 8 workers in test mode"
    echo "  $0 4 --stop-time 10.0    # Launch 4 workers with custom stop time"
    exit 1
fi

N=$1
shift  # Remove N from arguments, remaining args passed to python script

# Validate N is a positive integer
if ! [[ "$N" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: N must be a positive integer, got '$N'"
    exit 1
fi

# Get the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Activate conda environment
eval "$(/home/krosenfeld/miniforge3/bin/conda shell.bash hook)"
conda activate dedalus3

echo "Launching $N workers for generate_lyapunov.py..."
echo "Additional arguments: $@"
echo "Working directory: $PROJECT_DIR"
echo ""

# Array to store PIDs
declare -a PIDS

# Launch N workers
for ((i=0; i<N; i++)); do
    echo "Starting worker $i of $N..."
    python "$SCRIPT_DIR/generate_lyapunov.py" -ntot "$N" -tid "$i" "$@" &
    PIDS+=($!)
done

echo ""
echo "All $N workers launched. PIDs: ${PIDS[*]}"
echo "Waiting for all workers to complete..."

# Wait for all background processes
FAILED=0
for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then
        echo "Worker with PID $pid failed"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
if [ $FAILED -eq 0 ]; then
    echo "All $N workers completed successfully!"
else
    echo "$FAILED of $N workers failed"
    exit 1
fi
