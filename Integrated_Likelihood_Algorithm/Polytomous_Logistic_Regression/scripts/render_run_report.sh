#!/bin/bash

# Usage: ./render_report.sh path/to/run_dir
# Example: ./render_report.sh experiments/exp_v1.0.0/simulations/sim_01/iter_001

set -e

# Check for input argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 path/to/run_dir (e.g. experiments/exp_v1.0.0/simulations/sim_01/iter_001)"
  exit 1
fi

# Resolve absolute path to run_dir
INPUT_PATH="$1"
if command -v realpath >/dev/null 2>&1; then
  RUN_DIR=$(realpath "$INPUT_PATH")
else
  RUN_DIR=$(cd "$INPUT_PATH" && pwd)
fi

# Paths
OUTPUT_FILE="run_report.html"
REPORT_PATH="$RUN_DIR/$OUTPUT_FILE"

# Render with quarto
quarto render run_report.qmd \
  -P run_dir="$RUN_DIR" 

mv "run_report.html" "$RUN_DIR/"

echo "Report saved to: $REPORT_PATH"
