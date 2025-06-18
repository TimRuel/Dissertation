#!/bin/bash

# Usage: ./render_report.sh path/to/config_snapshot.yml
# Example: ./render_report.sh experiments/exp_v1.0.0/simulations/sim_01/iter_001/config_snapshot.yml

set -e

# Check for input argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 path/to/config_snapshot.yml"
  exit 1
fi

INPUT_PATH="$1"

# Resolve absolute path to config_snapshot.yml
if command -v realpath >/dev/null 2>&1; then
  CONFIG_PATH=$(realpath "$INPUT_PATH")
else
  CONFIG_PATH=$(cd "$(dirname "$INPUT_PATH")" && pwd)/$(basename "$INPUT_PATH")
fi

# Extract parent directory (run directory)
RUN_DIR=$(dirname "$CONFIG_PATH")

# Define a custom output filename
OUTPUT_FILE="run_report.html"

# Render the report
quarto render run_report.qmd \
  --execute-params "$CONFIG_PATH" \
  --output "$OUTPUT_FILE"

mv "run_report.html" "$RUN_DIR/"

echo "Report saved to: $RUN_DIR/$OUTPUT_FILE"