#!/bin/bash

# Set up test directories and paths
RUN_DIR="test_run_dir"
BIN_DIR="mock_testRun/bin"
CONFIG_DIR="configFiles"
SENS_CSV="inputs/sensitivity_analysis.csv"
EXPERIMENTAL_CSV="inputs/experimental_data.csv"

# Create the run directory if it doesn't exist
mkdir -p "$RUN_DIR"
if [ -f venv/bin/activate ]; then
    source venv/bin/activate
elif [ -f venv/Scripts/activate ]; then
    source venv/Scripts/activate
else
    echo "No virtual environment found!"
    exit 1
fi

# Run the rope script in test mode
python -m src.calibration.run_rope \
    --log-level DEBUG \
    --param-ranking random_forest \
    --param-num 3 \
    --num-iterations 500 \
    --run-dir-parent "$RUN_DIR" \
    --sensitivity-analysis-csv "$SENS_CSV" \
    --experimental-data-csv "$EXPERIMENTAL_CSV" \
    --config-file-dir "$CONFIG_DIR" \
    --bin-dir "$BIN_DIR" \
    --parallel "mpc" \
    2>&1 | tee "rope_test_output.log"
