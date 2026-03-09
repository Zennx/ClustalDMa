#!/bin/bash
# ClustalDMα CLI Wrapper Script
# Makes it easier to run the CLI with the correct Python environment

# Detect script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Activate conda environment (modify if your environment has a different name)
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate structural 2>/dev/null || conda activate base
fi

# Run the CLI script with all passed arguments
python "$SCRIPT_DIR/clustal_cli.py" "$@"
