#!/bin/bash
# Launch ClustalDM Shiny App

echo "🧬 Starting ClustalDM Interactive App..."
echo "Make sure you're in the 'structural' conda environment"
echo ""

# Check if in conda environment
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "⚠️  No conda environment activated!"
    echo "Please run: conda activate structural"
    echo "Then try again: ./launch_app.sh"
    exit 1
fi

echo "Environment: $CONDA_DEFAULT_ENV"
echo "Starting Shiny server..."
echo "Opening browser automatically..."
echo ""

shiny run --reload --launch-browser app_main.py
