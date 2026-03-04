#!/bin/bash
# Startup script for ClustalDM with auto-browser opening

echo "="
echo "🧬 ClustalDM Docker Container"
echo "="
echo ""
echo "Services:"
echo "  📊 Main App: http://localhost:8000"
echo "  📁 File Browser: http://localhost:8080"
echo ""

# Start the Shiny app
shiny run app_main.py --host 0.0.0.0 --port 8000
