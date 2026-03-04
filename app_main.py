#!/usr/bin/env python3
"""
ClustalDM - Main Interactive Application
Strategy: Compute both RMSD + Jaccard, cluster by Jaccard only, use RMSD as QC
"""
import webbrowser

import os
import sys
from pathlib import Path
import warnings
import tempfile
import shutil
import zipfile
from typing import List

# Suppress joblib resource tracker warnings (harmless cleanup warnings from sklearn)
warnings.filterwarnings('ignore', category=UserWarning, module='joblib')

print("Loading dependencies...")
try:
    from shiny import App, ui, render, reactive
    import pandas as pd
    import numpy as np
    print("✓ Shiny and data libraries loaded")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)

# Import modular components
try:
    from core import PDBClusterer, find_pdb_files, InterfaceAnalyzer
    from core.hdock_parser_docker import HDOCKParser
    from visualization import (
        create_scatter_multimethod,
        create_hotspot_histogram,
        create_molstar_viewer_html as create_mol_viewer_html,  # Use Molstar instead of py3Dmol
        create_distance_heatmap,
        create_cluster_size_distribution,
    )
    print("✓ Core and visualization modules loaded (using Molstar)")
except ImportError as e:
    print(f"✗ Module import failed: {e}")
    sys.exit(1)

print("\nAll imports successful! Starting app...\n")

# Global state
clust = None
pdb_dir_path = None


# UI Definition
app_ui = ui.page_fluid(
    # Add Plotly.js
    ui.head_content(
        ui.HTML('<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'),
        ui.tags.style("""
            /* Dark mode theme */
            body {
                background-color: #1a1a1a !important;
                color: #e0e0e0 !important;
            }
            
            * {
                color: #e0e0e0;
            }
            
            .container-fluid {
                background-color: #1a1a1a !important;
                color: #e0e0e0 !important;
            }
            
            .well, .card {
                background-color: #2d2d2d !important;
                border: 1px solid #404040 !important;
                color: #e0e0e0 !important;
            }
            
            .card-header {
                background-color: #353535 !important;
                border-bottom: 1px solid #404040 !important;
                color: #ffffff !important;
            }
            
            .sidebar {
                background-color: #252525 !important;
                border-right: 1px solid #404040 !important;
                color: #e0e0e0 !important;
            }
            
            h1, h2, h3, h4, h5, h6, p, span, div, label {
                color: #e0e0e0 !important;
            }
            
            .nav-tabs {
                border-bottom: 1px solid #404040 !important;
            }
            
            .nav-tabs .nav-link {
                background-color: #2d2d2d !important;
                border: 1px solid #404040 !important;
                color: #b0b0b0 !important;
            }
            
            .nav-tabs .nav-link:hover {
                background-color: #3d3d3d !important;
                color: #ffffff !important;
            }
            
            .nav-tabs .nav-link.active {
                background-color: #1a1a1a !important;
                border-bottom-color: transparent !important;
                color: #ffffff !important;
            }
            
            input[type="text"], input[type="number"], select, textarea {
                background-color: #2d2d2d !important;
                border: 1px solid #404040 !important;
                color: #e0e0e0 !important;
            }
            
            input[type="text"]:focus, input[type="number"]:focus, select:focus {
                background-color: #353535;
                border-color: #5a9fd4;
                color: #ffffff;
            }
            
            .btn-primary {
                background-color: #5a9fd4;
                border-color: #4a8fc4;
            }
            
            .btn-primary:hover {
                background-color: #4a8fc4;
                border-color: #3a7fb4;
            }
            
            .btn-outline-secondary {
                background-color: #2d2d2d;
                border-color: #404040;
                color: #b0b0b0;
            }
            
            .btn-outline-secondary:hover {
                background-color: #3d3d3d;
                color: #ffffff;
            }
            
            .shiny-notification {
                background-color: #2d2d2d !important;
                border: 1px solid #404040 !important;
                color: #e0e0e0 !important;
            }
            
            pre {
                background-color: #1a1a1a !important;
                border: 1px solid #404040 !important;
                color: #00ff00 !important;
            }
            
            /* DataGrid/Table styling */
            .shiny-data-grid {
                background-color: #2d2d2d !important;
            }
            
            .shiny-data-grid th {
                background-color: #4a5568 !important;
                color: #ffffff !important;
                border-color: #404040 !important;
            }
            
            .shiny-data-grid td {
                background-color: #2d2d2d !important;
                color: #e0e0e0 !important;
                border-color: #404040 !important;
            }
            
            .shiny-data-grid tr:hover {
                background-color: #3d3d3d !important;
            }
            
            /* Plotly with transparent backgrounds + default template - CSS controls bg/text/gridlines */
            .plotly-graph-div,
            .js-plotly-plot {
                background-color: #1a1a1a;
            }
            
            /* Dark mode: make text light */
            .js-plotly-plot text {
                fill: #e0e0e0 !important;
            }
            
            /* Dark mode: hover labels */
            .js-plotly-plot .hoverlayer .hovertext,
            .js-plotly-plot .hoverlayer text {
                fill: #212529 !important;
            }
            
            .js-plotly-plot .hoverlayer path {
                fill: #ffffff !important;
                stroke: #cccccc !important;
            }
            
            /* Dark mode: make gridlines visible */
            .js-plotly-plot .gridlayer line,
            .js-plotly-plot .gridlayer path {
                stroke: #404040 !important;
            }
            
            /* Dark mode: zero lines AND axis borders */
            .js-plotly-plot .zerolinelayer line,
            .js-plotly-plot line.xlines-above,
            .js-plotly-plot line.ylines-above,
            .js-plotly-plot .xlines-below line,
            .js-plotly-plot .ylines-below line,
            .js-plotly-plot .overaxes-above line {
                stroke: #808080 !important;
            }
            
            /* Dark mode: catch all for axis elements */
            .js-plotly-plot g.xaxislayer-above path,
            .js-plotly-plot g.yaxislayer-above path,
            .js-plotly-plot g.xaxislayer-above line,
            .js-plotly-plot g.yaxislayer-above line {
                stroke: #808080 !important;
            }
            
            /* Dark mode: catch-all for any white lines/paths */
            .js-plotly-plot line[stroke="white"],
            .js-plotly-plot line[stroke="#ffffff"],
            .js-plotly-plot line[stroke="rgb(255, 255, 255)"],
            .js-plotly-plot path[stroke="white"],
            .js-plotly-plot path[stroke="#ffffff"],
            .js-plotly-plot path[stroke="rgb(255, 255, 255)"] {
                stroke: #808080 !important;
            }
            
            /* 3Dmol viewer dark mode */
            canvas {
                background-color: #1a1a1a !important;
            }
            
            .viewer_3Dmoljs {
                background-color: #1a1a1a !important;
            }
            
            .viewer_3Dmoljs canvas {
                background-color: #1a1a1a !important;
            }
            
            .text-primary {
                color: #5a9fd4 !important;
            }
            
            .text-info {
                color: #4dd4f1 !important;
            }
            
            .text-success {
                color: #6abf69 !important;
            }
            
            .text-warning {
                color: #f9ca24 !important;
            }
            
            .text-danger {
                color: #e76f51 !important;
            }
            
            /* Light mode overrides when active */
            body.light-mode {
                background-color: #ffffff !important;
                color: #212529 !important;
            }
            
            body.light-mode * {
                color: #212529 !important;
            }
            
            /* Light mode colored text - keep colors visible */
            body.light-mode .text-primary {
                color: #007bff !important;
            }
            
            body.light-mode .text-info {
                color: #17a2b8 !important;
            }
            
            body.light-mode .text-success {
                color: #28a745 !important;
            }
            
            body.light-mode .text-warning {
                color: #ffc107 !important;
            }
            
            body.light-mode .text-danger {
                color: #dc3545 !important;
            }
            
            body.light-mode .container-fluid {
                background-color: #ffffff !important;
            }
            
            body.light-mode .well, body.light-mode .card {
                background-color: #f8f9fa !important;
                border: 1px solid #dee2e6 !important;
                color: #212529 !important;
            }
            
            body.light-mode .card-header {
                background-color: #e9ecef !important;
                border-bottom: 1px solid #dee2e6 !important;
                color: #212529 !important;
            }
            
            body.light-mode .sidebar {
                background-color: #f8f9fa !important;
                border-right: 1px solid #dee2e6 !important;
                color: #212529 !important;
            }
            
            body.light-mode h1, body.light-mode h2, body.light-mode h3, body.light-mode h4, body.light-mode h5, body.light-mode h6, 
            body.light-mode p, body.light-mode span, body.light-mode div, body.light-mode label {
                color: #212529 !important;
            }
            
            body.light-mode .nav-tabs {
                border-bottom: 1px solid #dee2e6 !important;
            }
            
            body.light-mode .nav-tabs .nav-link {
                background-color: #f8f9fa !important;
                border: 1px solid #dee2e6 !important;
                color: #495057 !important;
            }
            
            body.light-mode .nav-tabs .nav-link:hover {
                background-color: #e9ecef !important;
                color: #212529 !important;
            }
            
            body.light-mode .nav-tabs .nav-link.active {
                background-color: #ffffff !important;
                border-bottom-color: transparent !important;
                color: #212529 !important;
            }
            
            body.light-mode input[type="text"], body.light-mode input[type="number"], body.light-mode select, body.light-mode textarea {
                background-color: #ffffff !important;
                border: 1px solid #ced4da !important;
                color: #212529 !important;
            }
            
            body.light-mode input[type="text"]:focus, body.light-mode input[type="number"]:focus, body.light-mode select:focus {
                background-color: #ffffff !important;
                border-color: #80bdff !important;
                color: #212529 !important;
            }
            
            body.light-mode .btn-primary {
                background-color: #007bff !important;
                border-color: #007bff !important;
                color: #ffffff !important;
            }
            
            body.light-mode .btn-outline-secondary {
                background-color: #ffffff !important;
                border-color: #6c757d !important;
                color: #6c757d !important;
            }
            
            body.light-mode .shiny-notification {
                background-color: #ffffff !important;
                border: 1px solid #dee2e6 !important;
                color: #212529 !important;
            }
            
            body.light-mode .shiny-data-grid {
                background-color: #ffffff !important;
            }
            
            body.light-mode .shiny-data-grid th {
                background-color: #e9ecef !important;
                color: #212529 !important;
                border-color: #dee2e6 !important;
            }
            
            body.light-mode .shiny-data-grid td {
                background-color: #ffffff !important;
                color: #212529 !important;
                border-color: #dee2e6 !important;
            }
            
            body.light-mode .shiny-data-grid tr:hover {
                background-color: #f8f9fa !important;
            }
            
            body.light-mode pre {
                background-color: #f8f9fa !important;
                border: 1px solid #dee2e6 !important;
                color: #28a745 !important;
            }
            
            /* Light mode Plotly - override bg, text, and gridlines */
            body.light-mode .plotly-graph-div,
            body.light-mode .js-plotly-plot {
                background-color: #ffffff !important;
            }
            
            /* Light mode: make text dark */
            body.light-mode .js-plotly-plot text {
                fill: #212529 !important;
            }
            
            /* Light mode: hover labels - keep them readable */
            body.light-mode .js-plotly-plot .hoverlayer .hovertext,
            body.light-mode .js-plotly-plot .hoverlayer text {
                fill: #212529 !important;
            }
            
            body.light-mode .js-plotly-plot .hoverlayer path {
                fill: #ffffff !important;
                stroke: #666666 !important;
            }
            
            /* Light mode: make gridlines visible on white */
            body.light-mode .js-plotly-plot .gridlayer line,
            body.light-mode .js-plotly-plot .gridlayer path {
                stroke: #d0d0d0 !important;
            }
            
            /* Light mode: zero lines AND axis borders */
            body.light-mode .js-plotly-plot .zerolinelayer line,
            body.light-mode .js-plotly-plot line.xlines-above,
            body.light-mode .js-plotly-plot line.ylines-above,
            body.light-mode .js-plotly-plot .xlines-below line,
            body.light-mode .js-plotly-plot .ylines-below line,
            body.light-mode .js-plotly-plot .overaxes-above line {
                stroke: #606060 !important;
            }
            
            /* Light mode: catch all for axis elements */
            body.light-mode .js-plotly-plot g.xaxislayer-above path,
            body.light-mode .js-plotly-plot g.yaxislayer-above path,
            body.light-mode .js-plotly-plot g.xaxislayer-above line,
            body.light-mode .js-plotly-plot g.yaxislayer-above line {
                stroke: #606060 !important;
            }
            
            body.light-mode canvas {
                background-color: #ffffff !important;
            }
            
            body.light-mode .viewer_3Dmoljs {
                background-color: #ffffff !important;
            }
            
            /* Theme toggle button */
            #theme_toggle {
                position: fixed;
                top: 10px;
                right: 20px;
                z-index: 9999;
                border-radius: 20px;
                padding: 8px 16px;
                font-size: 16px;
            }
            
            /* Structure info box - adapt to theme */
            #structure_info > div {
                background-color: #2d2d2d !important;
                border-left-color: #5a9fd4 !important;
            }
            
            #structure_info > div * {
                color: #e0e0e0 !important;
            }
            
            body.light-mode #structure_info > div {
                background-color: #f8f9fa !important;
                border-left-color: #0d6efd !important;
            }
            
            body.light-mode #structure_info > div * {
                color: #212529 !important;
            }
            
            body.light-mode #structure_info > div span {
                color: #0d6efd !important;
            }
        """),
        ui.tags.script("""
            // Check system preference on load
            function getSystemTheme() {
                return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
            }
            
            // Apply theme
            function applyTheme(theme) {
                if (theme === 'light') {
                    document.body.classList.add('light-mode');
                } else {
                    document.body.classList.remove('light-mode');
                }
                localStorage.setItem('theme', theme);
                
                // Update button text
                const btn = document.getElementById('theme_toggle');
                if (btn) {
                    btn.textContent = theme === 'dark' ? '☀️ Light' : '🌙 Dark';
                }
            }
            
            // Initialize theme from localStorage or system preference
            window.addEventListener('DOMContentLoaded', function() {
                const savedTheme = localStorage.getItem('theme') || getSystemTheme();
                applyTheme(savedTheme);
            });
            
            // Listen for theme changes WITHOUT reload
            Shiny.addCustomMessageHandler('theme-change', function(message) {
                applyTheme(message.theme);
            });
            
            // Listen for open-url messages (e.g., file browser)
            Shiny.addCustomMessageHandler('open-url', function(message) {
                window.open(message.url, '_blank');
            });
        """)
    ),
    # Header with title, theme toggle, and download
    ui.div(
        ui.h1("🧬 ClustalDM - Interactive Clustering Analysis", style="display: inline-block; margin: 20px;"),
        ui.download_button("download_all", "📥 Download Results", class_="btn btn-sm btn-success", style="margin-right: 10px;"),
        ui.input_action_button("theme_toggle", "🌙 Dark", class_="btn btn-sm btn-outline-secondary"),
        style="position: relative;"
    ),
    
    ui.layout_sidebar(
        ui.sidebar(
            ui.h4("Upload HDOCK Results"),
            ui.div(
                ui.input_file("hdock_out_files", "1. HDOCK Output (.out):", 
                             multiple=True,
                             accept=[".out"]),
                ui.tags.small("Upload one or more HDOCK .out files", 
                             style="color: #6c757d; display: block; margin-top: 5px; margin-bottom: 15px;"),
                style="margin-bottom: 10px;"
            ),
            ui.div(
                ui.input_file("receptor_pdb", "2. Receptor PDB:", 
                             accept=[".pdb"]),
                ui.tags.small("Upload receptor structure (protein/DNA)", 
                             style="color: #6c757d; display: block; margin-top: 5px; margin-bottom: 15px;"),
                style="margin-bottom: 10px;"
            ),
            ui.div(
                ui.input_file("ligand_pdb", "3. Ligand PDB:", 
                             accept=[".pdb"]),
                ui.tags.small("Upload ligand structure (protein/DNA/RNA)", 
                             style="color: #6c757d; display: block; margin-top: 5px; margin-bottom: 10px;"),
                style="margin-bottom: 10px;"
            ),
            ui.output_ui("upload_status"),
            ui.hr(),
            ui.h5("Analysis Options"),
            ui.input_text("motif_residues", "Motif Residues (optional):", 
                         value="",
                         placeholder="e.g., A:10-20,B:30-35"),
            ui.tags.small("Format: Chain:Start-End, comma-separated", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 5px;"),
            ui.input_radio_buttons(
                "sasa_mode",
                "ΔSASA Calculation:",
                choices={
                    "combined": "Combined (Total Interface)",
                    "split": "Split by Chains"
                },
                selected="combined"
            ),
            ui.tags.small("Atomic-level calculation with 1.4 Å probe radius (water molecule). Split mode requires motif residues. Negative ΔSASA indicates conformational changes or cavities.", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 10px;"),
            ui.hr(),
            ui.h5("Pose Extraction"),
            ui.input_slider("max_poses", "Max Poses per File:", 
                           min=10, max=4392, value=100, step=10),
            ui.tags.small("Limit poses extracted from each HDOCK output file (lower = faster testing)", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 10px;"),
            ui.hr(),
            ui.h5("Clustering (HDBSCAN)"),
            ui.input_slider("min_cluster_size", "Min Cluster Size:", 
                           min=2, max=20, value=5, step=1),
            ui.input_slider("min_samples_hdb", "Min Samples (noise threshold):", 
                           min=1, max=10, value=2, step=1),
            ui.input_checkbox("filter_duplicates", "Filter Near-Duplicates", value=True),
            ui.tags.small("Remove very similar structures before clustering (faster, cleaner trees)", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 5px;"),
            ui.input_numeric("duplicate_threshold", "Duplicate Threshold:", 
                           value=0.0001, min=0.00001, max=0.1, step=0.0001),
            ui.tags.small("HDBSCAN automatically finds optimal clustering without eps parameter", 
                         style="color: #6c757d; display: block; margin-top: 5px; margin-bottom: 10px;"),
            ui.hr(),
            ui.input_action_button("run_clustering", "🚀 Run Analysis", class_="btn-primary"),
            ui.input_action_button("recluster", "🔄 Re-Cluster (fast)", class_="btn-secondary", 
                                  style="margin-top: 10px; width: 100%;"),
            ui.tags.small("Re-Cluster uses cached distance matrices (instant!)", 
                         style="color: #6c757d; display: block; margin-top: 5px;"),
            ui.hr(),
            ui.h5("Status"),
            ui.output_text_verbatim("status"),
            width=320
        ),
        
        ui.navset_tab(
            ui.nav_panel(
                "📊 Overview",
                # Summary cards - 3 column layout
                ui.row(
                    ui.column(4,
                        ui.card(
                            ui.card_header("Cluster Summary"),
                            ui.output_ui("cluster_summary_cards")
                        )
                    ),
                    ui.column(4,
                        ui.card(
                            ui.card_header("Quality Metrics - Structure"),
                            ui.output_ui("structure_quality_cards")
                        )
                    ),
                    ui.column(4,
                        ui.card(
                            ui.card_header("Quality Metrics - Interface"),
                            ui.output_ui("interface_quality_cards")
                        )
                    )
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Cluster Summary Table with RMSD QC"),
                    ui.input_text("cluster_search", "Search table:", placeholder="Type to filter..."),
                    ui.output_data_frame("cluster_summary_table")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Cluster Size Distribution"),
                    ui.output_ui("cluster_size_plot")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Cluster Stability"),
                    ui.markdown("""
                    **HDBSCAN Cluster Persistence Scores**  
                    - Higher values = more stable/persistent clusters  
                    - Bar color intensity = cluster size (log scale)  
                    - Based on λ-integral from HDBSCAN's condensed tree
                    """),
                    ui.output_ui("cluster_stability_plot"),
                    ui.output_ui("cluster_stability_info")  # Warning text from plot metadata
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Motif Match Score"),
                    ui.markdown("""
                    **Motif Overlap Quality**  
                    - Shows % of motif residues making contact with DNA/RNA
                    - Higher values = better match to specified motif
                    - Clusters sorted by average match score (best first)
                    - Color: Green (high) to Red (low)
                    """),
                    ui.output_ui("motif_match_plot")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Pose Diversity within Clusters"),
                    ui.markdown("""
                    **Intra-Cluster RMSD Distribution**  
                    - Each point = one pose  
                    - X-axis: RMSD from cluster medoid (Å)  
                    - Y-axis: Cluster (sorted by size)  
                    - Color: ΔSASA (buried surface area)
                    """),
                    ui.output_ui("cluster_diversity_plot")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("HDBSCAN Condensed Tree"),
                    ui.markdown("""
                    **Cluster Formation & Stability** (Right-click image to save)  
                    - Vertical: Density (λ) | Horizontal: Cluster Size  
                    - Longer bars = more stable clusters  
                    - Color-coded by cluster
                    """),
                    ui.output_ui("hdbscan_condensed_tree_plot")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("HDBSCAN Single Linkage Tree"),
                    ui.markdown("""
                    **Minimum Spanning Tree** (Right-click image to save)  
                    - HDBSCAN's hierarchical structure visualization  
                    - Scales well to 1000+ poses  
                    - Color gradient shows merge distances
                    """),
                    ui.output_ui("hdbscan_linkage_tree_plot")
                )
            ),
            ui.nav_panel(
                "📊 Interactive Maps",
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.h5("Visualization Method"),
                        ui.input_radio_buttons(
                            "scatter_method",
                            "Dimensionality Reduction:",
                            choices={
                                "mds": "MDS (Multidimensional Scaling)",
                                "tsne": "t-SNE (t-Distributed Stochastic Neighbor Embedding)",
                                "pca": "PCA (Principal Component Analysis)"
                            },
                            selected="mds"
                        ),
                        ui.hr(),
                        ui.markdown("""
                        **Method Comparison:**
                        - **MDS**: Preserves pairwise distances (best for distance matrices)
                        - **t-SNE**: Emphasizes local structure, reveals clusters
                        - **PCA**: Linear projection, fast but may miss non-linear patterns
                        """),
                        width=300
                    ),
                    ui.card(
                        ui.card_header("Jaccard Distance Clustering"),
                        ui.output_ui("interactive_scatter")
                    ),
                    ui.hr(),
                    ui.card(
                        ui.card_header("Jaccard Distance Matrix"),
                        ui.input_radio_buttons(
                            "jaccard_heatmap_level",
                            "Heatmap Level:",
                            choices={
                                "pose": "Pose-Level (all structures)",
                                "cluster": "Cluster-Level (inter-cluster distances)",
                                "intra": "Intra-Cluster (cohesiveness validation)"
                            },
                            selected="cluster",
                            inline=True
                        ),
                        ui.panel_conditional(
                            "input.jaccard_heatmap_level === 'intra'",
                            ui.input_select("jaccard_intra_cluster", "Select Cluster:", choices={})
                        ),
                        ui.input_checkbox("jaccard_dendrogram", "Show Dendrogram", value=False),
                        ui.output_ui("jaccard_heatmap")
                    ),
                    ui.hr(),
                    ui.card(
                        ui.card_header("RMSD Distance Matrix"),
                        ui.input_radio_buttons(
                            "rmsd_heatmap_level",
                            "Heatmap Level:",
                            choices={
                                "pose": "Pose-Level (all structures)",
                                "cluster": "Cluster-Level (inter-cluster distances)",
                                "intra": "Intra-Cluster (cohesiveness validation)"
                            },
                            selected="cluster",
                            inline=True
                        ),
                        ui.panel_conditional(
                            "input.rmsd_heatmap_level === 'intra'",
                            ui.input_select("rmsd_intra_cluster", "Select Cluster:", choices={})
                        ),
                        ui.input_checkbox("rmsd_dendrogram", "Show Dendrogram", value=False),
                        ui.output_ui("rmsd_heatmap")
                    )
                )
            ),
            ui.nav_panel(
                "🔥 Binding Hotspots",
                ui.card(
                    ui.card_header("Consensus Residues by Cluster"),
                    ui.input_text("consensus_search", "Search Consensus Residues:", placeholder="Search cluster, residues, or frequency...", width="100%"),
                    ui.output_data_frame("consensus_residues_table")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("All Chains Combined (Color-coded: Red=A, Blue=B)"),
                    ui.output_ui("hotspot_plot_combined")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Per-Chain Hotspot Analysis"),
                    ui.output_ui("hotspot_plots_by_chain")
                ),
                ui.hr(),
                ui.output_ui("hotspot_note"),
                ui.hr(),
                ui.card(
                    ui.card_header("Contact Residue Heatmap by Cluster"),
                    ui.markdown("""
                        **Intra-cluster Contact Pattern Analysis**: Shows which protein residues each structure within a cluster contacts.
                        - Colors represent actual contact distances in Ångströms (red = closer, blue = farther)
                        - **Vertical patterns** = conserved residues, **Diagonal** = shifted binding modes
                        - Structures automatically ordered by contact pattern similarity
                    """),
                    ui.input_select(
                        "contact_heatmap_cluster",
                        "Select Cluster:",
                        choices=[],
                        width="100%"
                    ),
                    ui.output_ui("contact_residue_heatmap")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Per-Cluster Hotspot Analysis"),
                    ui.input_select("cluster_select_hotspot", "Select Cluster:", choices=[], width="300px"),
                    ui.output_ui("cluster_hotspot_plot")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Hotspot Data"),
                    ui.output_data_frame("hotspot_table")
                )
            ),
            ui.nav_panel(
                "📋 Structure Details",
                ui.card(
                    ui.card_header("Interface statistics for all structures"),
                    ui.div(
                        ui.input_select("filter_cluster", "Filter by Cluster:", 
                                       choices={}, 
                                       selected="",
                                       width="300px"),
                        ui.tags.small("Select a cluster to view only structures in that cluster. Default: Cluster 1",
                                     style="color: #6c757d; margin-left: 10px;")
                    ),
                    ui.output_ui("stats_table_info"),
                    ui.output_data_frame("stats_table")
                ),
                ui.hr(),
                ui.card(
                    ui.card_header("Export PDB Structures"),
                    ui.div(
                        ui.input_select("export_cluster", "Select Cluster to Export:",
                                       choices={},
                                       selected="",
                                       width="300px"),
                        ui.download_button("export_cluster_pdbs", "📦 Export Cluster PDBs", class_="btn-success", style="margin-top: 10px;"),
                        ui.tags.small("Downloads a ZIP file containing all PDB files for the selected cluster",
                                     style="color: #6c757d; display: block; margin-top: 10px;")
                    ),
                    ui.output_ui("export_status")
                )
            ),
            ui.nav_panel(
                "🔬 3D Viewer",
                ui.card(
                    ui.card_header("Interactive 3D Structure Viewer (Molstar)"),
                    ui.layout_columns(
                        ui.input_select("viewer_cluster", "Filter by Cluster:",
                                       choices={},
                                       selected="",
                                       width="200px"),
                        ui.input_radio_buttons("structure_type", "Show:",
                                              choices={"all": "All Structures", "medoid": "Medoids Only"},
                                              selected="all",
                                              inline=True),
                        ui.input_select("structure_select", "Select Structure:", 
                                       choices=[], 
                                       width="300px"),
                        col_widths=[3, 4, 5]
                        ),
                        ui.output_ui("structure_info"),
                        ui.output_ui("structure_viewer")
                    )
                )
        )
    )
)


# ============================================================================
# Helper function for parallel SASA computation (must be at module level for pickling)
# ============================================================================

def _compute_sasa_worker(args):
    """Worker function for parallel SASA computation (module-level for pickling)"""
    idx, pdb_file, contact_pairs, split_mode_val, motif_dict_val, protein_pdb_path, nucleic_pdb_path = args
    
    # Import here to ensure each process has the module
    from core.analysis import InterfaceAnalyzer
    import os
    import numpy as np
    import freesasa
    
    try:
        # Get protein and nucleic chains from contact data
        contacts = contact_pairs[idx]
        protein_chains = set()
        nucleic_chains = set()
        
        for contact in contacts:
            prot_res = contact['protein_residue']
            nuc_res = contact['nucleic_residue']
            
            if ':' in prot_res:
                protein_chains.add(prot_res.split(':')[0])
            if ':' in nuc_res:
                nucleic_chains.add(nuc_res.split(':')[0])
        
        # Compute reference SASA in this worker (freesasa objects can't be pickled)
        precomputed_protein = None
        precomputed_nucleic = None
        precomputed_pdbs = None
        
        if protein_pdb_path and nucleic_pdb_path:
            freesasa.setVerbosity(freesasa.silent)
            params = freesasa.Parameters()
            params.setProbeRadius(1.4)
            
            structure_protein = freesasa.Structure(protein_pdb_path)
            result_protein = freesasa.calc(structure_protein, params)
            sasa_protein = result_protein.totalArea()
            
            structure_nucleic = freesasa.Structure(nucleic_pdb_path)
            result_nucleic = freesasa.calc(structure_nucleic, params)
            sasa_nucleic = result_nucleic.totalArea()
            
            precomputed_protein = {
                'sasa': sasa_protein,
                'structure': structure_protein,
                'result': result_protein
            }
            precomputed_nucleic = {
                'sasa': sasa_nucleic,
                'structure': structure_nucleic,
                'result': result_nucleic
            }
            precomputed_pdbs = {
                'protein_pdb': protein_pdb_path,
                'nucleic_pdb': nucleic_pdb_path
            }
        
        # Calculate ΔSASA
        sasa_result = InterfaceAnalyzer.calculate_delta_sasa(
            pdb_file,
            list(protein_chains),
            list(nucleic_chains),
            motif_residues=motif_dict_val,
            split_by_chain=split_mode_val,
            precomputed_protein_sasa=precomputed_protein,
            precomputed_nucleic_sasa=precomputed_nucleic,
            precomputed_structures=precomputed_pdbs
        )
        
        result_dict = {
            'structure': os.path.basename(pdb_file),
            'total_delta_sasa': sasa_result['total_delta_sasa'],
            'motif_delta_sasa': sasa_result.get('motif_delta_sasa'),
            'motif_sasa_protein': sasa_result.get('motif_sasa_protein'),
            'sasa_complex': sasa_result['sasa_complex'],
            'sasa_protein': sasa_result['sasa_protein'],
            'sasa_nucleic': sasa_result['sasa_nucleic']
        }
        
        # Add per-chain ΔSASA if in split mode
        if split_mode_val and 'chain_sasa' in sasa_result:
            for chain_id, chain_sasa_val in sasa_result['chain_sasa'].items():
                result_dict[f'chain_{chain_id}_sasa'] = chain_sasa_val
        
        return result_dict
        
    except Exception as e:
        import traceback
        error_detail = f"{type(e).__name__}: {str(e)}"
        return {
            'structure': os.path.basename(pdb_file),
            'total_delta_sasa': np.nan,
            'motif_delta_sasa': np.nan,
            'motif_sasa_protein': np.nan,
            'sasa_complex': np.nan,
            'sasa_protein': np.nan,
            'sasa_nucleic': np.nan,
            'error': error_detail
        }


# ============================================================================
# Server Logic
# ============================================================================

def server(input, output, session):
    global clust, pdb_dir_path
    
    # Reactive values
    interface_stats = reactive.Value(None)
    hotspots = reactive.Value(None)
    cluster_summary = reactive.Value(None)
    rmsd_stats = reactive.Value(None)
    rmsd_matrix = reactive.Value(None)
    sasa_data = reactive.Value(None)  # ΔSASA results
    current_motif = reactive.Value(None)  # Store motif for SASA visualization
    analysis_complete = reactive.Value(False)  # Track if analysis is done
    status_log = reactive.Value([])  # Terminal-style output log
    current_theme = reactive.Value('dark')  # Track current theme
    
    # Upload-specific reactive values
    temp_dir = reactive.Value(None)  # Temporary directory for uploaded files
    uploaded_files = reactive.Value([])  # List of uploaded .out files
    receptor_pdb_file = reactive.Value(None)  # Uploaded receptor PDB
    ligand_pdb_file = reactive.Value(None)  # Uploaded ligand PDB
    parsed_pdbs = reactive.Value([])  # List of generated PDB files
    auto_run_trigger = reactive.Value(0)  # Trigger for auto-run (increment to trigger)
    
    def log_status(message):
        """Add message to status log"""
        current_log = status_log.get()
        current_log.append(message)
        status_log.set(current_log)
    
    @output
    @render.ui
    def upload_status():
        """Show upload status"""
        out_files = uploaded_files.get()
        receptor = receptor_pdb_file.get()
        ligand = ligand_pdb_file.get()
        pdbs = parsed_pdbs.get()
        
        status_items = []
        
        # Check uploads
        if out_files:
            status_items.append(ui.tags.li(f"✓ {len(out_files)} HDOCK .out file(s)", style="color: #28a745;"))
        else:
            status_items.append(ui.tags.li("⚠ No .out files uploaded", style="color: #6c757d;"))
        
        if receptor:
            status_items.append(ui.tags.li("✓ Receptor PDB", style="color: #28a745;"))
        else:
            status_items.append(ui.tags.li("⚠ No receptor PDB", style="color: #6c757d;"))
        
        if ligand:
            status_items.append(ui.tags.li("✓ Ligand PDB", style="color: #28a745;"))
        else:
            status_items.append(ui.tags.li("⚠ No ligand PDB", style="color: #6c757d;"))
        
        # Show ready status
        if out_files and receptor and ligand:
            if pdbs:
                status_items.append(ui.tags.li(f"✓ {len(pdbs)} poses generated - Ready for analysis!", 
                                              style="color: #007bff; font-weight: bold;"))
            else:
                status_items.append(ui.tags.li("👉 Ready! Click 'Start Analysis' to extract ALL poses and cluster", 
                                              style="color: #007bff; font-weight: bold;"))
        
        return ui.tags.ul(*status_items, style="list-style: none; padding-left: 0;")
    
    @reactive.effect
    @reactive.event(input.hdock_out_files)
    def handle_out_upload():
        """Handle HDOCK .out file uploads"""
        files_info = input.hdock_out_files()
        if not files_info:
            return
        
        try:
            # Create temporary directory if it doesn't exist
            if temp_dir.get() is None:
                tmp = tempfile.mkdtemp(prefix="clustaldm_")
                temp_dir.set(tmp)
                log_status(f"Created temp directory: {tmp}")
            
            tmp = temp_dir.get()
            out_files = []
            
            for file_info in files_info:
                file_path = os.path.join(tmp, file_info['name'])
                shutil.copy(file_info['datapath'], file_path)
                out_files.append(file_path)
                log_status(f"Uploaded .out: {file_info['name']}")
            
            uploaded_files.set(out_files)
            # Clear any previous generated PDBs
            parsed_pdbs.set([])
            
            ui.notification_show(
                f"✓ Uploaded {len(out_files)} HDOCK .out file(s)",
                type="success",
                duration=3
            )
            
        except Exception as e:
            log_status(f"Upload error: {str(e)}")
            ui.notification_show(f"Upload failed: {str(e)}", type="error", duration=5)
    
    @reactive.effect
    @reactive.event(input.receptor_pdb)
    def handle_receptor_upload():
        """Handle receptor PDB upload"""
        file_info = input.receptor_pdb()
        if not file_info:
            return
        
        try:
            if temp_dir.get() is None:
                tmp = tempfile.mkdtemp(prefix="clustaldm_")
                temp_dir.set(tmp)
            
            tmp = temp_dir.get()
            # Save with a standard name
            receptor_path = os.path.join(tmp, "receptor.pdb")
            shutil.copy(file_info[0]['datapath'], receptor_path)
            receptor_pdb_file.set(receptor_path)
            log_status(f"Uploaded receptor: {file_info[0]['name']}")
            
            ui.notification_show(
                f"✓ Receptor PDB uploaded: {file_info[0]['name']}",
                type="success",
                duration=3
            )
            
        except Exception as e:
            log_status(f"Receptor upload error: {str(e)}")
            ui.notification_show(f"Receptor upload failed: {str(e)}", type="error", duration=5)
    
    @reactive.effect
    @reactive.event(input.ligand_pdb)
    def handle_ligand_upload():
        """Handle ligand PDB upload"""
        file_info = input.ligand_pdb()
        if not file_info:
            return
        
        try:
            if temp_dir.get() is None:
                tmp = tempfile.mkdtemp(prefix="clustaldm_")
                temp_dir.set(tmp)
            
            tmp = temp_dir.get()
            # Save with a standard name
            ligand_path = os.path.join(tmp, "ligand.pdb")
            shutil.copy(file_info[0]['datapath'], ligand_path)
            ligand_pdb_file.set(ligand_path)
            log_status(f"Uploaded ligand: {file_info[0]['name']}")
            
            ui.notification_show(
                f"✓ Ligand PDB uploaded: {file_info[0]['name']}",
                type="success",
                duration=3
            )
            
        except Exception as e:
            log_status(f"Ligand upload error: {str(e)}")
            ui.notification_show(f"Ligand upload failed: {str(e)}", type="error", duration=5)
    
    @reactive.effect
    def auto_run_validation():
        """Auto-run analysis when all three files are uploaded and validated"""
        # Watch the reactive values
        out_files = uploaded_files.get()
        receptor = receptor_pdb_file.get()
        ligand = ligand_pdb_file.get()
        
        # Check if all files are uploaded
        if not out_files or not receptor or not ligand:
            return  # Not ready yet
        
        # Check if analysis is already running/complete to avoid re-triggering
        if analysis_complete.get():
            return
        
        # Extract job IDs from filenames to validate they match
        def extract_job_id(filepath):
            """Extract HDOCK job ID from filename (e.g., hdock_69083e48ce879.out -> 69083e48ce879)"""
            import re
            basename = os.path.basename(filepath)
            # Match pattern: hdock_JOBID or rec_JOBID or lig_JOBID
            match = re.search(r'(?:hdock|rec|lig)_([a-f0-9]+)', basename)
            return match.group(1) if match else None
        
        # Get job IDs
        out_job_ids = set()
        for out_file in out_files:
            job_id = extract_job_id(out_file)
            if job_id:
                out_job_ids.add(job_id)
        
        receptor_job_id = extract_job_id(receptor)
        ligand_job_id = extract_job_id(ligand)
        
        # Validate all files are from the same job
        if receptor_job_id and ligand_job_id and out_job_ids:
            # Check if receptor and ligand have the same job ID
            if receptor_job_id == ligand_job_id:
                # Check if at least one .out file matches
                if receptor_job_id in out_job_ids:
                    log_status("")
                    log_status("✓ All files validated: Same HDOCK job ID detected")
                    log_status(f"  Job ID: {receptor_job_id}")
                    log_status("🚀 Auto-starting analysis...")
                    ui.notification_show(
                        "✓ Files validated! Auto-starting analysis...",
                        type="success",
                        duration=3
                    )
                    # Trigger analysis by incrementing the trigger value
                    # This will fire the run_analysis reactive effect
                    auto_run_trigger.set(auto_run_trigger.get() + 1)
                else:
                    log_status("⚠ Warning: .out file job ID doesn't match receptor/ligand")
                    ui.notification_show(
                        "⚠️ Job ID mismatch. Click 'Run Analysis' manually if files are correct.",
                        type="warning",
                        duration=5
                    )
            else:
                log_status("⚠ Warning: Receptor and ligand have different job IDs")
                ui.notification_show(
                    "⚠️ Receptor and ligand may be from different jobs. Verify before running.",
                    type="warning",
                    duration=5
                )
    
    # Import visualization functions for export
    from visualization.interactive import (
        create_interactive_scatter,
        create_cluster_size_distribution,
        create_cluster_stability_plot,
        create_cluster_diversity_plot,
        create_distance_heatmap,
        create_distance_heatmap_with_dendrogram,
        create_hdbscan_condensed_tree,
        create_hdbscan_tree,
        compute_cluster_distance_matrix,
        compute_intra_cluster_matrices,
        create_cluster_distance_heatmap,
        create_intra_cluster_heatmaps,
        create_single_intra_cluster_heatmap
    )
    
    @reactive.effect
    @reactive.event(input.theme_toggle)
    async def toggle_theme():
        """Toggle between light and dark themes - CSS only, no reload"""
        theme = current_theme.get()
        new_theme = 'light' if theme == 'dark' else 'dark'
        current_theme.set(new_theme)
        
        # Apply theme via CSS only (no page reload)
        await session.send_custom_message('theme-change', {'theme': new_theme})
    
    # Cleanup temporary directory on session end
    @reactive.effect
    def cleanup():
        """Clean up temporary files when session ends"""
        def cleanup_temp():
            tmp = temp_dir.get()
            if tmp and os.path.exists(tmp):
                try:
                    shutil.rmtree(tmp)
                    print(f"Cleaned up temp directory: {tmp}")
                except Exception as e:
                    print(f"Error cleaning up temp directory: {e}")
        
        session.on_ended(cleanup_temp)
    
    # Download all results as ZIP
    @render.download(filename=lambda: f"clustaldm_results_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.zip")
    def download_all():
        """Generate ZIP file with all analysis results"""
        if not analysis_complete.get():
            return None
        
        # Create temporary ZIP file
        zip_path = os.path.join(tempfile.gettempdir(), f"clustaldm_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.zip")
        
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Add cluster info text
            info_content = "=== CLUSTER INFORMATION ===\n\n"
            unique_labels = sorted(set(clust.labels))
            for label in unique_labels:
                cluster_name = f"Cluster {label}" if label != -1 else "Noise"
                members = [os.path.basename(clust.pdb_files[i]) for i, l in enumerate(clust.labels) if l == label]
                info_content += f"{cluster_name}: {len(members)} structures\n"
                for member in members:
                    info_content += f"  - {member}\n"
                info_content += "\n"
            zipf.writestr("cluster_info.txt", info_content)
            
            # Add CSV files
            if interface_stats.get() is not None:
                csv_buffer = interface_stats.get().to_csv(index=False)
                zipf.writestr("cluster_statistics.csv", csv_buffer)
            
            if cluster_summary.get() is not None:
                csv_buffer = cluster_summary.get().to_csv(index=False)
                zipf.writestr("cluster_summary.csv", csv_buffer)
            
            if rmsd_stats.get() is not None:
                csv_buffer = rmsd_stats.get().to_csv(index=False)
                zipf.writestr("cluster_rmsd_summary.csv", csv_buffer)
            
            if sasa_data.get() is not None:
                csv_buffer = sasa_data.get().to_csv(index=False)
                zipf.writestr("delta_sasa_results.csv", csv_buffer)
            
            # Add plots as HTML
            try:
                from visualization.interactive import (
                    create_interactive_scatter,
                    create_cluster_size_distribution,
                    create_distance_heatmap,
                    create_hdbscan_condensed_tree,
                    create_hdbscan_tree,
                    compute_cluster_distance_matrix,
                    compute_intra_cluster_matrices,
                    create_cluster_distance_heatmap,
                    create_intra_cluster_heatmaps,
                    create_single_intra_cluster_heatmap
                )
                
                # Cluster scatter plot
                if interface_stats.get() is not None:
                    fig = create_interactive_scatter(
                        clust,
                        interface_stats.get(),
                        rmsd_matrix.get(),
                        pdb_dir_path
                    )
                    zipf.writestr("cluster_visualization.html", fig.to_html())
                
                # Size distribution
                fig = create_cluster_size_distribution(
                    clust.labels,
                    duplicate_groups=getattr(clust, 'duplicate_groups', None),
                    rescued_cluster_ids=getattr(clust, 'rescued_cluster_ids', None)
                )
                zipf.writestr("cluster_size_distribution.html", fig.to_html())
                
                
                # Heatmaps
                fig = create_distance_heatmap(
                    clust.distance_matrix,
                    clust.labels,
                    title="Jaccard Distance Matrix",
                    pdb_files=clust.pdb_files
                )
                zipf.writestr("jaccard_heatmap.html", fig.to_html())
                
                if rmsd_matrix.get() is not None:
                    fig = create_distance_heatmap(
                        rmsd_matrix.get(),
                        clust.labels,
                        title="RMSD Distance Matrix",
                        pdb_files=clust.pdb_files
                    )
                    zipf.writestr("rmsd_heatmap.html", fig.to_html())
                    
            except Exception as e:
                print(f"Error generating plots for ZIP: {e}")
        
        return zip_path
    
    # Export cluster PDBs
    @render.download(
        filename=lambda: f"cluster_{input.export_cluster()}_pdbs_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.zip"
    )
    def export_cluster_pdbs():
        """Export all PDB files for a specific cluster as ZIP"""
        if not analysis_complete.get() or clust is None:
            return None
        
        selected_cluster = input.export_cluster()
        if not selected_cluster:
            return None
        
        try:
            cluster_id = int(selected_cluster)
        except ValueError:
            return None
        
        # Get structures in this cluster
        cluster_indices = [i for i, label in enumerate(clust.labels) if label == cluster_id]
        
        if not cluster_indices:
            ui.notification_show(f"No structures found in cluster {cluster_id}", type="warning")
            return None
        
        # Create ZIP file
        zip_path = os.path.join(
            tempfile.gettempdir(), 
            f"cluster_{cluster_id}_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.zip"
        )
        
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Add cluster info
            info_lines = [
                f"Cluster {cluster_id} - {len(cluster_indices)} structures\n",
                "="*50 + "\n\n"
            ]
            
            # Add statistics if available
            stats_df = interface_stats.get()
            if stats_df is not None:
                cluster_stats = stats_df[stats_df['cluster'] == cluster_id]
                if len(cluster_stats) > 0:
                    info_lines.append("Structure Statistics:\n")
                    info_lines.append("-"*50 + "\n")
                    for _, row in cluster_stats.iterrows():
                        info_lines.append(f"\n{row['structure']}:\n")
                        info_lines.append(f"  Protein contacts: {row['n_protein_residues']}\n")
                        info_lines.append(f"  Nucleic contacts: {row['n_nucleic_residues']}\n")
                        if 'total_delta_sasa' in row and pd.notna(row['total_delta_sasa']):
                            info_lines.append(f"  Total ΔSASA: {row['total_delta_sasa']:.2f} Å²\n")
                        if 'motif_delta_sasa' in row and pd.notna(row['motif_delta_sasa']):
                            info_lines.append(f"  Motif ΔSASA: {row['motif_delta_sasa']:.2f} Å²\n")
                        if row.get('representative', False):
                            info_lines.append(f"  ** CLUSTER REPRESENTATIVE **\n")
            
            zipf.writestr("cluster_info.txt", "".join(info_lines))
            
            # Add PDB files
            for idx in cluster_indices:
                pdb_path = clust.pdb_files[idx]
                if os.path.exists(pdb_path):
                    zipf.write(pdb_path, arcname=os.path.basename(pdb_path))
            
            ui.notification_show(
                f"Exported {len(cluster_indices)} PDB files from cluster {cluster_id}",
                type="message",
                duration=3
            )
        
        return zip_path
    
    @output
    @render.ui
    def export_status():
        """Show export status message"""
        if not analysis_complete.get():
            return ui.p("Run analysis first to enable export", style="color: #6c757d; margin-top: 10px;")
        
        selected = input.export_cluster()
        if selected:
            try:
                cluster_id = int(selected)
                n_structs = sum(1 for l in clust.labels if l == cluster_id)
                return ui.p(f"Ready to export {n_structs} PDB files", 
                           style="color: #28a745; margin-top: 10px; font-weight: bold;")
            except:
                pass
        return ui.p("Select a cluster above", style="color: #6c757d; margin-top: 10px;")
    
    @reactive.effect
    @reactive.event(input.run_clustering, auto_run_trigger)
    def run_analysis():
        global clust, pdb_dir_path
        
        # Import necessary modules at function start
        import os
        import MDAnalysis as mda
        import freesasa
        
        # Clear previous log
        status_log.set([])
        
        # Check if we need to generate poses first
        uploaded = uploaded_files.get()
        receptor = receptor_pdb_file.get()
        ligand = ligand_pdb_file.get()
        pdb_files = parsed_pdbs.get()
        
        # Validate uploads
        if not uploaded:
            ui.notification_show("⚠️ Please upload HDOCK .out file(s) first!", type="error", duration=5)
            return
        if not receptor:
            ui.notification_show("⚠️ Please upload receptor PDB file!", type="error", duration=5)
            return
        if not ligand:
            ui.notification_show("⚠️ Please upload ligand PDB file!", type="error", duration=5)
            return
        
        # Generate ALL poses if not already done
        if not pdb_files:
            log_status("=" * 50)
            log_status(f"📦 Extracting poses from HDOCK results (max {input.max_poses()} per file)...")
            log_status("=" * 50)
            
            try:
                tmp = temp_dir.get()
                all_pdbs = []
                
                for out_path in uploaded:
                    filename = os.path.basename(out_path)
                    log_status(f"Parsing: {filename}")
                    
                    try:
                        parser = HDOCKParser(out_path)
                        total_poses = parser.num_poses
                        log_status(f"  → Found {total_poses} poses in {filename}")
                        
                        output_dir = os.path.join(tmp, f"poses_{Path(filename).stem}")
                        os.makedirs(output_dir, exist_ok=True)
                        
                        # Get max_poses parameter from slider
                        max_poses_limit = input.max_poses()
                        
                        # Generate PDB files up to max_poses limit
                        generated = parser.generate_all_poses(
                            output_dir=output_dir,
                            max_poses=max_poses_limit,
                            receptor_pdb=receptor,
                            ligand_pdb=ligand
                        )
                        
                        all_pdbs.extend(generated)
                        log_status(f"  ✓ Extracted {len(generated)} poses")
                        
                    except Exception as e:
                        log_status(f"  ✗ Error parsing {filename}: {str(e)}")
                        ui.notification_show(
                            f"Error parsing {filename}: {str(e)}",
                            type="error",
                            duration=5
                        )
                        return
                
                parsed_pdbs.set(all_pdbs)
                pdb_files = all_pdbs
                log_status(f"✓ Total poses extracted: {len(all_pdbs)}")
                log_status("")
                
            except Exception as e:
                log_status(f"Pose generation error: {str(e)}")
                ui.notification_show(f"Pose generation failed: {str(e)}", type="error", duration=5)
                return
        
        # Show immediate notification for job submission
        ui.notification_show("🚀 Running clustering analysis...", type="message", duration=3)
        
        # Get HDBSCAN parameters
        min_cluster_size = input.min_cluster_size()
        min_samples_hdb = input.min_samples_hdb()
        
        # Show input acknowledgement
        log_status("="*50)
        log_status("🚀 Starting Clustering Analysis")
        log_status("="*50)
        log_status(f"Total structures: {len(pdb_files)}")
        log_status(f"HDBSCAN min_cluster_size: {min_cluster_size}")
        log_status(f"HDBSCAN min_samples: {min_samples_hdb}")
        log_status("")
        
        pdb_dir_path = temp_dir.get()  # Use temp directory
        
        try:
            log_status(f"✓ Analyzing {len(pdb_files)} structures")
            log_status("")
            
            # Initialize clusterer
            log_status("Initializing clusterer...")
            ui.notification_show("📊 Initializing clusterer...", type="message", duration=2)
            from core.analysis import InterfaceAnalyzer
            clust = PDBClusterer(pdb_files)
            log_status("✓ Clusterer initialized")
            log_status("")
            
            # Parse motif residues early for Jaccard screening
            motif_input = input.motif_residues().strip()
            motif_dict = InterfaceAnalyzer.parse_motif_residues(motif_input) if motif_input else None
            
            # STEP 1: Compute Jaccard matrix (chain-aware!) with smart filtering
            log_status("STEP 1: Computing Jaccard contact matrix...")
            ui.notification_show("🧬 STEP 1: Computing contact matrix...", type="message", duration=3)
            if motif_dict:
                log_status(f"  \ud83c\udfaf MOTIF MODE: Local search enabled")
                log_status(f"  Motif: {motif_input}")
                log_status(f"  \u2192 Off-target poses (no motif contact) will be filtered")
            else:
                log_status(f"  \ud83c\udf0d GLOBAL MODE: SASA pre-filtering enabled")
                log_status(f"  \u2192 Non-interacting poses (\u0394SASA < 0.01 \u00c5\u00b2) will be filtered")
            clust.compute_jaccard_contact_matrix(n_jobs=-1, motif_residues=motif_dict)  # Use all CPU cores
            log_status("✓ Jaccard matrix computed")
            ui.notification_show("✓ Contact matrix done", type="message", duration=2)
            log_status("")
            
            # STEP 2: Compute RMSD matrix for QC (nucleic P atoms, no alignment)
            log_status("STEP 2: Computing RMSD matrix (for QC)...")
            ui.notification_show("📏 STEP 2: Computing RMSD...", type="message", duration=3)
            # Load structures for RMSD calculation
            clust.load_structures()
            # Temporarily change selection for RMSD calculation
            original_selection = clust.selection
            clust.selection = 'nucleic and name P'
            clust.compute_distance_matrix(align_structures=False, n_jobs=-1)  # Parallel RMSD!
            rmsd_mat = clust.distance_matrix.copy()  # Save RMSD matrix
            rmsd_matrix.set(rmsd_mat)
            clust.selection = original_selection  # Restore original
            log_status("✓ RMSD matrix computed")
            print("\n=== CHECKPOINT: RMSD computation complete ===", flush=True)
            ui.notification_show("✓ RMSD done", type="message", duration=2)
            log_status("")
            
            # STEP 2.5: Compute ΔSASA for each structure
            print("\n=== CHECKPOINT: Starting SASA computation ===", flush=True)
            log_status("STEP 2.5: Computing ΔSASA (buried surface area)...")
            ui.notification_show("💧 STEP 2.5: Computing ΔSASA...", type="message", duration=3)
            
            # Store motif for later use in visualization (already parsed above)
            current_motif.set(motif_dict)
            
            # Only use motif and split mode if motif is actually provided
            use_motif = motif_dict is not None and len(motif_dict) > 0
            
            if use_motif:
                log_status(f"  Using motif: {motif_input}")
            
            # PERFORMANCE OPTIMIZATION: Compute protein and nucleic SASA only once
            # (they're identical across all poses, only complex changes)
            log_status("  Pre-computing isolated protein and nucleic acid SASA (once)...")
            # Extract reference protein/nucleic structures for SASA optimization
            # (workers will compute SASA from these files, not pickle C structures)
            protein_pdb_path = None
            nucleic_pdb_path = None
            
            if len(clust.pdb_files) > 0:
                try:
                    # Use first structure to extract protein/nucleic for reference
                    first_pdb = clust.pdb_files[0]
                    u = mda.Universe(first_pdb)
                    
                    # Create temp directory for separated structures
                    basename = os.path.basename(first_pdb).replace('.pdb', '')
                    debug_dir = os.path.join(os.path.dirname(first_pdb), 'sasa_debug')
                    os.makedirs(debug_dir, exist_ok=True)
                    
                    protein_pdb_path = os.path.join(debug_dir, f"{basename}_protein_reference.pdb")
                    nucleic_pdb_path = os.path.join(debug_dir, f"{basename}_nucleic_reference.pdb")
                    
                    # Extract protein and nucleic structures
                    protein_atoms = u.select_atoms("protein")
                    nucleic_atoms = u.select_atoms("nucleic")
                    protein_atoms.write(protein_pdb_path)
                    nucleic_atoms.write(nucleic_pdb_path)
                    
                    log_status(f"  ✓ Reference structures extracted (workers will compute SASA in parallel)")
                    
                except Exception as e:
                    log_status(f"  ⚠️ Could not extract reference structures, will compute per-structure: {e}")
                    protein_pdb_path = None
                    nucleic_pdb_path = None
            
            log_status("")
            log_status("STEP 2: Computing ΔSASA for all structures...")
            log_status(f"Mode: {'Motif + Split by chain' if input.sasa_mode() == 'split' and use_motif else 'Motif only' if use_motif else 'Total interface'}")
            log_status(f"Processing {len(clust.pdb_files)} structures in parallel...")
            ui.notification_show(f"⚙️ Processing {len(clust.pdb_files)} structures...", type="message", duration=3)
            
            # Read SASA mode parameter BEFORE parallel execution (can't pickle reactive context)
            split_mode = input.sasa_mode() == "split" and use_motif
            
            # Prepare arguments for parallel processing (only picklable objects)
            sasa_args = []
            for idx, pdb_file in enumerate(clust.pdb_files):
                sasa_args.append((
                    idx, pdb_file, clust.contact_residue_pairs, split_mode,
                    motif_dict if use_motif else None,
                    protein_pdb_path, nucleic_pdb_path
                ))
            
            # Run parallel computation with multiprocessing (worker function defined at module level)
            from multiprocessing import Pool, cpu_count
            log_status(f"  Using {cpu_count()} CPU cores for parallel SASA computation...")
            print(f"\n=== CHECKPOINT: About to create Pool with {cpu_count()} cores ===", flush=True)
            print(f"=== Processing {len(sasa_args)} structures ===", flush=True)
            
            # Use imap_unordered to get results as they complete (allows progress updates)
            sasa_results = []
            print("=== CHECKPOINT: Creating Pool ===", flush=True)
            pool = Pool(processes=cpu_count())
            try:
                print("=== CHECKPOINT: Pool created, starting imap_unordered ===", flush=True)
                total = len(sasa_args)
                # Use larger chunksize for better efficiency with large datasets
                # This reduces IPC overhead by batching work to processes
                chunksize = max(1, total // (cpu_count() * 4))
                print(f"=== Using chunksize={chunksize} for {total} structures ===", flush=True)
                for i, result in enumerate(pool.imap_unordered(_compute_sasa_worker, sasa_args, chunksize=chunksize), 1):
                    sasa_results.append(result)
                    # Show progress every 10% or every 20 structures (whichever is smaller)
                    progress_interval = min(20, max(1, total // 10))
                    if i % progress_interval == 0 or i == total:
                        percent = int(100 * i / total)
                        log_status(f"  Progress: {i}/{total} structures ({percent}%)")
                        print(f"SASA Progress: {percent}% ({i}/{total})", flush=True)
                
                print(f"\n=== CHECKPOINT: Loop complete, closing pool ===", flush=True)
            finally:
                pool.close()
                print("=== CHECKPOINT: Pool closed, joining ===", flush=True)
                pool.join()
                print("=== CHECKPOINT: Pool joined ===", flush=True)
            
            print(f"=== CHECKPOINT: Got {len(sasa_results)} results, creating DataFrame ===", flush=True)
            log_status(f"✓ SASA computation complete for {total} structures")
            ui.notification_show("✓ SASA computation complete", type="message", duration=2)
            
            # Store ΔSASA results
            print("=== CHECKPOINT: Calling pd.DataFrame() ===", flush=True)
            sasa_df = pd.DataFrame(sasa_results)
            print(f"=== CHECKPOINT: DataFrame created with {len(sasa_df)} rows ===", flush=True)
            
            # DEBUG: Check what we actually got
            log_status(f"DEBUG: Created SASA DataFrame with {len(sasa_df)} rows")
            log_status(f"DEBUG: SASA DataFrame columns: {sasa_df.columns.tolist()}")
            if len(sasa_df) > 0:
                log_status(f"DEBUG: First row total_delta_sasa: {sasa_df['total_delta_sasa'].iloc[0]}")
                log_status(f"DEBUG: First row structure: {sasa_df['structure'].iloc[0]}")
            
            sasa_data.set(sasa_df)  # Save to reactive value
            
            # Count errors
            n_errors = sum(1 for r in sasa_results if 'error' in r)
            if n_errors > 0:
                log_status(f"⚠️  ΔSASA computation completed with {n_errors} errors")
            else:
                log_status(f"✓ ΔSASA computed successfully for all {len(sasa_results)} structures")
            
            # Log motif ΔSASA range only if motif was used
            if use_motif and 'motif_delta_sasa' in sasa_df.columns:
                motif_values = sasa_df['motif_delta_sasa'].dropna()
                if len(motif_values) > 0:
                    log_status(f"  Motif ΔSASA range: {motif_values.min():.1f} - {motif_values.max():.1f} Å²")
            
            log_status(f"  Total ΔSASA range: {sasa_df['total_delta_sasa'].min():.1f} - {sasa_df['total_delta_sasa'].max():.1f} Å²")
            log_status("")
            
            # STEP 3: Cluster by Jaccard only
            log_status("STEP 3: Clustering by Jaccard distance...")
            ui.notification_show("🎯 STEP 3: Running HDBSCAN clustering...", type="message", duration=3)
            # Restore Jaccard matrix for clustering
            clust.distance_matrix = clust.jaccard_distance_matrix
            clust.metric_name = "Jaccard"
            
            filter_dupes = input.filter_duplicates()
            dup_threshold = input.duplicate_threshold()
            clust.cluster_hdbscan(
                min_cluster_size=min_cluster_size, 
                min_samples=min_samples_hdb,
                filter_duplicates=filter_dupes,
                duplicate_threshold=dup_threshold
            )
            
            n_clusters = len(set(clust.labels)) - (1 if -1 in clust.labels else 0)
            n_noise = sum(1 for l in clust.labels if l == -1)
            log_status(f"✓ Clustering complete: {n_clusters} clusters, {n_noise} noise points")
            ui.notification_show(f"✓ Found {n_clusters} clusters", type="message", duration=2)
            log_status("")
            
            # STEP 4: Compute RMSD stats per cluster
            log_status("STEP 4: Computing inter-cluster RMSD statistics...")
            ui.notification_show("📊 STEP 4: Computing statistics...", type="message", duration=2)
            rmsd_qc = compute_rmsd_qc(clust.labels, rmsd_mat)
            rmsd_stats.set(rmsd_qc)
            log_status("✓ RMSD QC computed")
            log_status("")
            
            # Get statistics
            log_status("Computing interface statistics...")
            stats_df = clust.get_interface_stats()
            
            # Mark representative structures based on medoid
            if 'medoid_idx' in rmsd_qc.columns:
                # Create a mapping of cluster -> medoid_idx
                medoid_map = dict(zip(rmsd_qc['cluster'], rmsd_qc['medoid_idx']))
                
                # Add structure index to stats_df
                stats_df['structure_idx'] = range(len(stats_df))
                
                # Mark representatives
                def is_representative(row):
                    cluster = row['cluster']
                    if cluster == -1:  # Noise points
                        return False
                    medoid_idx = medoid_map.get(cluster, -1)
                    return row['structure_idx'] == medoid_idx
                
                stats_df['representative'] = stats_df.apply(is_representative, axis=1)
                stats_df = stats_df.drop(columns=['structure_idx'])
            
            # Merge ΔSASA data if available
            sasa_df_current = sasa_data.get()
            if sasa_df_current is not None:
                log_status(f"DEBUG: SASA data has {len(sasa_df_current)} rows")
                log_status(f"DEBUG: Stats data has {len(stats_df)} rows")
                log_status(f"DEBUG: SASA columns: {sasa_df_current.columns.tolist()}")
                log_status(f"DEBUG: Stats columns before merge: {stats_df.columns.tolist()}")
                
                # Check structure names match
                log_status(f"DEBUG: Sample SASA structures: {sasa_df_current['structure'].head(3).tolist()}")
                log_status(f"DEBUG: Sample Stats structures: {stats_df['structure'].head(3).tolist()}")
                
                # Merge all SASA columns except 'structure' (which is the key)
                # Keep 'error' column if present to show which SASA computations failed
                sasa_cols_to_merge = [col for col in sasa_df_current.columns if col != 'structure']
                stats_df = stats_df.merge(
                    sasa_df_current[['structure'] + sasa_cols_to_merge], 
                    on='structure', 
                    how='left'
                )
                log_status(f"DEBUG: Stats columns after merge: {stats_df.columns.tolist()}")
                log_status(f"DEBUG: Sample total ΔSASA values: {stats_df['total_delta_sasa'].head(3).tolist()}")
                if 'motif_delta_sasa' in stats_df.columns:
                    log_status(f"DEBUG: Sample motif ΔSASA values: {stats_df['motif_delta_sasa'].head(3).tolist()}")
            
            # Add motif match scores if available
            if hasattr(clust, 'motif_match_scores') and clust.motif_match_scores is not None:
                match_scores = clust.motif_match_scores
                if len(match_scores) > 0 and len(match_scores) == len(stats_df):
                    stats_df['motif_match_pct'] = match_scores
                    log_status("✓ Motif match scores added to interface stats")
            
            interface_stats.set(stats_df)
            
            # Get hotspots
            hotspots_df = clust.get_binding_hotspots()
            hotspots.set(hotspots_df)
            log_status("✓ Hotspot analysis complete")
            log_status("")
            
            # Get cluster summary withwe RMSD QC
            summary_df = clust.get_cluster_summary(threshold=50)
            
            # Add RMSD stats to summary
            summary_with_rmsd = summary_df.merge(
                rmsd_qc, on='cluster', how='left'
            )
            
            # Add ΔSASA stats to summary (mean and std per cluster)
            sasa_df_current = sasa_data.get()
            if sasa_df_current is not None and 'total_delta_sasa' in sasa_df_current.columns:
                # Merge cluster labels with SASA data
                sasa_with_labels = sasa_df_current.copy()
                sasa_with_labels['cluster'] = clust.labels
                
                # Calculate cluster-level ΔSASA statistics
                agg_dict = {'total_delta_sasa': ['mean', 'std', 'max']}
                if 'motif_delta_sasa' in sasa_with_labels.columns:
                    agg_dict['motif_delta_sasa'] = ['mean', 'std']
                
                # Add per-chain ΔSASA columns if available
                chain_cols = [col for col in sasa_with_labels.columns if col.startswith('chain_') and col.endswith('_sasa')]
                for col in chain_cols:
                    agg_dict[col] = ['mean', 'std']
                
                sasa_cluster_stats = sasa_with_labels.groupby('cluster').agg(agg_dict).reset_index()
                
                # Flatten column names
                new_cols = ['cluster', 'mean_delta_sasa', 'std_delta_sasa', 'max_delta_sasa']
                if 'motif_delta_sasa' in agg_dict:
                    new_cols.extend(['mean_motif_sasa', 'std_motif_sasa'])
                
                # Add chain-specific column names
                for col in chain_cols:
                    chain_id = col.replace('chain_', '').replace('_sasa', '')
                    new_cols.extend([f'mean_chain_{chain_id}_sasa', f'std_chain_{chain_id}_sasa'])
                
                sasa_cluster_stats.columns = new_cols
                
                # Merge with summary
                summary_with_rmsd = summary_with_rmsd.merge(
                    sasa_cluster_stats, on='cluster', how='left'
                )
                log_status("✓ ΔSASA stats added to cluster summary")
            
            # Add motif match scores if available
            if hasattr(clust, 'motif_match_scores') and clust.motif_match_scores is not None:
                try:
                    match_scores = clust.motif_match_scores
                    if len(match_scores) > 0 and any(s is not None for s in match_scores):
                        # Compute average match score per cluster
                        match_score_by_cluster = {}
                        for cid in summary_with_rmsd['cluster'].values:
                            if cid == -1:
                                match_score_by_cluster[cid] = 0.0
                            else:
                                cluster_indices = [i for i, l in enumerate(clust.labels) if l == cid]
                                cluster_scores = [match_scores[i] for i in cluster_indices if match_scores[i] is not None]
                                if len(cluster_scores) > 0:
                                    match_score_by_cluster[cid] = np.mean(cluster_scores)
                                else:
                                    match_score_by_cluster[cid] = 0.0
                        
                        summary_with_rmsd['avg_motif_match'] = summary_with_rmsd['cluster'].map(match_score_by_cluster)
                        log_status("✓ Motif match scores added to cluster summary")
                except Exception as e:
                    log_status(f"⚠️ Could not compute motif match scores: {str(e)}")
            
            # Add cluster stability scores if HDBSCAN clustering was used
            if hasattr(clust, 'hdbscan_clusterer') and clust.hdbscan_clusterer is not None:
                try:
                    from visualization.interactive import create_cluster_stability_plot
                    # Get stability scores (returns figure, but we extract the scores)
                    persistence_scores = clust.hdbscan_clusterer.cluster_persistence_
                    cluster_ids_for_stability = [idx for idx in summary_with_rmsd['cluster'].values if idx != -1]
                    
                    # Check if we need fallback calculation
                    if np.all(persistence_scores == 1.0) or np.any(np.isinf(persistence_scores)):
                        # Use hybrid fallback to get stability scores
                        tree_df = clust.hdbscan_clusterer.condensed_tree_.to_pandas()
                        tree_df_finite = tree_df[np.isfinite(tree_df['lambda_val'])].copy()
                        
                        stability_scores = {}
                        for cid in cluster_ids_for_stability:
                            cluster_mask = np.array(clust.labels) == cid
                            cluster_indices = np.where(cluster_mask)[0]
                            
                            if len(cluster_indices) > 1:
                                # Compute median distance
                                cluster_submatrix = clust.distance_matrix[np.ix_(cluster_indices, cluster_indices)]
                                upper_tri_indices = np.triu_indices_from(cluster_submatrix, k=1)
                                cluster_dists = cluster_submatrix[upper_tri_indices]
                                median_dist = max(np.median(cluster_dists), 1e-6)
                                
                                # Get lambda range
                                cluster_point_rows = tree_df_finite[tree_df_finite['child'].isin(cluster_indices)]
                                if len(cluster_point_rows) > 0:
                                    lambda_vals = cluster_point_rows['lambda_val'].values
                                    lambda_vals_finite = lambda_vals[np.isfinite(lambda_vals)]
                                    if len(lambda_vals_finite) > 1:
                                        lambda_range = max(lambda_vals_finite.max() - lambda_vals_finite.min(), 1e-6)
                                        stability = (lambda_range * len(cluster_indices)) / median_dist
                                    else:
                                        stability = float(len(cluster_indices)) / median_dist
                                else:
                                    stability = float(len(cluster_indices)) / median_dist
                            else:
                                stability = 0.0
                            
                            stability_scores[cid] = stability
                    else:
                        # Use primary persistence scores
                        stability_scores = {cid: score for cid, score in zip(cluster_ids_for_stability, persistence_scores)}
                    
                    # Add to summary
                    summary_with_rmsd['stability_score'] = summary_with_rmsd['cluster'].map(
                        lambda c: stability_scores.get(c, 0.0) if c != -1 else 0.0
                    )
                    log_status("✓ Stability scores added to cluster summary")
                except Exception as e:
                    log_status(f"⚠️ Could not compute stability scores: {str(e)}")
            
            cluster_summary.set(summary_with_rmsd)
            log_status("✓ Cluster summary computed")
            log_status("")
            
            # Update cluster filter dropdowns
            cluster_ids = [str(c) for c in sorted(set(clust.labels))]
            cluster_choices = {}
            for cid in cluster_ids:
                n_structs = sum(1 for l in clust.labels if l == int(cid))
                if int(cid) == -1:
                    cluster_choices[cid] = f"Noise ({n_structs} structures)"
                else:
                    cluster_choices[cid] = f"Cluster {cid} ({n_structs} structures)"
            # Add "All Clusters" option at the end
            cluster_choices["all"] = "All Clusters (use with caution!)"
            
            # Build intra-cluster selector choices (exclude noise)
            intra_cluster_choices = {}
            for cid in cluster_ids:
                if int(cid) >= 0:  # Exclude noise (-1)
                    n_structs = sum(1 for l in clust.labels if l == int(cid))
                    intra_cluster_choices[cid] = f"Cluster {cid} ({n_structs} structures)"
            
            # Default to Cluster 1 (second cluster, usually smaller)
            # Fall back to cluster 0 if only 1 cluster, or "all" if no clusters
            if len(cluster_ids) > 1:
                default_cluster = cluster_ids[1]  # Cluster 1
            elif len(cluster_ids) == 1:
                default_cluster = cluster_ids[0]  # Only cluster 0
            else:
                default_cluster = "all"  # No clusters
            
            # Default intra-cluster selection (first non-noise cluster)
            default_intra = next((cid for cid in cluster_ids if int(cid) >= 0), "0")
            
            ui.update_select("filter_cluster", choices=cluster_choices, selected=default_cluster)
            ui.update_select("viewer_cluster", choices=cluster_choices, selected=default_cluster)
            ui.update_select("cluster_select_hotspot", choices=cluster_ids)
            ui.update_select("contact_heatmap_cluster", choices=cluster_ids, selected=default_intra)
            ui.update_select("jaccard_intra_cluster", choices=intra_cluster_choices, selected=default_intra)
            ui.update_select("rmsd_intra_cluster", choices=intra_cluster_choices, selected=default_intra)
            
            # Update export cluster dropdown (only specific clusters, no "all" option)
            export_choices = {}
            for cid in cluster_ids:
                n_structs = sum(1 for l in clust.labels if l == int(cid))
                export_choices[cid] = f"Cluster {cid} ({n_structs} structures)"
            ui.update_select("export_cluster", choices=export_choices, selected=default_cluster if default_cluster != "all" else cluster_ids[0] if cluster_ids else "")
            
            # Initialize structure dropdown for 3D viewer with default cluster
            structure_names = [os.path.basename(f) for f in clust.pdb_files]
            if default_cluster != "all":
                # Filter by default cluster
                cluster_id = int(default_cluster)
                filtered_indices = [i for i, label in enumerate(clust.labels) if label == cluster_id]
                filtered_names = [structure_names[i] for i in filtered_indices]
                ui.update_select("structure_select", choices=filtered_names)
            else:
                # Show all structures
                ui.update_select("structure_select", choices=structure_names)
            
            # Set analysis complete flag
            analysis_complete.set(True)
            log_status("="*50)
            log_status("✅ Analysis Complete!")
            log_status(f"Results: {n_clusters} clusters, {n_noise} noise points")
            log_status("="*50)
            
            ui.notification_show("✅ Analysis complete!", type="message")
            
        except Exception as e:
            analysis_complete.set(False)
            ui.notification_show(f"Error: {str(e)}", type="error")
            log_status("❌ ERROR: Analysis failed")
            log_status(f"Error details: {str(e)}")
            import traceback
            traceback.print_exc()
    
    @reactive.effect
    @reactive.event(input.recluster)
    def recluster_only():
        """
        Re-run clustering with different parameters using cached distance matrices
        This is FAST because it skips all expensive computations!
        """
        global clust
        
        if not analysis_complete.get():
            ui.notification_show("Run full analysis first!", type="warning")
            return
        
        if clust is None:
            ui.notification_show("No data available - run analysis first", type="error")
            return
        
        # Get new parameters
        min_cluster_size = input.min_cluster_size()
        min_samples_hdb = input.min_samples_hdb()
        filter_dupes = input.filter_duplicates()
        dup_threshold = input.duplicate_threshold()
        
        log_status("")
        log_status("="*50)
        log_status("🔄 Re-Clustering with new parameters...")
        log_status(f"New min_cluster_size: {min_cluster_size}, New min_samples: {min_samples_hdb}")
        log_status(f"Filter duplicates: {filter_dupes}, Threshold: {dup_threshold}")
        log_status("="*50)
        
        try:
            # Re-cluster using CACHED Jaccard matrix (instant!)
            clust.distance_matrix = clust.jaccard_distance_matrix
            clust.metric_name = "Jaccard"
            clust.cluster_hdbscan(
                min_cluster_size=min_cluster_size, 
                min_samples=min_samples_hdb,
                filter_duplicates=filter_dupes,
                duplicate_threshold=dup_threshold
            )
            
            n_clusters = len(set(clust.labels)) - (1 if -1 in clust.labels else 0)
            n_noise = sum(1 for l in clust.labels if l == -1)
            log_status(f"✓ Re-clustering complete: {n_clusters} clusters, {n_noise} noise points")
            
            # Recompute cluster statistics with cached RMSD
            rmsd_mat = rmsd_matrix.get()
            if rmsd_mat is not None:
                rmsd_qc = compute_rmsd_qc(clust.labels, rmsd_mat)
                rmsd_stats.set(rmsd_qc)
                log_status("✓ RMSD stats updated")
            
            # Recompute interface stats (fast - just regrouping)
            stats_df = clust.get_interface_stats()
            
            # Mark representatives
            if rmsd_qc is not None:
                medoid_map = dict(zip(rmsd_qc['cluster'], rmsd_qc['medoid_idx']))
                stats_df['structure_idx'] = range(len(stats_df))
                def is_representative(row):
                    cluster = row['cluster']
                    if cluster == -1:
                        return False
                    medoid_idx = medoid_map.get(cluster, -1)
                    return row['structure_idx'] == medoid_idx
                stats_df['representative'] = stats_df.apply(is_representative, axis=1)
                stats_df = stats_df.drop(columns=['structure_idx'])
            
            # Merge SASA data
            sasa_df_current = sasa_data.get()
            if sasa_df_current is not None:
                sasa_cols_to_merge = [col for col in sasa_df_current.columns if col != 'structure']
                stats_df = stats_df.merge(
                    sasa_df_current[['structure'] + sasa_cols_to_merge],
                    on='structure',
                    how='left'
                )
            
            # Store updated stats
            interface_stats.set(stats_df)
            
            # Recompute cluster summary
            log_status("Recomputing cluster summary...")
            summary = stats_df.groupby('cluster').agg({
                'n_protein_residues': ['count', 'mean', 'std'],
                'protein_residues': lambda x: len(set().union(*[set(r.split(',')) for r in x if r]))
            }).reset_index()
            
            summary.columns = ['cluster', 'count', 'mean_contacts', 'std_contacts', 'total_unique_residues']
            
            # Add RMSD stats
            summary_with_rmsd = summary.merge(rmsd_qc[['cluster', 'mean_rmsd', 'std_rmsd', 'medoid_idx']], 
                                               on='cluster', how='left')
            
            # Add SASA stats
            if sasa_df_current is not None:
                sasa_with_labels = sasa_df_current.copy()
                sasa_with_labels['cluster'] = clust.labels
                agg_dict = {'total_delta_sasa': ['mean', 'std', 'max']}
                if 'motif_delta_sasa' in sasa_with_labels.columns:
                    agg_dict['motif_delta_sasa'] = ['mean', 'std']
                chain_cols = [col for col in sasa_with_labels.columns if col.startswith('chain_') and col.endswith('_sasa')]
                for col in chain_cols:
                    agg_dict[col] = ['mean', 'std']
                sasa_cluster_stats = sasa_with_labels.groupby('cluster').agg(agg_dict).reset_index()
                new_cols = ['cluster', 'mean_delta_sasa', 'std_delta_sasa', 'max_delta_sasa']
                if 'motif_delta_sasa' in agg_dict:
                    new_cols.extend(['mean_motif_sasa', 'std_motif_sasa'])
                for col in chain_cols:
                    chain_id = col.replace('chain_', '').replace('_sasa', '')
                    new_cols.extend([f'mean_chain_{chain_id}_sasa', f'std_chain_{chain_id}_sasa'])
                sasa_cluster_stats.columns = new_cols
                summary_with_rmsd = summary_with_rmsd.merge(sasa_cluster_stats, on='cluster', how='left')
            
            cluster_summary.set(summary_with_rmsd)
            
            # Update selectors
            cluster_ids = [str(c) for c in sorted(set(clust.labels)) if c != -1]
            cluster_choices = {}
            for cid in cluster_ids:
                n_structs = sum(1 for l in clust.labels if l == int(cid))
                cluster_choices[cid] = f"Cluster {cid} ({n_structs} structures)"
            # Add "All Clusters" option at the end
            cluster_choices["all"] = "All Clusters (use with caution!)"
            
            # Default to Cluster 1 (second cluster, usually smaller)
            if len(cluster_ids) > 1:
                default_cluster = cluster_ids[1]  # Cluster 1
            elif len(cluster_ids) == 1:
                default_cluster = cluster_ids[0]  # Only cluster 0
            else:
                default_cluster = "all"  # No clusters
            
            ui.update_select("filter_cluster", choices=cluster_choices, selected=default_cluster)
            ui.update_select("viewer_cluster", choices=cluster_choices, selected=default_cluster)
            ui.update_select("cluster_select_hotspot", choices=cluster_ids)
            ui.update_select("contact_heatmap_cluster", choices=cluster_ids, selected=cluster_ids[0] if cluster_ids else "")
            
            # Update export cluster dropdown
            export_choices = {}
            for cid in cluster_ids:
                n_structs = sum(1 for l in clust.labels if l == int(cid))
                export_choices[cid] = f"Cluster {cid} ({n_structs} structures)"
            ui.update_select("export_cluster", choices=export_choices, selected=default_cluster if default_cluster != "all" else cluster_ids[0] if cluster_ids else "")
            
            # Update structure dropdown for 3D viewer
            structure_names = [os.path.basename(f) for f in clust.pdb_files]
            if default_cluster != "all":
                cluster_id = int(default_cluster)
                filtered_indices = [i for i, label in enumerate(clust.labels) if label == cluster_id]
                filtered_names = [structure_names[i] for i in filtered_indices]
                ui.update_select("structure_select", choices=filtered_names)
            else:
                ui.update_select("structure_select", choices=structure_names)
            
            log_status("✓ All statistics updated")
            log_status("="*50)
            log_status(f"✅ Re-Clustering Complete!")
            log_status(f"Results: {n_clusters} clusters, {n_noise} noise points")
            log_status("="*50)
            
            ui.notification_show(f"✅ Re-clustered: {n_clusters} clusters!", type="message")
            
        except Exception as e:
            ui.notification_show(f"Re-clustering error: {str(e)}", type="error")
            log_status(f"❌ Re-clustering failed: {str(e)}")
            import traceback
            traceback.print_exc()
    
    @reactive.effect
    @reactive.event(input.viewer_cluster, input.structure_type)
    def update_structure_choices():
        """Update structure dropdown when cluster or structure type filter changes in 3D viewer"""
        if clust is None:
            return
        
        selected_cluster = input.viewer_cluster()
        structure_type = input.structure_type()
        if not selected_cluster:
            return
        
        structure_names = [os.path.basename(f) for f in clust.pdb_files]
        
        # Get medoid indices for filtering
        medoid_map = {}
        stats_df = interface_stats.get()
        if stats_df is not None and 'representative' in stats_df.columns:
            # Find medoid for each cluster from representative column
            for cluster_id in set(clust.labels):
                if cluster_id == -1:
                    continue
                cluster_stats = stats_df[stats_df['cluster'] == cluster_id]
                rep_stats = cluster_stats[cluster_stats['representative'] == True]
                if len(rep_stats) > 0:
                    # Get the index of the representative
                    rep_structure = rep_stats.iloc[0]['structure']
                    for i, f in enumerate(clust.pdb_files):
                        if os.path.basename(f) == rep_structure:
                            medoid_map[cluster_id] = i
                            break
        
        if selected_cluster == "all":
            # Show all structures or filter by type
            if structure_type == "all":
                filtered_indices = list(range(len(structure_names)))
            elif structure_type == "medoid" and medoid_map:
                # Show only medoids from all clusters
                filtered_indices = list(medoid_map.values())
            elif structure_type == "non_medoid" and medoid_map:
                # Show only non-medoids from all clusters
                medoid_set = set(medoid_map.values())
                filtered_indices = [i for i in range(len(structure_names)) if i not in medoid_set]
            else:
                filtered_indices = list(range(len(structure_names)))
        else:
            # Filter by cluster
            cluster_id = int(selected_cluster)
            cluster_indices = [i for i, label in enumerate(clust.labels) if label == cluster_id]
            
            if structure_type == "all":
                filtered_indices = cluster_indices
            elif structure_type == "medoid" and cluster_id in medoid_map:
                # Show only medoid of this cluster
                medoid_idx = medoid_map[cluster_id]
                filtered_indices = [medoid_idx] if medoid_idx in cluster_indices else []
            elif structure_type == "non_medoid" and cluster_id in medoid_map:
                # Show only non-medoids of this cluster
                medoid_idx = medoid_map[cluster_id]
                filtered_indices = [i for i in cluster_indices if i != medoid_idx]
            else:
                filtered_indices = cluster_indices
        
        filtered_names = [structure_names[i] for i in filtered_indices]
        ui.update_select("structure_select", choices=filtered_names)
    
    @output
    @render.text
    def status():
        # Display accumulated terminal-style log
        log_messages = status_log.get()
        
        if not log_messages:
            return "No analysis run yet.\nEnter directory and click 'Run Analysis'."
        
        return "\n".join(log_messages)
    
    # ========== OVERVIEW TAB ==========
    
    @output
    @render.ui
    def cluster_summary_cards():
        """Summary cards showing cluster counts"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see summary", style="padding: 20px; text-align: center;")
        
        n_clusters = len(set(clust.labels)) - (1 if -1 in clust.labels else 0)
        n_noise = sum(1 for l in clust.labels if l == -1)
        n_total = len(clust.pdb_files)
        
        return ui.TagList(
            ui.div(
                ui.h2(f"{n_clusters}", class_="text-primary", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Clusters Found", style="margin: 5px 0;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            ),
            ui.div(
                ui.h2(f"{n_total}", class_="text-info", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Total Structures", style="margin: 5px 0;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            ),
            ui.div(
                ui.h2(f"{n_noise}", class_="text-warning", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Noise Points", style="margin: 5px 0;"),
                class_="text-center",
                style="padding: 10px;"
            )
        )
    
    @output
    @render.ui
    def structure_quality_cards():
        """Quality metric cards for structure metrics (DBCV, RMSD)"""
        global clust  # Make sure we're accessing the global clust
        
        rmsd_mat = rmsd_matrix.get()
        
        print(f"\n=== DEBUG structure_quality_cards() ===", flush=True)
        print(f"analysis_complete: {analysis_complete.get()}", flush=True)
        print(f"rmsd_mat is None: {rmsd_mat is None}", flush=True)
        print(f"clust is None: {clust is None}", flush=True)
        sys.stdout.flush()
        
        if not analysis_complete.get() or rmsd_mat is None or clust is None:
            print("DEBUG: Early return - not ready", flush=True)
            return ui.p("Run analysis to see metrics", style="padding: 20px; text-align: center;")
        
        # Get clustering quality score from HDBSCAN
        # We use median GLOSH (Global-Local Outlier Score from Hierarchies)
        # GLOSH measures how outlier-like each point is
        # Lower scores = better clustering (points are more inlier-like)
        quality_score = 0.0
        metric_label = "Cluster Quality"
        
        print(f"DEBUG: clust object type: {type(clust)}", flush=True)
        print(f"DEBUG: clust is None? {clust is None}", flush=True)
        sys.stdout.flush()
        
        try:
            if clust is None:
                print("DEBUG: clust is None!", flush=True)
            elif not hasattr(clust, 'hdbscan_clusterer'):
                print("DEBUG: clust has no hdbscan_clusterer attribute", flush=True)
            elif clust.hdbscan_clusterer is None:
                print("DEBUG: clust.hdbscan_clusterer is None", flush=True)
            elif not hasattr(clust.hdbscan_clusterer, 'outlier_scores_'):
                print("DEBUG: hdbscan_clusterer has no outlier_scores_ attribute", flush=True)
            else:
                outlier_scores = clust.hdbscan_clusterer.outlier_scores_
                print(f"DEBUG: Raw GLOSH scores shape: {outlier_scores.shape}", flush=True)
                finite_scores = outlier_scores[np.isfinite(outlier_scores)]
                print(f"DEBUG: Finite GLOSH scores: min={finite_scores.min():.3f}, max={finite_scores.max():.3f}", flush=True)
                if len(finite_scores) > 0:
                    quality_score = np.median(finite_scores)
                    metric_label = "Median GLOSH"
                    print(f"DEBUG: Computed quality_score (median GLOSH) = {quality_score}", flush=True)
                else:
                    print("DEBUG: No finite GLOSH values!", flush=True)
            sys.stdout.flush()
        except Exception as e:
            print(f"DEBUG: Error getting clustering quality: {e}", flush=True)
            import traceback
            traceback.print_exc()
            sys.stdout.flush()
            quality_score = 0.0
        
        # Get sequential RMSD from cluster summary
        rmsd_qc_df = rmsd_stats.get()
        if rmsd_qc_df is not None and 'mean_rmsd_sequential' in rmsd_qc_df.columns:
            # Average sequential RMSD across all real clusters (exclude noise)
            real_clusters = rmsd_qc_df[rmsd_qc_df['cluster'] != -1]
            if len(real_clusters) > 0:
                mean_rmsd_seq = real_clusters['mean_rmsd_sequential'].mean()
            else:
                mean_rmsd_seq = rmsd_mat.mean()
        else:
            mean_rmsd_seq = rmsd_mat.mean()
        
        max_rmsd = rmsd_mat.max()
        
        # 3 structure quality cards
        return ui.TagList(
            ui.div(
                ui.h2(f"{quality_score:.3f}", class_="text-success", style="font-size: 3em; margin: 10px 0;"),
                ui.p(metric_label, style="margin: 5px 0;"),
                #ui.p("Lower = better clustering", style="margin: 0; font-size: 0.8em; color: #666;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            ),
            ui.div(
                ui.h2(f"{mean_rmsd_seq:.1f} Å", class_="text-info", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Mean RMSD", style="margin: 5px 0; font-size: 0.9em;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            ),
            ui.div(
                ui.h2(f"{max_rmsd:.1f} Å", class_="text-danger", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Max RMSD", style="margin: 5px 0;"),
                class_="text-center",
                style="padding: 10px;"
            )
        )
    
    @output
    @render.ui
    def interface_quality_cards():
        """Quality metric cards for interface metrics (ΔSASA)"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see metrics", style="padding: 20px; text-align: center;")
        
        # Get ΔSASA statistics
        sasa_df = sasa_data.get()
        if sasa_df is not None and 'total_delta_sasa' in sasa_df.columns:
            mean_sasa = sasa_df['total_delta_sasa'].mean()
            max_sasa = sasa_df['total_delta_sasa'].max()
            # Get motif SASA on protein (should be constant, take first non-null value)
            if 'motif_sasa_protein' in sasa_df.columns:
                motif_sasa_protein = sasa_df['motif_sasa_protein'].dropna().iloc[0] if len(sasa_df['motif_sasa_protein'].dropna()) > 0 else None
            else:
                motif_sasa_protein = None
        else:
            # No SASA data
            return ui.p("ΔSASA not calculated", style="padding: 20px; text-align: center; color: #999;")
        
        # 3 interface quality cards
        cards = [
            ui.div(
                ui.h2(f"{mean_sasa:.0f} Å²", class_="text-warning", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Mean Buried ΔSASA", style="margin: 5px 0; font-size: 0.9em;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            ),
            ui.div(
                ui.h2(f"{max_sasa:.0f} Å²", class_="text-warning", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Max Buried ΔSASA", style="margin: 5px 0; font-size: 0.9em;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            )
        ]
        
        # Add motif SASA if available
        if motif_sasa_protein is not None:
            cards.append(
                ui.div(
                    ui.h2(f"{motif_sasa_protein:.0f} Å²", class_="text-warning", style="font-size: 3em; margin: 10px 0;"),
                    ui.p("Motif SASA", style="margin: 5px 0;"),
                    class_="text-center",
                    style="padding: 10px;"
                )
            )
        else:
            # Placeholder for alignment
            cards.append(
                ui.div(
                    ui.h2("—", class_="text-muted", style="font-size: 3em; margin: 10px 0;"),
                    ui.p("Motif SASA", style="margin: 5px 0; color: #999;"),
                    class_="text-center",
                    style="padding: 10px;"
                )
            )
        
        return ui.TagList(*cards)
    
    @output
    @render.data_frame
    def cluster_summary_table():
        summary_df = cluster_summary.get()
        if summary_df is None or summary_df.empty:
            return pd.DataFrame()
        
        # Add medoid structure name if we have it
        if 'medoid_idx' in summary_df.columns and clust is not None:
            summary_df['medoid_structure'] = summary_df['medoid_idx'].apply(
                lambda idx: os.path.basename(clust.pdb_files[idx]) if idx >= 0 and idx < len(clust.pdb_files) else 'N/A'
            )
        
        # Show key columns including intra-cluster RMSD metrics, stability, motif match, and medoid info
        display_cols = ['cluster', 'binding_mode', 'n_structures', 'stability_score', 'avg_motif_match', 'n_consensus', 
                       'mean_rmsd', 'std_rmsd', 'median_rmsd', 'max_rmsd', 'min_rmsd',
                       'mean_rmsd_sequential', 'medoid_avg_rmsd', 'medoid_structure']
        
        # Add ΔSASA columns if available
        if 'mean_delta_sasa' in summary_df.columns:
            display_cols.extend(['mean_delta_sasa', 'std_delta_sasa', 'max_delta_sasa'])
        if 'mean_motif_sasa' in summary_df.columns:
            display_cols.extend(['mean_motif_sasa', 'std_motif_sasa'])
        
        # Add per-chain ΔSASA columns if available
        chain_sasa_cols = [col for col in summary_df.columns if col.startswith('mean_chain_') and col.endswith('_sasa')]
        for col in chain_sasa_cols:
            display_cols.append(col)
            # Also add std columns
            std_col = col.replace('mean_', 'std_')
            if std_col in summary_df.columns:
                display_cols.append(std_col)
        
        available_cols = [c for c in display_cols if c in summary_df.columns]
        display_df = summary_df[available_cols].copy()

        # Apply cluster search/filter if provided
        try:
            search_val = input.cluster_search() if hasattr(input, 'cluster_search') else None
            if search_val:
                s = str(search_val).strip()
                if s:
                    # Case-insensitive search across all displayed columns
                    mask = display_df.astype(str).apply(lambda row: row.str.contains(s, case=False, na=False).any(), axis=1)
                    display_df = display_df[mask]
        except Exception:
            # If input isn't available in this context, skip filtering
            pass
        
        # Rename for clarity with proper units
        rename_map = {}
        if 'stability_score' in display_df.columns:
            rename_map['stability_score'] = 'Stability Score'
        if 'avg_motif_match' in display_df.columns:
            rename_map['avg_motif_match'] = 'Avg Motif Match (%)'
        if 'mean_rmsd' in display_df.columns:
            rename_map['mean_rmsd'] = 'Intra-Cluster Mean RMSD (Å)'
        if 'std_rmsd' in display_df.columns:
            rename_map['std_rmsd'] = 'Intra-Cluster Std RMSD (Å)'
        if 'median_rmsd' in display_df.columns:
            rename_map['median_rmsd'] = 'Median RMSD (Å)'
        if 'max_rmsd' in display_df.columns:
            rename_map['max_rmsd'] = 'Max RMSD (Å)'
        if 'min_rmsd' in display_df.columns:
            rename_map['min_rmsd'] = 'Min RMSD (Å)'
        if 'mean_rmsd_sequential' in display_df.columns:
            rename_map['mean_rmsd_sequential'] = 'Sequential RMSD (Å)'
        if 'medoid_avg_rmsd' in display_df.columns:
            rename_map['medoid_avg_rmsd'] = 'Medoid Avg RMSD (Å)'
        if 'mean_delta_sasa' in display_df.columns:
            rename_map['mean_delta_sasa'] = 'Mean Buried ΔSASA (Å²)'
        if 'std_delta_sasa' in display_df.columns:
            rename_map['std_delta_sasa'] = 'Std Buried ΔSASA (Å²)'
        if 'max_delta_sasa' in display_df.columns:
            rename_map['max_delta_sasa'] = 'Max Buried ΔSASA (Å²)'
        if 'mean_motif_sasa' in display_df.columns:
            rename_map['mean_motif_sasa'] = 'Mean Motif ΔSASA (Å²)'
        if 'std_motif_sasa' in display_df.columns:
            rename_map['std_motif_sasa'] = 'Std Motif ΔSASA (Å²)'
        
        # Rename per-chain ΔSASA columns
        for col in display_df.columns:
            if col.startswith('mean_chain_') and col.endswith('_sasa'):
                chain_id = col.replace('mean_chain_', '').replace('_sasa', '')
                rename_map[col] = f'Mean Chain {chain_id} ΔSASA (Å²)'
            elif col.startswith('std_chain_') and col.endswith('_sasa'):
                chain_id = col.replace('std_chain_', '').replace('_sasa', '')
                rename_map[col] = f'Std Chain {chain_id} ΔSASA (Å²)'
        
        if rename_map:
            display_df = display_df.rename(columns=rename_map)
        
        # Format numeric columns to 4 significant figures
        numeric_cols = ['Avg Motif Match (%)', 'Stability Score',
                       'Intra-Cluster Mean RMSD (Å)', 'Intra-Cluster Std RMSD (Å)', 
                       'Median RMSD (Å)', 'Max RMSD (Å)', 'Min RMSD (Å)',
                       'Sequential RMSD (Å)', 'Medoid Avg RMSD (Å)',
                       'Mean Buried ΔSASA (Å²)', 'Std Buried ΔSASA (Å²)', 'Max Buried ΔSASA (Å²)',
                       'Mean Motif ΔSASA (Å²)', 'Std Motif ΔSASA (Å²)']
        
        # Add chain ΔSASA columns to numeric formatting
        for col in display_df.columns:
            if 'Chain' in col and 'ΔSASA' in col:
                numeric_cols.append(col)
        
        for col in numeric_cols:
            if col in display_df.columns:
                display_df[col] = display_df[col].apply(lambda x: f"{x:.4g}" if pd.notna(x) else x)
        
        return render.DataGrid(display_df, width="100%", height=400)
    
    @output
    @render.ui

    def cluster_size_plot():
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see cluster sizes", style="padding: 20px; text-align: center;")
        
        try:
            fig = create_cluster_size_distribution(
                clust.labels,
                duplicate_groups=getattr(clust, 'duplicate_groups', None),
                rescued_cluster_ids=getattr(clust, 'rescued_cluster_ids', None)
            )
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id='cluster_size_dist')
            wrapped_html = f'<div style="width: 100%; height: 400px;">{html}</div>'
            return ui.HTML(wrapped_html)
        except Exception as e:
            return ui.HTML(f"<pre style='color: red;'>Error: {e}</pre>")
    
    @output
    @render.ui
    def cluster_stability_plot():
        """Display HDBSCAN cluster stability (persistence) scores"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see cluster stability", style="padding: 20px; text-align: center;")
        
        # Check if HDBSCAN clusterer is available
        if not hasattr(clust, 'hdbscan_clusterer') or clust.hdbscan_clusterer is None:
            return ui.p("Cluster stability requires HDBSCAN clustering", 
                       style="padding: 20px; text-align: center; color: #999;")
        
        try:
            # Pass distance matrix for fallback calculations
            distance_mat = clust.distance_matrix if hasattr(clust, 'distance_matrix') else None
            fig = create_cluster_stability_plot(
                clust.hdbscan_clusterer, 
                clust.labels, 
                distance_mat,
                duplicate_groups=getattr(clust, 'duplicate_groups', None),
                rescued_cluster_ids=getattr(clust, 'rescued_cluster_ids', None)
            )
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id='cluster_stability')
            wrapped_html = f'<div style="width: 100%; height: 400px;">{html}</div>'
            return ui.HTML(wrapped_html)
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def cluster_stability_info():
        """Display stability plot warning text below the plot"""
        if not analysis_complete.get() or clust is None:
            return None
        
        if not hasattr(clust, 'hdbscan_clusterer') or clust.hdbscan_clusterer is None:
            return None
        
        try:
            # Recreate the plot to get the metadata (lightweight operation)
            distance_mat = clust.distance_matrix if hasattr(clust, 'distance_matrix') else None
            fig = create_cluster_stability_plot(
                clust.hdbscan_clusterer, 
                clust.labels, 
                distance_mat,
                duplicate_groups=getattr(clust, 'duplicate_groups', None),
                rescued_cluster_ids=getattr(clust, 'rescued_cluster_ids', None)
            )
            warning_text = fig.layout.meta.get('warning_text', '')
            
            if warning_text:
                # Style the info box with blue accent
                return ui.div(
                    ui.p(warning_text, style="margin: 0; white-space: pre-line;"),
                    style="""
                        padding: 10px; 
                        margin: 10px 0; 
                        background-color: #e7f3ff; 
                        border-left: 4px solid #2196F3; 
                        font-size: 0.9em;
                        color: #333;
                        border-radius: 4px;
                    """
                )
            return None
        except Exception:
            return None
    
    @output
    @render.ui
    def motif_match_plot():
        """Display motif match score bar chart"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see motif match scores", style="padding: 20px; text-align: center;")
        
        # Check if motif match scores are available
        if not hasattr(clust, 'motif_match_scores') or clust.motif_match_scores is None:
            return ui.p("Motif match scores not available (motif residues not specified)", 
                       style="padding: 20px; text-align: center; color: #999;")
        
        try:
            from visualization.interactive import create_motif_match_score_plot
            fig = create_motif_match_score_plot(clust.labels, clust.motif_match_scores, clust.pdb_files)
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id='motif_match_score')
            wrapped_html = f'<div style="width: 100%; height: 400px;">{html}</div>'
            return ui.HTML(wrapped_html)
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def cluster_diversity_plot():
        """Display pose diversity within clusters (RMSD from medoid)"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see cluster diversity", style="padding: 20px; text-align: center;")
        
        rmsd_mat = rmsd_matrix.get()
        if rmsd_mat is None:
            return ui.p("RMSD matrix not available", style="padding: 20px; text-align: center; color: #999;")
        
        try:
            fig = create_cluster_diversity_plot(
                clust.labels, 
                rmsd_mat, 
                clust.pdb_files,
                sasa_data=sasa_data.get()
            )
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id='cluster_diversity')
            # Dynamic height based on number of clusters
            n_clusters = len(set(clust.labels)) - (1 if -1 in clust.labels else 0)
            height = max(400, n_clusters * 40)
            wrapped_html = f'<div style="width: 100%; height: {height}px;">{html}</div>'
            return ui.HTML(wrapped_html)
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def hdbscan_condensed_tree_plot():
        """Display HDBSCAN condensed tree showing cluster hierarchy and stability"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see HDBSCAN condensed tree", style="padding: 20px; text-align: center;")
        
        # Check if HDBSCAN clusterer is available
        if not hasattr(clust, 'hdbscan_clusterer') or clust.hdbscan_clusterer is None:
            return ui.p("HDBSCAN tree not available (may be using different clustering method)", 
                       style="padding: 20px; text-align: center; color: #999;")
        
        try:
            # Pass the PDBClusterer object - it has the cached tree from clustering time
            html = create_hdbscan_condensed_tree(
                clust,  # Pass the whole clusterer object
                labels=clust.labels,
                labels_for_tree=None,
                pdb_files=clust.pdb_files
            )
            return ui.HTML(html)
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def hdbscan_linkage_tree_plot():
        """Display HDBSCAN single linkage tree (dendrogram)"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see HDBSCAN linkage tree", style="padding: 20px; text-align: center;")
        
        # Check if HDBSCAN clusterer is available
        if not hasattr(clust, 'hdbscan_clusterer') or clust.hdbscan_clusterer is None:
            return ui.p("HDBSCAN tree not available (may be using different clustering method)", 
                       style="padding: 20px; text-align: center; color: #999;")
        
        try:
            html = create_hdbscan_tree(clust.hdbscan_clusterer, clust.labels, clust.pdb_files)
            return ui.HTML(html)
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def jaccard_heatmap():
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see Jaccard distance matrix", style="padding: 20px; text-align: center;")
        
        try:
            jmat = clust.jaccard_distance_matrix
            labels = clust.labels
            level = input.jaccard_heatmap_level()
            show_dendrogram = input.jaccard_dendrogram()
            
            if level == "cluster":
                # Cluster-level heatmap
                cluster_matrix, cluster_ids, cluster_sizes = compute_cluster_distance_matrix(jmat, labels)
                if show_dendrogram:
                    # Create cluster labels for dendrogram
                    cluster_labels = [f"Cluster {cid} (n={cluster_sizes[cid]})" for cid in cluster_ids]
                    # Create a temporary array with cluster labels as "pdb_files" for the dendrogram
                    fig = create_distance_heatmap_with_dendrogram(cluster_matrix, labels=None, metric_name='Jaccard Distance', pdb_files=cluster_labels)
                else:
                    fig = create_cluster_distance_heatmap(cluster_matrix, cluster_ids, cluster_sizes, metric_name='Jaccard Distance')
                return ui.HTML(fig.to_html(full_html=False, include_plotlyjs='require'))
            
            elif level == "intra":
                # Single cluster intra-cluster heatmap
                selected_cluster = input.jaccard_intra_cluster()
                if not selected_cluster:
                    return ui.p("Select a cluster to view intra-cluster distances", style="padding: 20px; text-align: center;")
                
                cluster_id = int(selected_cluster)
                fig = create_single_intra_cluster_heatmap(jmat, labels, cluster_id, clust.pdb_files, metric_name='Jaccard Distance')
                return ui.HTML(fig.to_html(full_html=False, include_plotlyjs='require'))
            
            else:  # pose level
                # Limit to first 500 structures to prevent browser crash
                max_display = 500
                if jmat.shape[0] > max_display:
                    jmat_display = jmat[:max_display, :max_display]
                    labels_display = labels[:max_display]
                    pdb_files_display = clust.pdb_files[:max_display] if clust.pdb_files else None
                    note = f"<div style='padding: 10px; background: #fff3cd; margin-bottom: 10px;'>⚠️ Showing first {max_display} of {clust.jaccard_distance_matrix.shape[0]} structures</div>"
                else:
                    jmat_display = jmat
                    labels_display = labels
                    pdb_files_display = clust.pdb_files
                    note = ""
                
                if show_dendrogram:
                    fig = create_distance_heatmap_with_dendrogram(jmat_display, labels_display, 
                                                                  metric_name='Jaccard Distance',
                                                                  pdb_files=pdb_files_display)
                else:
                    fig = create_distance_heatmap(jmat_display, labels_display, metric_name='Jaccard Distance')
                
                return ui.HTML(note + fig.to_html(full_html=False, include_plotlyjs='require'))
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def rmsd_heatmap():
        rmsd_mat = rmsd_matrix.get()
        if not analysis_complete.get() or rmsd_mat is None or clust is None:
            return ui.p("Run analysis to see RMSD distance matrix", style="padding: 20px; text-align: center;")
        
        try:
            labels = clust.labels
            level = input.rmsd_heatmap_level()
            show_dendrogram = input.rmsd_dendrogram()
            
            if level == "cluster":
                # Cluster-level heatmap
                cluster_matrix, cluster_ids, cluster_sizes = compute_cluster_distance_matrix(rmsd_mat, labels)
                if show_dendrogram:
                    # Create cluster labels for dendrogram
                    cluster_labels = [f"Cluster {cid} (n={cluster_sizes[cid]})" for cid in cluster_ids]
                    fig = create_distance_heatmap_with_dendrogram(cluster_matrix, labels=None, metric_name='RMSD (Å)', pdb_files=cluster_labels)
                else:
                    fig = create_cluster_distance_heatmap(cluster_matrix, cluster_ids, cluster_sizes, metric_name='RMSD (Å)')
                return ui.HTML(fig.to_html(full_html=False, include_plotlyjs='require'))
            
            elif level == "intra":
                # Single cluster intra-cluster heatmap
                selected_cluster = input.rmsd_intra_cluster()
                if not selected_cluster:
                    return ui.p("Select a cluster to view intra-cluster distances", style="padding: 20px; text-align: center;")
                
                cluster_id = int(selected_cluster)
                fig = create_single_intra_cluster_heatmap(rmsd_mat, labels, cluster_id, clust.pdb_files, metric_name='RMSD (Å)')
                return ui.HTML(fig.to_html(full_html=False, include_plotlyjs='require'))
            
            else:  # pose level
                # Limit to first 500 structures to prevent browser crash
                max_display = 500
                if rmsd_mat.shape[0] > max_display:
                    rmsd_mat_display = rmsd_mat[:max_display, :max_display]
                    labels_display = labels[:max_display]
                    pdb_files_display = clust.pdb_files[:max_display] if clust.pdb_files else None
                    note = f"<div style='padding: 10px; background: #fff3cd; margin-bottom: 10px;'>⚠️ Showing first {max_display} of {rmsd_mat.shape[0]} structures</div>"
                else:
                    rmsd_mat_display = rmsd_mat
                    labels_display = labels
                    pdb_files_display = clust.pdb_files
                    note = ""
                
                if show_dendrogram:
                    fig = create_distance_heatmap_with_dendrogram(rmsd_mat_display, labels_display,
                                                                  metric_name='RMSD (Å)',
                                                                  pdb_files=pdb_files_display)
                else:
                    fig = create_distance_heatmap(rmsd_mat_display, labels_display, metric_name='RMSD (Å)')
                
                return ui.HTML(note + fig.to_html(full_html=False, include_plotlyjs='require'))
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    # ========== INTERACTIVE MAPS TAB ==========
    
    @output
    @render.ui
    def interactive_scatter():
        stats_df = interface_stats.get()
        if stats_df is None or clust is None:
            return ui.p("Run analysis to see interactive plot", style="padding: 20px; text-align: center;")
        
        try:
            method = input.scatter_method()
            # print(f"\n=== Creating {method.upper()} scatter plot ===")
            
            # Use Jaccard matrix for visualization
            fig = create_scatter_multimethod(stats_df, clust.jaccard_distance_matrix, method=method)
            
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id='scatter_plot')
            wrapped_html = f'<div style="width: 100%; height: 600px;">{html}</div>'
            
            return ui.HTML(wrapped_html)
            
        except Exception as e:
            import traceback
            error_msg = f"Error: {e}\n{traceback.format_exc()}"
            # print(error_msg)
            return ui.HTML(f"<pre style='color: red;'>{error_msg}</pre>")
    
    # ========== HOTSPOTS TAB ==========
    
    @output
    @render.ui
    def hotspot_plot_combined():
        hotspots_df = hotspots.get()
        if hotspots_df is None:
            return ui.p("Run analysis to see hotspot analysis", style="padding: 20px; text-align: center;")
        
        try:
            fig = create_hotspot_histogram(hotspots_df, top_n=None, smooth=True, split_by_chain=False)
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id='hotspot_combined')
            wrapped_html = f'<div style="width: 100%; height: 450px;">{html}</div>'
            return ui.HTML(wrapped_html)
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def hotspot_plots_by_chain():
        hotspots_df = hotspots.get()
        if hotspots_df is None:
            return ui.p("Run analysis to see per-chain plots", style="padding: 20px; text-align: center;")
        
        try:
            figs_dict = create_hotspot_histogram(hotspots_df, top_n=None, smooth=True, split_by_chain=True)
            
            if not figs_dict:
                return ui.p("No chain data available")
            
            html_parts = []
            for chain in sorted(figs_dict.keys()):
                fig = figs_dict[chain]
                html = fig.to_html(full_html=False, include_plotlyjs='require', div_id=f'hotspot_chain_{chain}')
                html_parts.append(f'<div style="width: 100%; height: 380px; margin-bottom: 20px;">{html}</div>')
            
            return ui.HTML(''.join(html_parts))
            
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def contact_residue_heatmap():
        """Contact residue heatmap showing intra-cluster interaction patterns"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see contact patterns", style="padding: 20px; text-align: center;")
        
        # Get selected cluster
        selected_cluster = input.contact_heatmap_cluster()
        if not selected_cluster:
            return ui.p("Select a cluster to view contact patterns", style="padding: 20px; text-align: center;")
        
        try:
            cluster_id = int(selected_cluster)
        except (ValueError, TypeError):
            return ui.p("Invalid cluster selection", style="padding: 20px; text-align: center;")
        
        try:
            from visualization.interactive import create_contact_residue_heatmap
            
            # Check if we have contact residue pairs
            if not hasattr(clust, 'contact_residue_pairs') or clust.contact_residue_pairs is None:
                return ui.p("Contact residue data not available", style="padding: 20px; text-align: center; color: #999;")
            
            fig = create_contact_residue_heatmap(
                clust.contact_residue_pairs,
                clust.labels,
                clust.pdb_files,
                cluster_id=cluster_id
            )
            
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id='contact_residue_heatmap')
            wrapped_html = f'<div style="width: 100%; height: 600px;">{html}</div>'
            return ui.HTML(wrapped_html)
            
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error creating contact heatmap: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.ui
    def cluster_hotspot_plot():
        """Per-cluster hotspot plot"""
        if clust is None or not input.cluster_select_hotspot():
            return ui.p("Select a cluster", style="padding: 20px; text-align: center;")
        
        try:
            cluster_id = int(input.cluster_select_hotspot())
            cluster_hotspots = clust.get_binding_hotspots(cluster_id=cluster_id)
            
            if cluster_hotspots.empty:
                return ui.p(f"No data for cluster {cluster_id}")
            
            n_structures = sum(1 for l in clust.labels if l == cluster_id)
            
            # Small cluster warning
            if n_structures <= 3:
                warning_html = f"""
                <div style="padding: 10px; background-color: #fff3cd; border-left: 4px solid #ffc107; 
                           border-radius: 5px; margin-bottom: 10px; color: #212529;">
                    <b>⚠️ Small Cluster</b>: Only {n_structures} structure(s). 
                    Frequencies will be high (multiples of {100/n_structures:.1f}%).
                </div>
                """
            else:
                warning_html = ""
            
            fig = create_hotspot_histogram(cluster_hotspots, top_n=None, smooth=True, split_by_chain=False)
            fig.update_layout(title=f'Cluster {cluster_id} Binding Pattern ({n_structures} structures)')
            
            html = fig.to_html(full_html=False, include_plotlyjs='require', div_id=f'cluster_{cluster_id}_hotspot')
            wrapped_html = f'<div style="width: 100%; height: 450px;">{html}</div>'
            
            return ui.HTML(warning_html + wrapped_html)
            
        except Exception as e:
            import traceback
            return ui.HTML(f"<pre style='color: red;'>Error: {e}\n{traceback.format_exc()}</pre>")
    
    @output
    @render.data_frame
    def consensus_residues_table():
        """Table showing consensus residues for each cluster"""
        summary_df = cluster_summary.get()
        if summary_df is None:
            return pd.DataFrame()
        
        # Extract just cluster and consensus_residues columns
        if 'consensus_residues' in summary_df.columns:
            consensus_df = summary_df[['cluster', 'binding_mode', 'n_structures', 'n_consensus', 'consensus_residues']].copy()
            
            # Rename for clarity
            consensus_df = consensus_df.rename(columns={
                'n_structures': 'Structures',
                'n_consensus': 'Consensus Count',
                'consensus_residues': 'Consensus Residues',
                'binding_mode': 'Binding Mode'
            })
            
            # Apply search filter
            search_term = input.consensus_search()
            if search_term:
                search_lower = search_term.lower()
                mask = (
                    consensus_df['cluster'].astype(str).str.lower().str.contains(search_lower, na=False) |
                    consensus_df['Binding Mode'].astype(str).str.lower().str.contains(search_lower, na=False) |
                    consensus_df['Consensus Residues'].astype(str).str.lower().str.contains(search_lower, na=False) |
                    consensus_df['Structures'].astype(str).str.lower().str.contains(search_lower, na=False) |
                    consensus_df['Consensus Count'].astype(str).str.lower().str.contains(search_lower, na=False)
                )
                consensus_df = consensus_df[mask]
            
            return render.DataGrid(consensus_df, width="100%", height=400)
        else:
            return pd.DataFrame()
    
    @output
    @render.ui
    def hotspot_note():
        note_html = """
        <div style="padding: 15px; background-color: #fff3cd; border-left: 4px solid #ffc107; border-radius: 5px; color: #212529;">
            <b>ℹ️ Note:</b> Different residues may appear at the same position because this shows 
            <b>all structures combined</b>. Use <b>Cmd/Ctrl+F</b> to search for specific residues.
            Check the <b>Overview</b> tab for consensus residues within each binding mode.
        </div>
        """
        return ui.HTML(note_html)
    
    @output
    @render.data_frame
    def hotspot_table():
        hotspots_df = hotspots.get()
        if hotspots_df is None or hotspots_df.empty:
            return pd.DataFrame()
        
        return render.DataGrid(hotspots_df, width="100%", height=400)
    
    # ========== STRUCTURE DETAILS TAB ==========
    
    @output
    @render.ui
    def stats_table_info():
        """Show info about how many rows will be displayed"""
        stats_df = interface_stats.get()
        if stats_df is None or stats_df.empty:
            return ui.p("")
        
        selected_cluster = input.filter_cluster()
        if not selected_cluster:
            return ui.p("")
        
        if selected_cluster == "all":
            n_rows = len(stats_df)
            return ui.HTML(f"""
                <div style='padding: 10px; background: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; 
                           border-radius: 5px; margin: 10px 0;'>
                    <strong>⚠️ Warning:</strong> Displaying ALL {n_rows} structures may freeze your browser!
                    <br><small>Recommend selecting a specific cluster instead.</small>
                </div>
            """)
        else:
            cluster_id = int(selected_cluster)
            n_rows = len(stats_df[stats_df['cluster'] == cluster_id])
            return ui.HTML(f"""
                <div style='padding: 5px 10px; background: #d1ecf1; color: #0c5460; border-radius: 5px; margin: 10px 0;'>
                    Showing {n_rows} structures from Cluster {cluster_id}
                </div>
            """)
    
    @output
    @render.data_frame
    def stats_table():
        stats_df = interface_stats.get()
        if stats_df is None or stats_df.empty:
            return pd.DataFrame()
        
        # Only render if a cluster is selected (prevents initial massive load)
        selected_cluster = input.filter_cluster()
        if not selected_cluster:
            return pd.DataFrame({"Info": ["Select a cluster to view structures"]})
        
        # Filter by cluster if selected
        if selected_cluster != "all":
            cluster_id = int(selected_cluster)
            stats_df = stats_df[stats_df['cluster'] == cluster_id].copy()
        else:
            # Hard limit for "All Clusters" to prevent crashes
            MAX_ROWS = 1000
            if len(stats_df) > MAX_ROWS:
                stats_df = stats_df.head(MAX_ROWS).copy()
                # Add a note in the dataframe
        
        # Create detailed table with contact lists
        display_df = stats_df.copy()
        
        # Keep full residue lists for display (don't truncate)
        # The DataGrid will handle scrolling and the user can see full text
        
        # Get motif residues for highlighting (if available)
        motif_dict = None
        if hasattr(input, 'motif_residues') and input.motif_residues():
            motif_input = input.motif_residues().strip()
            if motif_input:
                try:
                    motif_dict = parse_motif_residues(motif_input)
                except:
                    pass
        
        # Start with basic columns (contacts will be added at the end)
        cols = ['structure', 'cluster', 'n_protein_residues', 'n_nucleic_residues', 'representative']
        
        # Add motif match % BEFORE ΔSASA columns (as requested)
        if 'motif_match_pct' in display_df.columns:
            display_df['Motif Match (%)'] = display_df['motif_match_pct'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "")
            cols.append('Motif Match (%)')
        
        # Add ΔSASA columns if available - format to 4 s.f.
        if 'total_delta_sasa' in display_df.columns:
            display_df['Buried ΔSASA (Å²)'] = display_df['total_delta_sasa'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "")
            cols.append('Buried ΔSASA (Å²)')
        if 'motif_delta_sasa' in display_df.columns:
            display_df['Motif ΔSASA (Å²)'] = display_df['motif_delta_sasa'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "")
            cols.append('Motif ΔSASA (Å²)')
        
        # Add per-chain ΔSASA columns if available (split mode) - format to 4 s.f.
        chain_sasa_cols = [col for col in display_df.columns if col.startswith('chain_') and col.endswith('_sasa')]
        for col in sorted(chain_sasa_cols):  # Sort to ensure consistent order
            # Extract chain ID from column name (e.g., 'chain_A_sasa' -> 'A')
            chain_id = col.replace('chain_', '').replace('_sasa', '')
            new_name = f'Chain {chain_id} ΔSASA (Å²)'
            display_df[new_name] = display_df[col].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "")
            cols.append(new_name)
        
        # Add FULL contact lists at the end (no truncation - DataGrid handles display)
        # Bold motif residues in protein contacts if motif is specified
        if 'protein_residues' in display_df.columns:
            if motif_dict:
                # Bold motif residues
                def bold_motif_residues(res_str):
                    if pd.isna(res_str) or not res_str:
                        return res_str
                    residues = res_str.split(', ')
                    bolded_residues = []
                    for res in residues:
                        # Parse residue format: "CHAIN:RESNAME123"
                        is_motif = False
                        if ':' in res:
                            chain_id, res_part = res.split(':', 1)
                            # Extract residue number from res_part (e.g., "ARG123" -> 123)
                            res_num_str = ''.join(filter(str.isdigit, res_part))
                            if res_num_str and chain_id in motif_dict:
                                res_num = int(res_num_str)
                                if res_num in motif_dict[chain_id]:
                                    is_motif = True
                        
                        if is_motif:
                            bolded_residues.append(f'**{res}**')
                        else:
                            bolded_residues.append(res)
                    return ', '.join(bolded_residues)
                
                display_df['Protein Contacts'] = display_df['protein_residues'].apply(bold_motif_residues)
            else:
                display_df['Protein Contacts'] = display_df['protein_residues']
            cols.append('Protein Contacts')
        if 'nucleic_residues' in display_df.columns:
            display_df['Nucleic Contacts'] = display_df['nucleic_residues']
            cols.append('Nucleic Contacts')
        
        available_cols = [c for c in cols if c in display_df.columns]
        display_df = display_df[available_cols]
        
        # Rename for clarity
        rename_map = {
            'structure': 'Structure',
            'cluster': 'Cluster',
            'n_protein_residues': 'N Protein',
            'n_nucleic_residues': 'N Nucleic',
            'representative': 'Representative'
        }
        display_df = display_df.rename(columns=rename_map)
        
        return render.DataGrid(display_df, width="100%", height=600, filters=True)
    
    # ========== 3D VIEWER TAB ==========
    
    @output
    @render.ui
    def structure_info():
        if clust is None or not input.structure_select():
            return ui.p("")
        
        try:
            structure_name = input.structure_select()
            stats_df = interface_stats.get()
            if stats_df is not None:
                struct_info = stats_df[stats_df['structure'] == structure_name].iloc[0]
                
                info_html = f"""
                <div style="padding: 10px; background-color: #f8f9fa; margin-bottom: 10px; border-radius: 5px; border-left: 4px solid #0d6efd;">
                    <b>{structure_name}</b> - <span style="color: #0d6efd;">Cluster {struct_info['cluster']}</span><br>
                    <small>Protein residues: {struct_info['n_protein_residues']} | DNA/RNA residues: {struct_info['n_nucleic_residues']}</small><br>
                    <small>Interface: {struct_info.get('protein_residues', '')[:80]}...</small>
                </div>
                """
                return ui.HTML(info_html)
        except Exception as e:
            pass
        return ui.p("")
    
    @output
    @render.ui
    def structure_viewer():
        # Cache clust reference at start to prevent it becoming None during render
        clust_ref = clust
        
        if clust_ref is None:
            return ui.p("Run analysis first", style="padding: 20px; text-align: center;")
        
        if not input.structure_select():
            return ui.p("Select a structure to view in 3D", style="padding: 20px; text-align: center;")
        
        try:
            structure_name = input.structure_select()
            pdb_file = None
            for f in clust_ref.pdb_files:
                if os.path.basename(f) == structure_name:
                    pdb_file = f
                    break
            
            if not pdb_file or not os.path.exists(pdb_file):
                return ui.p(f"Structure not found: {structure_name}")
            
            # Always create a SASA-colored PDB for Molstar (it can use B-factor data)
            pdb_to_display = pdb_file
            sasa_min = 0
            sasa_max = 100
            
            # Create a PDB with ΔSASA values in B-factor column
            try:
                log_status(f"Creating ΔSASA visualization for {structure_name}...")
                
                # Get contact data for this structure to identify chains correctly
                structure_idx = None
                for i, f in enumerate(clust_ref.pdb_files):
                    if os.path.basename(f) == structure_name:
                        structure_idx = i
                        break
                
                contacts = None
                if structure_idx is not None and clust_ref.contact_residue_pairs:
                    contacts = clust_ref.contact_residue_pairs[structure_idx]
                
                motif = current_motif.get()
                
                sasa_pdb = InterfaceAnalyzer.create_sasa_viewer_pdb(
                    pdb_file, 
                    motif_residues=motif,
                    contacts=contacts
                )
                
                if sasa_pdb and os.path.exists(sasa_pdb):
                    pdb_to_display = sasa_pdb
                    log_status(f"✓ ΔSASA PDB created: {os.path.basename(sasa_pdb)}")
                    
                    # Read metadata for color range
                    meta_file = sasa_pdb.replace('.pdb', '_meta.txt')
                    if os.path.exists(meta_file):
                        with open(meta_file, 'r') as f:
                            for line in f:
                                if line.startswith('min='):
                                    sasa_min = float(line.split('=')[1])
                                elif line.startswith('max='):
                                    sasa_max = float(line.split('=')[1])
                        log_status(f"✓ Color range: {sasa_min:.2f} - {sasa_max:.2f} Å²")
                else:
                    log_status("⚠️ Could not create ΔSASA visualization - using original PDB")
            
            except Exception as e:
                import traceback
                error_msg = f"Error creating ΔSASA visualization: {str(e)}"
                log_status(f"❌ {error_msg}")
                # Fall back to original PDB
                pdb_to_display = pdb_file
            
            # Create 3D viewer with Molstar (display settings now handled by Molstar UI)
            # SASA coloring is enabled by default if SASA data is present
            viewer_html = create_mol_viewer_html(
                pdb_to_display, 
                width=900, 
                height=650,
                protein_style='cartoon',  # Default - Molstar allows runtime changes
                protein_color='chain',     # Default to chain coloring
                dna_style='stick',         # Default - Molstar allows runtime changes
                show_surface=False,        # Default - can be toggled in Molstar
                surface_opacity=0.7,       # Default
                zoom_level=1.0,            # Default
                color_by_sasa=True,        # Always show SASA if available (B-factor data)
                sasa_min=sasa_min,
                sasa_max=sasa_max
            )
            
            return ui.HTML(viewer_html)
        
        except Exception as e:
            import traceback
            log_status(f"❌ Error loading structure: {e}")
            return ui.p(f"Error loading structure: {e}\n\n{traceback.format_exc()}")


# Helper function: Compute RMSD QC stats
def compute_rmsd_qc(labels, rmsd_matrix):
    """
    Compute inter-cluster RMSD statistics for QC
    
    Uses a sequential arrangement approach (greedy nearest-neighbor path)
    to better represent cluster tightness, rather than averaging ALL pairwise distances.
    
    Think of it like: A → B → C → D (chain) instead of averaging all pairs (AB, AC, AD, BC, BD, CD)
    This gives a better sense of how structures are "arranged" in sequence.
    
    Also computes the medoid: the structure with minimum average RMSD to all others in the cluster.
    The medoid is the most representative structure.
    """
    unique_clusters = sorted(set(labels))
    rows = []
    
    for cluster_id in unique_clusters:
        indices = [i for i, l in enumerate(labels) if l == cluster_id]
        
        if len(indices) < 2:
            rows.append({
                'cluster': cluster_id,
                'mean_rmsd': 0.0,
                'std_rmsd': 0.0,
                'mean_rmsd_sequential': 0.0,
                'median_rmsd': 0.0,
                'max_rmsd': 0.0,
                'min_rmsd': 0.0,
                'medoid_idx': indices[0] if indices else -1,
                'medoid_avg_rmsd': 0.0
            })
            continue
        
        # Method 1: ALL pairwise RMSD (old method - for comparison)
        cluster_rmsds_all = []
        for i in range(len(indices)):
            for j in range(i+1, len(indices)):
                cluster_rmsds_all.append(rmsd_matrix[indices[i], indices[j]])
        
        # Method 2: SEQUENTIAL arrangement (greedy nearest-neighbor path)
        # Start with first structure, build a path by always picking nearest unvisited neighbor
        visited = [indices[0]]
        unvisited = set(indices[1:])
        sequential_rmsds = []
        
        current = indices[0]
        while unvisited:
            # Find nearest unvisited neighbor
            min_dist = float('inf')
            nearest = None
            for candidate in unvisited:
                dist = rmsd_matrix[current, candidate]
                if dist < min_dist:
                    min_dist = dist
                    nearest = candidate
            
            if nearest is not None:
                sequential_rmsds.append(min_dist)
                visited.append(nearest)
                unvisited.remove(nearest)
                current = nearest
        
        # Method 3: MEDOID computation
        # Find the structure with minimum average RMSD to all others in cluster
        medoid_idx = None
        min_avg_rmsd = float('inf')
        
        for idx in indices:
            avg_rmsd = np.mean([rmsd_matrix[idx, other] for other in indices if other != idx])
            if avg_rmsd < min_avg_rmsd:
                min_avg_rmsd = avg_rmsd
                medoid_idx = idx
        
        rows.append({
            'cluster': cluster_id,
            'mean_rmsd': np.mean(cluster_rmsds_all),  # Mean of all pairwise RMSDs within cluster
            'std_rmsd': np.std(cluster_rmsds_all),    # Std dev of pairwise RMSDs (cluster variability)
            'mean_rmsd_sequential': np.mean(sequential_rmsds) if sequential_rmsds else 0.0,  # Sequential path
            'median_rmsd': np.median(cluster_rmsds_all),
            'max_rmsd': np.max(cluster_rmsds_all),
            'min_rmsd': np.min(cluster_rmsds_all),
            'medoid_idx': medoid_idx,
            'medoid_avg_rmsd': min_avg_rmsd
        })
    
    return pd.DataFrame(rows)


# Create app
app = App(app_ui, server)

if __name__ == "__main__":
    print("=" * 70)
    print("🧬 ClustalDM - Main Interactive Application")
    print("=" * 70)
    print("\n✨ Strategy:")
    print("  1. Compute BOTH Jaccard (chain-aware) + RMSD")
    print("  2. Cluster by Jaccard ONLY (binding mode)")
    print("  3. Use RMSD as QC metric (inter-cluster similarity)")
    print("\n🚀 Starting server...")
    print("=" * 70)
    print()
