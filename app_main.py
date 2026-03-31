#!/usr/bin/env python3
"""
ClustalDM - Main Interactive Application
Strategy: Compute both RMSD + Jaccard, cluster by Jaccard only, use RMSD as QC
"""
import webbrowser

import os
import sys
import re
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
    from visualization import (
        create_scatter_multimethod,
        create_hotspot_histogram,
        create_molstar_viewer_html as create_mol_viewer_html,  # Use Molstar instead of py3Dmol
        create_distance_heatmap,
        create_cluster_size_distribution,
    )
    from visualization.sequence_alignment import (
        create_alignment_visualization_medoids,
        create_alignment_visualization_within_cluster
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
        ui.h1("🧬 ClustalDMɑ - AlphaFold2 Model Clustering", style="display: inline-block; margin: 20px;"),
        ui.download_button("download_all", "📥 Download Results", class_="btn btn-sm btn-success", style="margin-right: 10px;"),
        ui.input_action_button("theme_toggle", "🌙 Dark", class_="btn btn-sm btn-outline-secondary"),
        style="position: relative;"
    ),
    
    ui.layout_sidebar(
        ui.sidebar(
            ui.h4("Upload AlphaFold2 Models"),
            ui.div(
                ui.input_text("af2_directory", "📁 AlphaFold2 Model Directory:",
                             value="",
                             placeholder="/path/to/alphafold/models"),
                ui.tags.small("Paste directory path containing PDB/CIF files (searches subdirectories automatically)", 
                             style="color: #6c757d; display: block; margin-top: 5px; margin-bottom: 15px;"),
                style="margin-bottom: 10px;"
            ),
            ui.div(
                ui.input_file("reference_model", "Reference Model (optional):", 
                             accept=[".pdb", ".cif"]),
                ui.tags.small("Upload reference structure for superposition/RMSD calculation", 
                             style="color: #6c757d; display: block; margin-top: 5px; margin-bottom: 10px;"),
                ui.input_checkbox(
                    "reference_is_alphafold",
                    "Reference is AlphaFold (use B-factor as pLDDT in reference lane)",
                    value=False
                ),
                ui.tags.small(
                    "If unchecked, reference lane hides pLDDT/B-factor to avoid misinterpreting crystal B-factors.",
                    style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 10px;"
                ),
                style="margin-bottom: 10px;"
            ),
            ui.div(
                ui.input_file("reference_sequence", "Reference Sequence (optional):", 
                             accept=[".fasta", ".fa", ".txt", ".pdb", ".cif"]),
                ui.tags.small("Upload reference FASTA, PDB, or CIF for residue number validation/correction (recommended for chopped AlphaFold models)", 
                             style="color: #6c757d; display: block; margin-top: 5px; margin-bottom: 10px;"),
                style="margin-bottom: 10px;"
            ),
            ui.output_ui("upload_status"),
            ui.hr(),
            ui.h5("Analysis Options"),
            ui.input_text("protein_chains", "Protein Chain ID(s):", 
                         value="A",
                         placeholder="e.g., A or A,B"),
            ui.tags.small("Provide receptor protein chain labels (comma-separated). All other chains are treated as partner/ligand.", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 10px;"),
            ui.input_text("motif_residues", "Motif Residues (optional):", 
                         value="",
                         placeholder="e.g., A:10-20,B:30-35"),
            ui.tags.small("Format: Chain:Start-End, comma-separated", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 5px;"),
            ui.input_checkbox("apply_residue_offset", "Apply Residue Number Corrections", value=False),
            ui.tags.small("Enable for 'chopped' AlphaFold models (extracts offset from filename/alignment)", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 10px;"),
            ui.hr(),
            ui.h5("Clustering (HDBSCAN)"),
            ui.input_slider("min_cluster_size", "Min Cluster Size:", 
                           min=2, max=20, value=5, step=1),
            ui.input_slider("min_samples_hdb", "Min Samples (noise threshold):", 
                           min=1, max=10, value=2, step=1),
            ui.input_checkbox("filter_duplicates", "Filter Near-Duplicates", value=True),
            ui.tags.small("Remove very similar structures before clustering (faster, cleaner trees)", 
                         style="color: #6c757d; display: block; margin-top: -5px; margin-bottom: 10px;"),
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
                    - Color: By cluster assignment
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
                "🧬 Sequence QC",
                ui.card(
                    ui.card_header("Sequence Alignment Quality Control"),
                    ui.markdown("""
                        **Visualize secondary structure and pLDDT confidence for:**
                        - **Medoid Comparison**: Compare secondary structure across cluster medoids
                        - **Within Cluster**: Examine all structures within a single cluster
                        
                        *Colour Key*  
                        *pLDDT score*: Dark blue (>90), Light blue (70-90), Yellow (50-70), Orange (<50)  
                        *Secondary Structure*:  
                            α-helix (H, red), 3-10 helix (G, orange), π-helix (I, pink),  
                            β-sheet (E, blue), β-bridge (B, light blue),  
                            Turn (T, green), Bend (S, yellow), Coil/Undefined (-, gray).  
                    """),
                    ui.input_radio_buttons(
                        "sequence_viz_mode",
                        "Visualization Mode:",
                        choices={
                            "medoids": "Medoid Comparison (across clusters)",
                            "within": "Within Cluster (all structures)"
                        },
                        selected="medoids",
                        inline=True
                    ),
                    ui.panel_conditional(
                        "input.sequence_viz_mode === 'within'",
                        ui.input_select(
                            "sequence_viz_cluster",
                            "Select Cluster:",
                            choices={},
                            width="300px"
                        )
                    ),
                    ui.output_ui("sequence_alignment_plot")
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
    tm_score_stats = reactive.Value(None)
    current_motif = reactive.Value(None)  # Store motif for visualization
    analysis_complete = reactive.Value(False)  # Track if analysis is done
    status_log = reactive.Value([])  # Terminal-style output log
    current_theme = reactive.Value('dark')  # Track current theme
    
    # AlphaFold2-specific reactive values
    af2_models_dir = reactive.Value(None)  # Directory containing AF2 models
    reference_model_file = reactive.Value(None)  # Optional reference model for RMSD
    reference_sequence_file = reactive.Value(None)  # Optional reference FASTA for residue correction
    af2_model_files = reactive.Value([])  # List of found AF2 model files
    
    def log_status(message):
        """Add message to status log"""
        current_log = status_log.get()
        current_log.append(message)
        status_log.set(current_log)

    def parse_protein_chain_labels(chain_text):
        """Parse comma/space separated chain labels from user input."""
        if not chain_text:
            return ['A']
        labels = [c.strip() for c in re.split(r'[\s,;]+', chain_text.strip()) if c.strip()]
        return [c.upper() for c in labels] if labels else ['A']

    def build_chain_selections(pdb_file, protein_chain_labels):
        """
        Build MDTraj selections using user-provided protein chain labels.
        Partner/ligand is defined as all non-water atoms NOT in protein chains.
        """
        import mdtraj as md

        traj = md.load(pdb_file)
        topology = traj.topology

        chain_label_to_indices = {}
        all_chain_labels = []
        for chain in topology.chains:
            chain_label = str(chain.chain_id).strip()
            chain_key = chain_label.upper()
            chain_label_to_indices.setdefault(chain_key, []).append(chain.index)
            all_chain_labels.append(chain_label)

        # Resolve requested chain labels to MDTraj chain indices
        protein_chain_indices = []
        missing_labels = []
        for label in protein_chain_labels:
            key = label.upper()
            if key in chain_label_to_indices:
                protein_chain_indices.extend(chain_label_to_indices[key])
            else:
                missing_labels.append(label)

        protein_chain_indices = sorted(set(protein_chain_indices))

        # Fallback: if none matched, auto-detect protein chains from topology
        if len(protein_chain_indices) == 0:
            auto_protein_indices = []
            auto_protein_labels = []
            for chain in topology.chains:
                residues = list(chain.residues)
                if any(getattr(res, 'is_protein', False) for res in residues):
                    auto_protein_indices.append(chain.index)
                    auto_protein_labels.append(str(chain.chain_id).strip())

            if auto_protein_indices:
                protein_chain_indices = sorted(set(auto_protein_indices))
                protein_chain_labels = sorted(set([lbl.upper() for lbl in auto_protein_labels]))

        if len(protein_chain_indices) == 0:
            raise ValueError(
                f"Could not resolve protein chains from input {protein_chain_labels}. "
                f"Available chains: {sorted(set(all_chain_labels))}"
            )

        chain_clause = ' or '.join([f"chainid == {idx}" for idx in protein_chain_indices])
        protein_selection = f"protein and ({chain_clause})"
        partner_selection = f"not water and not ({chain_clause})"

        return {
            'protein_selection': protein_selection,
            'partner_selection': partner_selection,
            'protein_chain_indices': protein_chain_indices,
            'missing_labels': missing_labels,
            'available_labels': sorted(set(all_chain_labels)),
            'resolved_labels': protein_chain_labels
        }
    
    @output
    @render.ui
    def upload_status():
        """Show upload status"""
        models_dir = af2_models_dir.get()
        ref_model = reference_model_file.get()
        model_files = af2_model_files.get()
        
        status_items = []
        
        # Check directory input
        if models_dir and os.path.isdir(models_dir):
            status_items.append(ui.tags.li(f"✓ Models directory: {models_dir}", style="color: #28a745;"))
            if model_files:
                status_items.append(ui.tags.li(f"✓ Found {len(model_files)} AlphaFold2 model(s)", style="color: #28a745;"))
            else:
                status_items.append(ui.tags.li("⚠ No PDB/mmCIF files found in directory", style="color: #ffc107;"))
        else:
            status_items.append(ui.tags.li("⚠ No models directory specified", style="color: #6c757d;"))
        
        if ref_model:
            status_items.append(ui.tags.li("✓ Reference model uploaded", style="color: #28a745;"))
        else:
            status_items.append(ui.tags.li("ℹ Optional: Reference model for RMSD", style="color: #6c757d;"))
        
        # Show ready status
        if models_dir and model_files:
            status_items.append(ui.tags.li(f"👉 Ready! Click 'Start Analysis' to cluster {len(model_files)} models", 
                                          style="color: #007bff; font-weight: bold;"))
        
        return ui.tags.ul(*status_items, style="list-style: none; padding-left: 0;")
    
    @reactive.effect
    @reactive.event(input.af2_directory)
    def handle_directory_change():
        """Handle AF2 directory input change"""
        directory = input.af2_directory().strip()
        if not directory:
            af2_models_dir.set(None)
            af2_model_files.set([])
            return
        
        # Expand ~ to home directory
        directory = os.path.expanduser(directory)
        
        if not os.path.isdir(directory):
            log_status(f"⚠ Not a valid directory: {directory}")
            af2_models_dir.set(None)
            af2_model_files.set([])
            return
        
        try:
            # Find all PDB and mmCIF files recursively
            log_status(f"Scanning {directory} and subdirectories...")
            model_files = find_pdb_files(directory)
            
            if model_files:
                af2_models_dir.set(directory)
                af2_model_files.set(model_files)
                log_status(f"✓ Found {len(model_files)} AlphaFold2 models in {directory} (including subdirectories)")
                ui.notification_show(
                    f"✓ Found {len(model_files)} model file(s)",
                    type="success",
                    duration=3
                )
            else:
                af2_models_dir.set(directory)
                af2_model_files.set([])
                log_status(f"⚠ No PDB/mmCIF files found in {directory}")
                ui.notification_show(
                    "No PDB/mmCIF files found in directory",
                    type="warning",
                    duration=3
                )
            
        except Exception as e:
            log_status(f"Directory scan error: {str(e)}")
            ui.notification_show(f"Error scanning directory: {str(e)}", type="error", duration=5)
    
    @reactive.effect
    @reactive.event(input.reference_model)
    def handle_reference_upload():
        """Handle reference model upload"""
        file_info = input.reference_model()
        if not file_info:
            return
        
        try:
            # Create temporary directory if needed
            if not os.path.exists(tempfile.gettempdir()):
                os.makedirs(tempfile.gettempdir())
            
            # Save reference model
            ref_path = os.path.join(tempfile.gettempdir(), f"reference_{file_info[0]['name']}")
            shutil.copy(file_info[0]['datapath'], ref_path)
            reference_model_file.set(ref_path)
            log_status(f"Uploaded reference: {file_info[0]['name']}")
            
            ui.notification_show(
                f"✓ Reference model uploaded: {file_info[0]['name']}",
                type="success",
                duration=3
            )
            
        except Exception as e:
            log_status(f"Reference upload error: {str(e)}")
            ui.notification_show(f"Error uploading reference: {str(e)}", type="error", duration=5)
    
    @reactive.effect
    @reactive.event(input.reference_sequence)
    def handle_reference_sequence_upload():
        """Handle reference sequence (FASTA or PDB) upload"""
        file_info = input.reference_sequence()
        if not file_info:
            return
        
        try:
            # Create temporary directory if needed
            if not os.path.exists(tempfile.gettempdir()):
                os.makedirs(tempfile.gettempdir())
            
            # Save reference file
            filename = file_info[0]['name']
            ref_path = os.path.join(tempfile.gettempdir(), f"reference_sequence_{filename}")
            shutil.copy(file_info[0]['datapath'], ref_path)
            reference_sequence_file.set(ref_path)
            log_status(f"Uploaded reference: {filename}")
            
            ui.notification_show(
                f"✓ Reference sequence uploaded: {file_info[0]['name']}",
                type="success",
                duration=3
            )
            
        except Exception as e:
            log_status(f"Reference sequence upload error: {str(e)}")
            ui.notification_show(f"Error uploading reference sequence: {str(e)}", type="error", duration=5)
            ui.notification_show(f"Reference upload failed: {str(e)}", type="error", duration=5)
    
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

            if tm_score_stats.get() is not None:
                csv_buffer = tm_score_stats.get().to_csv(index=False)
                zipf.writestr("tm_score_vs_reference.csv", csv_buffer)
            
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
    @reactive.event(input.run_clustering)
    def run_analysis():
        global clust, pdb_dir_path
        
        # Import necessary modules at function start
        import os
        import MDAnalysis as mda
        
        # Clear previous log
        status_log.set([])
        
        # Get AF2 models directory and files
        models_dir = af2_models_dir.get()
        model_files = af2_model_files.get()
        ref_model = reference_model_file.get()
        
        # Validate inputs
        if not models_dir or not model_files:
            ui.notification_show("⚠️ Please specify a valid directory with AlphaFold2 models!", type="error", duration=5)
            return
        
        # Show immediate notification for job submission
        ui.notification_show("🚀 Running clustering analysis...", type="message", duration=3)
        
        # Get HDBSCAN parameters
        min_cluster_size = input.min_cluster_size()
        min_samples_hdb = input.min_samples_hdb()
        
        # Show input acknowledgement
        log_status("="*50)
        log_status("🚀 Starting Clustering Analysis (AlphaFold2 Models)")
        log_status("="*50)
        log_status(f"Total structures: {len(model_files)}")
        log_status(f"HDBSCAN min_cluster_size: {min_cluster_size}")
        log_status(f"HDBSCAN min_samples: {min_samples_hdb}")
        if ref_model:
            log_status(f"Reference model: {os.path.basename(ref_model)}")
        log_status("")
        
        pdb_dir_path = models_dir
        pdb_files = model_files
        
        try:
            log_status(f"✓ Analyzing {len(pdb_files)} structures")
            log_status("")
            
            # Initialize clusterer
            log_status("Initializing clusterer...")
            ui.notification_show("📊 Initializing clusterer...", type="message", duration=2)
            from core.analysis import InterfaceAnalyzer
            from core.io_utils import read_fasta_sequence, extract_sequence_from_pdb
            
            clust = PDBClusterer(pdb_files)
            
            # Load reference sequence if provided (FASTA or PDB)
            ref_seq_file = reference_sequence_file.get()
            ref_model_file = reference_model_file.get()
            
            print(f"\n=== Reference Loading ===")
            print(f"  ref_seq_file: {ref_seq_file}")
            print(f"  ref_model_file: {ref_model_file}")
            print(f"  ref_seq exists: {os.path.exists(ref_seq_file) if ref_seq_file else False}")
            print(f"  ref_model exists: {os.path.exists(ref_model_file) if ref_model_file else False}")
            
            # Initialize reference attributes
            clust.reference_sequence = None
            clust.reference_pdb = None  # Store PDB path if reference is a PDB
            
            # Priority 1: Reference sequence file (can be FASTA, PDB, or CIF)
            if ref_seq_file and os.path.exists(ref_seq_file):
                try:
                    # Check if it's a PDB/CIF or FASTA file
                    if ref_seq_file.lower().endswith(('.pdb', '.cif')):
                        # Extract sequence from user-provided protein chains
                        prot_chains = parse_protein_chain_labels(input.protein_chains())
                        ref_sequence = extract_sequence_from_pdb(ref_seq_file, chain_id=prot_chains)
                        log_status(f"✓ Reference sequence extracted from structure: {len(ref_sequence)} residues")
                        print(f"  Extracted from structure (protein chains {prot_chains}): {len(ref_sequence)} residues")
                        print(f"  First 50 chars: {ref_sequence[:50]}")
                        # Store structure path for visualization
                        clust.reference_pdb = ref_seq_file
                        print(f"  Stored reference_pdb path: {ref_seq_file}")
                    else:
                        ref_sequence = read_fasta_sequence(ref_seq_file)
                        log_status(f"✓ Reference sequence loaded from FASTA: {len(ref_sequence)} residues")
                        print(f"  Loaded from FASTA: {len(ref_sequence)} residues")
                        print(f"  First 50 chars: {ref_sequence[:50]}")
                    
                    # Explicitly set the reference sequence
                    clust.reference_sequence = ref_sequence
                    print(f"  Set clust.reference_sequence to string of length {len(ref_sequence)}")
                    print(f"  Verification: clust.reference_sequence is None = {clust.reference_sequence is None}")
                    log_status("  → Will validate and correct residue numbering")
                except Exception as e:
                    import traceback
                    log_status(f"Warning: Could not load reference sequence: {e}")
                    log_status("  → Proceeding without sequence validation")
                    print(f"  ERROR: {e}")
                    traceback.print_exc()
                    clust.reference_sequence = None
                    clust.reference_pdb = None
            
            # Priority 2: Reference model file (if no sequence file provided but model exists)
            if clust.reference_pdb is None and ref_model_file and os.path.exists(ref_model_file):
                print(f"  Using reference model for visualization: {ref_model_file}")
                clust.reference_pdb = ref_model_file
                # Also extract sequence from model if not already done
                if clust.reference_sequence is None:
                    try:
                        prot_chains = parse_protein_chain_labels(input.protein_chains())
                        ref_sequence = extract_sequence_from_pdb(ref_model_file, chain_id=prot_chains)
                        clust.reference_sequence = ref_sequence
                        print(f"  Extracted sequence from reference model (protein chains {prot_chains}): {len(ref_sequence)} residues")
                    except Exception as e:
                        print(f"  Could not extract sequence from reference model: {e}")
            
            if clust.reference_sequence is None and clust.reference_pdb is None:
                print(f"  No reference file provided or file doesn't exist")
            
            log_status("✓ Clusterer initialized")
            log_status("")
            
            # Parse motif residues early for Jaccard screening
            motif_input = input.motif_residues().strip()
            motif_dict = InterfaceAnalyzer.parse_motif_residues(motif_input) if motif_input else None
            
            # Get receptor protein chains from input
            protein_chains = parse_protein_chain_labels(input.protein_chains())
            log_status(f"Protein chain(s): {', '.join(protein_chains)}")

            # Build chain-aware selections from first structure
            chain_selection_info = build_chain_selections(pdb_files[0], protein_chains)
            protein_selection = chain_selection_info['protein_selection']
            ligand_selection = chain_selection_info['partner_selection']

            if chain_selection_info['missing_labels']:
                log_status(
                    f"  ⚠️ Chain label(s) not found: {', '.join(chain_selection_info['missing_labels'])} "
                    f"(available: {', '.join(chain_selection_info['available_labels'])})"
                )
            log_status(f"  Resolved protein chain indices: {chain_selection_info['protein_chain_indices']}")

            # Store protein chain config in clusterer for downstream views
            clust.protein_chains = protein_chains
            clust.protein_chain_indices = chain_selection_info['protein_chain_indices']
            
            # STEP 1: Compute Jaccard matrix (chain-aware!)
            log_status("STEP 1: Computing Jaccard contact matrix...")
            ui.notification_show("🧬 STEP 1: Computing contact matrix...", type="message", duration=3)
            if motif_dict:
                log_status(f"  🎯 MOTIF MODE: Local search enabled")
                log_status(f"  Motif: {motif_input}")
                log_status(f"  → Off-target poses (no motif contact) will be filtered")
            else:
                log_status(f"  🌍 GLOBAL MODE: All interface contacts")
            
            log_status(f"  Target (protein) selection: {protein_selection}")
            log_status(f"  Partner selection (NOT protein chains): {ligand_selection}")
            
            # Get residue offset toggle
            apply_offset = input.apply_residue_offset()
            if apply_offset:
                log_status(f"  📊 Residue offset correction: ENABLED (for chopped models)")
            else:
                log_status(f"  📊 Residue offset correction: DISABLED (full-length models)")
            
            clust.compute_jaccard_contact_matrix(
                n_jobs=-1, 
                motif_residues=motif_dict,
                protein_selection=protein_selection,
                nucleic_selection=ligand_selection,
                apply_offset=apply_offset
            )  # Use all CPU cores
            log_status("✓ Jaccard matrix computed")
            ui.notification_show("✓ Contact matrix done", type="message", duration=2)
            log_status("")

            # STEP 1.5: TM-score QC against provided reference structure (optional)
            tm_score_stats.set(None)
            if clust.reference_pdb and os.path.exists(clust.reference_pdb):
                log_status("STEP 1.5: Computing TM-score vs reference...")
                try:
                    tm_df = compute_tm_scores_to_reference(
                        reference_pdb=clust.reference_pdb,
                        pdb_files=clust.pdb_files,
                        protein_chains=protein_chains
                    )
                    if tm_df is not None and len(tm_df) > 0:
                        tm_score_stats.set(tm_df)
                        valid_tm = tm_df['tm_score'].dropna()
                        if len(valid_tm) > 0:
                            log_status(
                                f"✓ TM-score computed (mean={valid_tm.mean():.3f}, "
                                f"min={valid_tm.min():.3f}, max={valid_tm.max():.3f})"
                            )
                        else:
                            log_status("⚠️ TM-score produced no valid alignments")
                    else:
                        log_status("⚠️ TM-score skipped (no valid data)")
                except Exception as e:
                    log_status(f"⚠️ TM-score computation failed: {str(e)}")
                log_status("")
            
            # STEP 2: Compute RMSD matrix
            log_status("STEP 2: Computing RMSD matrix...")
            ui.notification_show("📏 STEP 2: Computing RMSD...", type="message", duration=3)
            # Load structures for RMSD calculation
            clust.load_structures()
            
            # Determine reference index for RMSD calculation
            ref_idx = 0  # Default to first structure
            if ref_model and os.path.exists(ref_model):
                try:
                    # Add reference model to the beginning of universes for superposition
                    from core.clusterer import load_universe_with_cif_support
                    ref_universe = load_universe_with_cif_support(ref_model)
                    log_status(f"  Using reference model for superposition: {os.path.basename(ref_model)}")
                    # We'll use it for alignment but not include in clustering
                except Exception as e:
                    log_status(f"  ⚠️ Could not load reference model, using first structure: {e}")
                    ref_universe = None
            
            # Compute RMSD with reference-based superposition if available
            # For AlphaFold, use 'protein and name CA' selection (standard)
            original_selection = clust.selection
            clust.selection = 'protein and name CA'
            
            clust.compute_distance_matrix(align_structures=True, reference_idx=ref_idx, n_jobs=-1)
            rmsd_mat = clust.distance_matrix.copy()  # Save RMSD matrix
            rmsd_matrix.set(rmsd_mat)
            clust.selection = original_selection  # Restore original
            log_status("✓ RMSD matrix computed")
            ui.notification_show("✓ RMSD done", type="message", duration=2)
            log_status("")
            
            # Store motif for later use in visualization
            current_motif.set(motif_dict)
            
            # STEP 3: Cluster by Jaccard only
            log_status("STEP 3: Clustering by Jaccard distance...")
            ui.notification_show("🎯 STEP 3: Running HDBSCAN clustering...", type="message", duration=3)
            # Restore Jaccard matrix for clustering
            clust.distance_matrix = clust.jaccard_distance_matrix
            clust.metric_name = "Jaccard"
            
            filter_dupes = input.filter_duplicates()
            clust.cluster_hdbscan(
                min_cluster_size=min_cluster_size, 
                min_samples=min_samples_hdb,
                filter_duplicates=filter_dupes
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
            
            # Add motif match scores if available
            if hasattr(clust, 'motif_match_scores') and clust.motif_match_scores is not None:
                match_scores = clust.motif_match_scores
                if len(match_scores) > 0 and len(match_scores) == len(stats_df):
                    stats_df['motif_match_pct'] = match_scores
                    log_status("✓ Motif match scores added to interface stats")

            # Add TM-score per structure if available
            tm_df = tm_score_stats.get()
            if tm_df is not None and len(tm_df) > 0:
                tm_map = dict(zip(tm_df['structure'], tm_df['tm_score']))
                stats_df['tm_score'] = stats_df['structure'].map(tm_map)
                log_status("✓ TM-score added to interface stats")
            
            interface_stats.set(stats_df)
            
            # Get hotspots
            hotspots_df = clust.get_binding_hotspots()
            hotspots.set(hotspots_df)
            log_status("✓ Hotspot analysis complete")
            log_status("")
            
            # Get cluster summary with RMSD QC
            summary_df = clust.get_cluster_summary(threshold=50)
            
            # Add RMSD stats to summary
            summary_with_rmsd = summary_df.merge(
                rmsd_qc, on='cluster', how='left'
            )
            
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

            # Add average TM-score per cluster (from per-structure stats)
            if 'tm_score' in stats_df.columns:
                cluster_tm = (
                    stats_df.groupby('cluster', as_index=False)['tm_score']
                    .mean()
                    .rename(columns={'tm_score': 'avg_tm_score'})
                )
                summary_with_rmsd = summary_with_rmsd.merge(cluster_tm, on='cluster', how='left')
                log_status("✓ Avg TM-score added to cluster summary")
            
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
            ui.update_select("sequence_viz_cluster", choices=cluster_ids, selected=default_intra)
            
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
        
        log_status("")
        log_status("="*50)
        log_status("🔄 Re-Clustering with new parameters...")
        log_status(f"New min_cluster_size: {min_cluster_size}, New min_samples: {min_samples_hdb}")
        log_status(f"Filter duplicates: {filter_dupes}")
        log_status("="*50)
        
        try:
            # Re-cluster using CACHED Jaccard matrix (instant!)
            clust.distance_matrix = clust.jaccard_distance_matrix
            clust.metric_name = "Jaccard"
            clust.cluster_hdbscan(
                min_cluster_size=min_cluster_size, 
                min_samples=min_samples_hdb,
                filter_duplicates=filter_dupes
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

            # Re-attach TM-score per structure if available
            tm_df = tm_score_stats.get()
            if tm_df is not None and len(tm_df) > 0:
                tm_map = dict(zip(tm_df['structure'], tm_df['tm_score']))
                stats_df['tm_score'] = stats_df['structure'].map(tm_map)
            
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

            # Add average TM-score per cluster
            if 'tm_score' in stats_df.columns:
                cluster_tm = (
                    stats_df.groupby('cluster', as_index=False)['tm_score']
                    .mean()
                    .rename(columns={'tm_score': 'avg_tm_score'})
                )
                summary_with_rmsd = summary_with_rmsd.merge(cluster_tm, on='cluster', how='left')
            
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
            ui.update_select("sequence_viz_cluster", choices=cluster_ids, selected=cluster_ids[0] if cluster_ids else "")
            
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
        """Quality metric cards for interface metrics"""
        if not analysis_complete.get() or clust is None:
            return ui.p("Run analysis to see metrics", style="padding: 20px; text-align: center;")
        
        # For AlphaFold2 models, we could show pLDDT-related metrics here in the future
        # For now, show contact-based metrics
        stats_df = interface_stats.get()
        if stats_df is not None:
            mean_protein_contacts = stats_df['n_protein_residues'].mean()
            mean_nucleic_contacts = stats_df['n_nucleic_residues'].mean()
            max_protein_contacts = stats_df['n_protein_residues'].max()
        else:
            return ui.p("Interface statistics not available", style="padding: 20px; text-align: center; color: #999;")
        
        # 3 interface quality cards
        cards = [
            ui.div(
                ui.h2(f"{mean_protein_contacts:.0f}", class_="text-warning", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Mean Protein Contacts", style="margin: 5px 0; font-size: 0.9em;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            ),
            ui.div(
                ui.h2(f"{max_protein_contacts:.0f}", class_="text-warning", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Max Protein Contacts", style="margin: 5px 0; font-size: 0.9em;"),
                class_="text-center mb-3",
                style="padding: 10px;"
            ),
            ui.div(
                ui.h2(f"{mean_nucleic_contacts:.0f}", class_="text-info", style="font-size: 3em; margin: 10px 0;"),
                ui.p("Mean Nucleic Contacts", style="margin: 5px 0;"),
                class_="text-center",
                style="padding: 10px;"
            )
        ]
        
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
                       'avg_tm_score',
                       'mean_rmsd_sequential', 'medoid_avg_rmsd', 'medoid_structure']
        
        # Only show columns that actually exist
        available_cols = [col for col in display_cols if col in summary_df.columns]
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
        if 'avg_tm_score' in display_df.columns:
            rename_map['avg_tm_score'] = 'Avg TM-score (vs Ref)'
        
        if rename_map:
            display_df = display_df.rename(columns=rename_map)
        
        # Format numeric columns to 4 significant figures
        numeric_cols = ['Avg Motif Match (%)', 'Stability Score',
                       'Intra-Cluster Mean RMSD (Å)', 'Intra-Cluster Std RMSD (Å)', 
                       'Median RMSD (Å)', 'Max RMSD (Å)', 'Min RMSD (Å)',
                       'Sequential RMSD (Å)', 'Medoid Avg RMSD (Å)',
                       'Avg TM-score (vs Ref)']
        
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
                clust.pdb_files
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
    @render.ui
    def sequence_alignment_plot():
        """Sequence alignment QC visualization with secondary structure and pLDDT"""
        if clust is None or not analysis_complete.get():
            return ui.p("Run analysis to view sequence QC", style="padding: 20px; text-align: center;")
        
        try:
            mode = input.sequence_viz_mode()
            
            if mode == "medoids":
                # Medoid comparison mode - get reference sequence if available
                print(f"\n=== Checking for reference sequence ===")
                print(f"  clust object exists: {clust is not None}")
                print(f"  hasattr(clust, 'reference_sequence'): {hasattr(clust, 'reference_sequence') if clust else 'N/A'}")
                print(f"  hasattr(clust, 'reference_pdb'): {hasattr(clust, 'reference_pdb') if clust else 'N/A'}")
                
                ref_seq = None
                ref_pdb = None
                
                if clust and hasattr(clust, 'reference_sequence'):
                    print(f"  clust.reference_sequence type: {type(clust.reference_sequence)}")
                    print(f"  clust.reference_sequence is None: {clust.reference_sequence is None}")
                    if clust.reference_sequence:
                        print(f"  clust.reference_sequence length: {len(clust.reference_sequence)}")
                        print(f"  First 50 chars: {clust.reference_sequence[:50]}")
                        ref_seq = clust.reference_sequence
                    else:
                        print(f"  clust.reference_sequence is None or empty")
                else:
                    print(f"  No reference_sequence attribute on clusterer")
                
                if clust and hasattr(clust, 'reference_pdb'):
                    print(f"  clust.reference_pdb: {clust.reference_pdb}")
                    ref_pdb = clust.reference_pdb
                
                print(f"  ref_seq to pass: {ref_seq is not None} (length={len(ref_seq) if ref_seq else 0})")
                print(f"  ref_pdb to pass: {ref_pdb}")
                
                # Check cluster_summary
                print(f"\n=== Checking cluster_summary ===")
                summary_df = cluster_summary.get()
                print(f"  cluster_summary is None: {summary_df is None}")
                if summary_df is not None:
                    print(f"  cluster_summary shape: {summary_df.shape}")
                    print(f"  cluster_summary columns: {list(summary_df.columns)}")
                    if 'consensus_residues' in summary_df.columns:
                        print(f"  Sample consensus_residues values:")
                        for idx, row in summary_df.head(3).iterrows():
                            print(f"    Cluster {row['cluster']}: '{row['consensus_residues']}'")
                    else:
                        print(f"  ⚠️ WARNING: 'consensus_residues' column NOT FOUND!")
                print(f"===================================\n")
                
                fig = create_alignment_visualization_medoids(
                    clust.pdb_files,
                    clust.labels,
                    reference_pdb=ref_pdb,
                    reference_sequence=ref_seq,
                    cluster_summary=summary_df,
                    protein_chains=clust.protein_chains if hasattr(clust, 'protein_chains') else ['A'],
                    reference_is_alphafold=input.reference_is_alphafold()
                )
            else:
                # Within cluster mode
                cluster_id_str = input.sequence_viz_cluster()
                if not cluster_id_str:
                    return ui.p("Select a cluster", style="padding: 20px; text-align: center;")
                
                cluster_id = int(cluster_id_str)
                
                fig = create_alignment_visualization_within_cluster(
                    clust.pdb_files,
                    clust.labels,
                    cluster_id,
                    protein_chains=clust.protein_chains if hasattr(clust, 'protein_chains') else ['A']
                )
            
            # Convert to HTML with explicit config
            print(f"\n=== FIGURE DEBUG ===")
            print(f"Number of traces: {len(fig.data)}")
            print(f"Figure height: {fig.layout.height}")
            print(f"Subplot rows: {len([d for d in fig.data if hasattr(d, 'xaxis')])}")
            if len(fig.data) > 0:
                print(f"First trace type: {type(fig.data[0]).__name__}")
                print(f"First trace x sample: {fig.data[0].x[:5] if hasattr(fig.data[0], 'x') and len(fig.data[0].x) > 0 else 'N/A'}")
            print(f"====================\n")
            
            # Return just the raw plotly div without extra wrapper
            config = {'displayModeBar': True, 'displaylogo': False}
            html = fig.to_html(
                full_html=False, 
                include_plotlyjs='cdn',
                config=config
            )
            return ui.HTML(html)
            
        except Exception as e:
            import traceback
            error_msg = f"Error creating sequence alignment visualization: {e}\n{traceback.format_exc()}"
            return ui.HTML(f"<pre style='color: red;'>{error_msg}</pre>")
    
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
                    motif_dict = InterfaceAnalyzer.parse_motif_residues(motif_input)
                except:
                    pass
        
        # Start with basic columns (contacts will be added at the end)
        cols = ['structure', 'cluster', 'n_protein_residues', 'n_nucleic_residues', 'representative']

        # Add per-structure TM-score vs reference (if available)
        if 'tm_score' in display_df.columns:
            display_df['TM-score (vs Ref)'] = display_df['tm_score'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "")
            cols.append('TM-score (vs Ref)')
        
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
                
                # Get TM-score if available
                tm_score_display = ""
                if 'tm_score' in struct_info and pd.notna(struct_info['tm_score']):
                    tm_score_val = struct_info['tm_score']
                    tm_score_display = f"<span style='float: right; font-weight: bold; color: #0d6efd;'>TM-score: {tm_score_val:.4g}</span>"
                
                # Get contact residues list (full, no truncation)
                contact_residues = struct_info.get('protein_residues', '')
                contact_residues_html = f"<div style='word-break: break-word; white-space: normal;'><small>Interface: {contact_residues}</small></div>" if contact_residues else "<small>Interface: N/A</small>"
                
                info_html = f"""
                <div style="padding: 10px; background-color: #f8f9fa; margin-bottom: 10px; border-radius: 5px; border-left: 4px solid #0d6efd;">
                    <b>{structure_name}</b> - <span style="color: #0d6efd;">Cluster {struct_info['cluster']}</span>{tm_score_display}<br style="clear: both;">
                    <small>Protein residues: {struct_info['n_protein_residues']} | DNA/RNA residues: {struct_info['n_nucleic_residues']}</small><br>
                    {contact_residues_html}
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
            
            # Use the original PDB file directly - no SASA preprocessing needed
            # AlphaFold2 models already contain pLDDT in B-factor column
            pdb_to_display = pdb_file
            
            # Create 3D viewer with Molstar (display settings handled by Molstar UI)
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
            )
            
            return ui.HTML(viewer_html)
        
        except Exception as e:
            import traceback
            log_status(f"❌ Error loading structure: {e}")
            return ui.p(f"Error loading structure: {e}\n\n{traceback.format_exc()}")


# Helper function: Compute RMSD QC stats
def compute_tm_scores_to_reference(reference_pdb, pdb_files, protein_chains=None):
    """
    Compute TM-score for each model against a provided reference structure.

    Notes
    -----
    - Uses CA atoms from selected protein chains.
    - Matches residues by (chain_instance, resSeq) across reference and target.
    - Returns NaN for structures without enough matched CA atoms.
    """
    import mdtraj as md

    def _get_chain_instance_map(topology):
        counts = {}
        chain_instance = {}
        for chain in topology.chains:
            label = str(chain.chain_id).strip().upper()
            counts[label] = counts.get(label, 0) + 1
            chain_instance[chain.index] = f"{label}#{counts[label]}"
        return chain_instance

    def _resolve_chain_indices(topology, chain_labels):
        if not chain_labels:
            # Auto-select all protein chains
            protein_indices = []
            for chain in topology.chains:
                residues = list(chain.residues)
                if any(getattr(res, 'is_protein', False) for res in residues):
                    protein_indices.append(chain.index)
            return sorted(set(protein_indices))

        requested = {c.strip().upper() for c in chain_labels if c and c.strip()}
        resolved = []
        for chain in topology.chains:
            label = str(chain.chain_id).strip().upper()
            if label in requested:
                resolved.append(chain.index)
        return sorted(set(resolved))

    def _extract_ca_map(traj, chain_indices):
        top = traj.topology
        chain_instance = _get_chain_instance_map(top)
        ca_map = {}
        for atom in top.atoms:
            if atom.name != 'CA':
                continue
            res = atom.residue
            if not getattr(res, 'is_protein', False):
                continue
            if res.chain.index not in chain_indices:
                continue
            key = (chain_instance[res.chain.index], int(res.resSeq))
            if key not in ca_map:
                ca_map[key] = atom.index
        return ca_map

    def _kabsch_align(mobile_xyz, target_xyz):
        mobile_center = mobile_xyz.mean(axis=0)
        target_center = target_xyz.mean(axis=0)

        mobile_centered = mobile_xyz - mobile_center
        target_centered = target_xyz - target_center

        cov = mobile_centered.T @ target_centered
        u, _, vt = np.linalg.svd(cov)
        rot = vt.T @ u.T

        # Ensure right-handed rotation
        if np.linalg.det(rot) < 0:
            vt[-1, :] *= -1
            rot = vt.T @ u.T

        return (mobile_centered @ rot) + target_center

    def _tm_score_from_distances(distances, length_norm):
        if length_norm <= 0:
            return np.nan
        if length_norm <= 21:
            d0 = 0.5
        else:
            d0 = 1.24 * ((length_norm - 15) ** (1.0 / 3.0)) - 1.8
            d0 = max(d0, 0.5)
        return float(np.mean(1.0 / (1.0 + (distances / d0) ** 2)))

    # Load reference
    ref_traj = md.load(reference_pdb)
    ref_chain_indices = _resolve_chain_indices(ref_traj.topology, protein_chains)
    if len(ref_chain_indices) == 0:
        raise ValueError(f"No matching protein chains found in reference for {protein_chains}")

    ref_ca_map = _extract_ca_map(ref_traj, ref_chain_indices)
    ref_ca_count = len(ref_ca_map)
    if ref_ca_count == 0:
        raise ValueError("No reference CA atoms found for selected protein chains")

    results = []
    for pdb_file in pdb_files:
        structure_name = os.path.basename(pdb_file)
        try:
            model_traj = md.load(pdb_file)
            model_chain_indices = _resolve_chain_indices(model_traj.topology, protein_chains)
            model_ca_map = _extract_ca_map(model_traj, model_chain_indices)

            common_keys = sorted(set(ref_ca_map.keys()) & set(model_ca_map.keys()))
            n_aligned = len(common_keys)
            coverage = (n_aligned / ref_ca_count) if ref_ca_count > 0 else np.nan

            if n_aligned < 3:
                results.append({
                    'structure': structure_name,
                    'tm_score': np.nan,
                    'n_aligned_ca': n_aligned,
                    'ca_coverage': coverage
                })
                continue

            ref_idx = np.array([ref_ca_map[k] for k in common_keys], dtype=int)
            model_idx = np.array([model_ca_map[k] for k in common_keys], dtype=int)

            ref_xyz = ref_traj.xyz[0, ref_idx, :]
            model_xyz = model_traj.xyz[0, model_idx, :]

            aligned_model_xyz = _kabsch_align(model_xyz, ref_xyz)
            distances_angstrom = np.linalg.norm(aligned_model_xyz - ref_xyz, axis=1) * 10.0

            tm_score = _tm_score_from_distances(distances_angstrom, length_norm=len(ref_xyz))
            results.append({
                'structure': structure_name,
                'tm_score': tm_score,
                'n_aligned_ca': n_aligned,
                'ca_coverage': coverage
            })

        except Exception:
            results.append({
                'structure': structure_name,
                'tm_score': np.nan,
                'n_aligned_ca': 0,
                'ca_coverage': 0.0
            })

    return pd.DataFrame(results)


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
