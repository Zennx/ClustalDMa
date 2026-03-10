#!/usr/bin/env python3
"""
ClustalDMα - Command Line Interface
Scalable protein-nucleic acid docking analysis for large datasets (1000+ models)

Performs the same analysis as the Shiny app but generates static HTML reports.
"""

import argparse
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime
import json

# Import core modules
from core.clusterer import PDBClusterer
from core.io_utils import parse_alphafold_filename, extract_sequence_from_pdb, read_fasta_sequence, find_pdb_files
from visualization import interactive


def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description='ClustalDMα: Protein-Nucleic Acid Docking Cluster Analysis (CLI)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  python clustal_cli.py -m models/ -o results/
  
  # With reference sequence and custom parameters
  python clustal_cli.py -m models/ -r reference.pdb -o results/ --min-cluster-size 10
  
  # Full parameter specification
  python clustal_cli.py -m models/ -r reference.pdb -o results/ \\
      --ligand-id L --apply-offset --min-cluster-size 5 --min-samples 2 \\
      --filter-duplicates --distance-cutoff 4.5
        """
    )
    
    # Required arguments
    parser.add_argument(
        '-m', '--models',
        required=True,
        help='Directory containing AlphaFold-Multimer models (required)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for results (required)'
    )
    
    # Optional arguments
    parser.add_argument(
        '-r', '--reference',
        default=None,
        help='Reference PDB file for sequence alignment (optional)'
    )
    
    parser.add_argument(
        '--ligand-id',
        default='B',
        help='Ligand chain ID (default: B for second chain)'
    )
    
    parser.add_argument(
        '--apply-offset',
        action='store_true',
        help='Apply residue number corrections for chopped models'
    )
    
    parser.add_argument(
        '--min-cluster-size',
        type=int,
        default=5,
        help='HDBSCAN minimum cluster size (default: 5)'
    )
    
    parser.add_argument(
        '--min-samples',
        type=int,
        default=2,
        help='HDBSCAN minimum samples (default: 2)'
    )
    
    parser.add_argument(
        '--filter-duplicates',
        action='store_true',
        help='Filter near-duplicate structures (RMSD < 0.5 Å)'
    )
    
    parser.add_argument(
        '--distance-cutoff',
        type=float,
        default=4.5,
        help='Contact distance cutoff in Ångströms (default: 4.5)'
    )
    
    parser.add_argument(
        '--n-jobs',
        type=int,
        default=-1,
        help='Number of parallel jobs (-1 for all CPUs, default: -1)'
    )
    
    return parser.parse_args()


def setup_output_directory(output_dir):
    """Create output directory structure"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    subdirs = {
        'html': output_path / 'html_reports',
        'csv': output_path / 'csv_exports',
        'plots': output_path / 'plots',
        'pdb': output_path / 'cluster_representatives'
    }
    
    for subdir in subdirs.values():
        subdir.mkdir(exist_ok=True)
    
    return subdirs


def save_html_report(fig, filepath, title):
    """Save plotly figure as standalone HTML"""
    html = fig.to_html(
        full_html=True,
        include_plotlyjs='cdn',
        config={'displayModeBar': True, 'displaylogo': False}
    )
    
    with open(filepath, 'w') as f:
        f.write(html)
    
    print(f"  ✓ Saved: {filepath.name}")


def create_multi_view_html(figures_dict, output_path, title, default_view=None):
    """
    Create HTML page with dropdown to switch between multiple Plotly figures
    
    Parameters:
    -----------
    figures_dict : dict
        Dictionary mapping view_name -> plotly figure
    output_path : Path
        Where to save the HTML
    title : str
        Page title
    default_view : str, optional
        Which view to show by default (first if None)
    """
    if not figures_dict:
        return
    
    view_names = list(figures_dict.keys())
    if default_view is None:
        default_view = view_names[0]
    
    # Generate HTML for each figure separately with full plot data
    figure_divs = {}
    for name, fig in figures_dict.items():
        # Get the HTML for just the plot div with embedded data
        html_content = fig.to_html(
            full_html=False,
            include_plotlyjs='cdn',
            div_id=f'plot-{name.replace(" ", "-")}'
        )
        figure_divs[name] = html_content
    
    # Build dropdown options
    dropdown_options = ''.join([
        f'<option value="{name}"{" selected" if name == default_view else ""}>{name}</option>'
        for name in view_names
    ])
    
    # Build plot divs (hidden by default)
    plot_divs_html = ''.join([
        f'<div id="view-{name.replace(" ", "-")}" class="plot-view" style="display: {"block" if name == default_view else "none"};">{html}</div>'
        for name, html in figure_divs.items()
    ])
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script src="https://cdn.plot.ly/plotly-2.18.0.min.js"></script>
    <style>
        body {{
            margin: 0;
            padding: 20px;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            background: #f5f5f5;
        }}
        .control-panel {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .control-panel label {{
            font-weight: bold;
            margin-right: 10px;
            color: #333;
        }}
        .control-panel select {{
            padding: 8px 12px;
            font-size: 14px;
            border: 1px solid #ddd;
            border-radius: 4px;
            background: white;
            cursor: pointer;
        }}
        #plot-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .plot-view {{
            width: 100%;
        }}
    </style>
</head>
<body>
    <h1>{title}</h1>
    
    <div class="control-panel">
        <label for="view-selector">Select View:</label>
        <select id="view-selector" onchange="updatePlot()">
            {dropdown_options}
        </select>
    </div>
    
    <div id="plot-container">
        {plot_divs_html}
    </div>
    
    <script>
        function updatePlot() {{
            const selector = document.getElementById('view-selector');
            const selectedView = selector.value;
            
            // Hide all views
            const views = document.querySelectorAll('.plot-view');
            views.forEach(view => view.style.display = 'none');
            
            // Show selected view
            const targetId = 'view-' + selectedView.replace(/ /g, '-');
            const targetView = document.getElementById(targetId);
            if (targetView) {{
                targetView.style.display = 'block';
            }}
        }}
    </script>
</body>
</html>'''
    
    with open(output_path, 'w') as f:
        f.write(html)
    print(f"  ✓ Saved: {output_path.name} ({len(figures_dict)} views)")


def generate_summary_html(output_dirs, params, cluster_summary, analysis_time):
    """Generate main summary HTML page with links to all reports"""
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>ClustalDMα Analysis Report</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            max-width: 1200px;
            margin: 40px auto;
            padding: 20px;
            background: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        .card {{
            background: white;
            padding: 25px;
            margin-bottom: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{ margin: 0; font-size: 2.5em; }}
        h2 {{ color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px; }}
        h3 {{ color: #555; }}
        .param-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .param-item {{
            padding: 15px;
            background: #f8f9fa;
            border-radius: 5px;
            border-left: 4px solid #667eea;
        }}
        .param-label {{
            font-weight: bold;
            color: #666;
            font-size: 0.9em;
            text-transform: uppercase;
        }}
        .param-value {{
            font-size: 1.2em;
            color: #333;
            margin-top: 5px;
        }}
        .cluster-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        .cluster-table th, .cluster-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        .cluster-table th {{
            background: #667eea;
            color: white;
            font-weight: bold;
        }}
        .cluster-table tr:hover {{
            background: #f5f5f5;
        }}
        .report-links {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .report-link {{
            display: block;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            text-decoration: none;
            border-radius: 8px;
            transition: transform 0.2s;
        }}
        .report-link:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }}
        .report-link h3 {{
            margin: 0 0 10px 0;
            color: white;
        }}
        .report-link p {{
            margin: 0;
            opacity: 0.9;
            font-size: 0.9em;
        }}
        .timestamp {{
            color: rgba(255,255,255,0.8);
            font-size: 0.9em;
            margin-top: 10px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 ClustalDMα Analysis Report</h1>
        <div class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
    </div>
    
    <div class="card">
        <h2>Analysis Parameters</h2>
        <div class="param-grid">
            <div class="param-item">
                <div class="param-label">Models Directory</div>
                <div class="param-value">{params['models']}</div>
            </div>
            <div class="param-item">
                <div class="param-label">Total Structures</div>
                <div class="param-value">{params['n_structures']}</div>
            </div>
            <div class="param-item">
                <div class="param-label">Ligand Chain ID</div>
                <div class="param-value">{params['ligand_id']}</div>
            </div>
            <div class="param-item">
                <div class="param-label">Distance Cutoff</div>
                <div class="param-value">{params['distance_cutoff']} Å</div>
            </div>
            <div class="param-item">
                <div class="param-label">Min Cluster Size</div>
                <div class="param-value">{params['min_cluster_size']}</div>
            </div>
            <div class="param-item">
                <div class="param-label">Min Samples</div>
                <div class="param-value">{params['min_samples']}</div>
            </div>
            <div class="param-item">
                <div class="param-label">Apply Residue Offset</div>
                <div class="param-value">{'Yes' if params['apply_offset'] else 'No'}</div>
            </div>
            <div class="param-item">
                <div class="param-label">Filter Duplicates</div>
                <div class="param-value">{'Yes' if params['filter_duplicates'] else 'No'}</div>
            </div>
            <div class="param-item">
                <div class="param-label">Analysis Time</div>
                <div class="param-value">{analysis_time:.1f} seconds</div>
            </div>
        </div>
    </div>
    
    <div class="card">
        <h2>Cluster Summary</h2>
        <table class="cluster-table">
            <thead>
                <tr>
                    <th>Cluster</th>
                    <th>Structures</th>
                    <th>Consensus Residues</th>
                    <th>Binding Mode</th>
                </tr>
            </thead>
            <tbody>
"""
    
    # Add cluster rows
    for _, row in cluster_summary.iterrows():
        if row['cluster'] != -1:
            html_content += f"""
                <tr>
                    <td><strong>Cluster {row['cluster']}</strong></td>
                    <td>{row['n_structures']}</td>
                    <td>{row['n_consensus']}</td>
                    <td>{row['binding_mode']}</td>
                </tr>
"""
    
    html_content += """
            </tbody>
        </table>
    </div>
    
    <div class="card">
        <h2>Interactive Reports</h2>
        <p>Click on any report to view interactive visualizations with dropdown menus:</p>
        <div class="report-links">
            <a href="html_reports/qc_summary.html" class="report-link" style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);">
                <h3>📊 Clustering QC Summary</h3>
                <p>Size distribution, stability, inter-cluster distances</p>
            </a>
            <a href="html_reports/dendrogram.html" class="report-link">
                <h3>🌲 Hierarchical Dendrogram</h3>
                <p>Full distance matrix with hierarchical clustering</p>
            </a>
            <a href="html_reports/dimensionality_reduction.html" class="report-link">
                <h3>🎯 Dimensionality Reduction</h3>
                <p>MDS, t-SNE, PCA, UMAP projections (dropdown menu)</p>
            </a>
            <a href="html_reports/jaccard_heatmap.html" class="report-link">
                <h3>🔥 Jaccard Distance Heatmap</h3>
                <p>Structure-to-structure similarity matrix</p>
            </a>
            <a href="html_reports/intra_cluster_heatmaps.html" class="report-link">
                <h3>📊 Intra-Cluster Heatmaps</h3>
                <p>Within-cluster distances (dropdown per cluster)</p>
            </a>
            <a href="html_reports/hotspots_overall.html" class="report-link">
                <h3>🔥 Overall Binding Hotspots</h3>
                <p>Contact residue frequencies across ALL structures</p>
            </a>
            <a href="html_reports/contact_hotspots_by_cluster.html" class="report-link">
                <h3>🔥 Binding Hotspots by Cluster</h3>
                <p>Contact residue patterns per cluster (dropdown menu)</p>
            </a>
            <a href="html_reports/condensed_tree.html" class="report-link">
                <h3>🌳 HDBSCAN Condensed Tree</h3>
                <p>Cluster formation hierarchy visualization</p>
            </a>
            <a href="html_reports/sequence_alignment.html" class="report-link">
                <h3>🧬 Sequence Alignment QC</h3>
                <p>Secondary structure and pLDDT confidence</p>
            </a>
        </div>
    </div>
    
    <div class="card">
        <h2>Exported Files</h2>
        <h3>CSV Data</h3>
        <ul>
            <li><code>csv_exports/cluster_summary.csv</code> - Cluster statistics</li>
            <li><code>csv_exports/structure_assignments.csv</code> - Structure-to-cluster mapping</li>
            <li><code>csv_exports/jaccard_matrix.csv</code> - Full similarity matrix</li>
            <li><code>csv_exports/contact_exports/</code> - Per-structure contact lists</li>
        </ul>
        
        <h3>Representative Structures</h3>
        <ul>
            <li><code>cluster_representatives/</code> - Medoid PDB files for each cluster</li>
        </ul>
    </div>
</body>
</html>
"""
    
    summary_path = Path(output_dirs['html']).parent / 'index.html'
    with open(summary_path, 'w') as f:
        f.write(html_content)
    
    print(f"\n✅ Main report saved: {summary_path}")
    return summary_path


def main():
    """Main CLI execution"""
    args = parse_arguments()
    
    print("=" * 80)
    print("🧬 ClustalDMα - Large-Scale Docking Analysis CLI")
    print("=" * 80)
    
    # Setup output directories
    print(f"\n📁 Setting up output directory: {args.output}")
    output_dirs = setup_output_directory(args.output)
    
    # Find model files (recursively searches subdirectories)
    print(f"\n🔍 Scanning for models in: {args.models}")
    pdb_files = find_pdb_files(args.models)
    
    if len(pdb_files) == 0:
        print(f"❌ ERROR: No PDB/mmCIF files found in {args.models}")
        print(f"   Searched recursively in all subdirectories.")
        sys.exit(1)
    
    print(f"   Found {len(pdb_files)} PDB/mmCIF files (including subdirectories)")
    
    # Show sample of found files
    if len(pdb_files) <= 5:
        for f in pdb_files:
            print(f"     - {os.path.relpath(f, args.models)}")
    else:
        for f in pdb_files[:3]:
            print(f"     - {os.path.relpath(f, args.models)}")
        print(f"     ... and {len(pdb_files) - 3} more files")
    
    # Load reference sequence if provided
    reference_sequence = None
    reference_pdb = None
    if args.reference:
        print(f"\n📖 Loading reference: {args.reference}")
        reference_pdb = args.reference
        if os.path.exists(reference_pdb):
            # Try to extract sequence from PDB
            if reference_pdb.endswith('.pdb'):
                reference_sequence = extract_sequence_from_pdb(reference_pdb, chain_id='A')
            elif reference_pdb.endswith('.fasta') or reference_pdb.endswith('.fa'):
                reference_sequence = read_fasta_sequence(reference_pdb)
            
            # Ensure we don't pass empty strings
            if reference_sequence and len(reference_sequence) > 0:
                print(f"   Reference sequence length: {len(reference_sequence)}")
            else:
                print(f"   ⚠️ Warning: Could not extract reference sequence")
                reference_sequence = None
    
    # Create clusterer
    print(f"\n⚙️  Initializing clusterer...")
    clusterer = PDBClusterer(pdb_files=[str(f) for f in pdb_files])
    clusterer.reference_sequence = reference_sequence
    
    # Determine residue offset
    residue_offset = 0
    if args.apply_offset:
        print(f"\n🔢 Residue offset correction: ENABLED")
        first_file = os.path.basename(str(pdb_files[0]))
        filename_info = parse_alphafold_filename(first_file)
        if filename_info:
            residue_offset = filename_info['offset']
            print(f"   Detected offset: {residue_offset}")
    else:
        print(f"\n🔢 Residue offset correction: DISABLED")
    
    # Run analysis
    print(f"\n🚀 Starting analysis pipeline...")
    print(f"   HDBSCAN parameters: min_cluster_size={args.min_cluster_size}, min_samples={args.min_samples}")
    print(f"   This may take several minutes for large datasets...")
    
    import time
    start_time = time.time()
    
    try:
        # Use chain-based selections (like Shiny app) for AlphaFold multimers
        # MDTraj uses 0-indexed chainid, so A=0, B=1, etc.
        ligand_chain_idx = ord(args.ligand_id) - ord('A')
        
        # Receptor: all protein atoms NOT in ligand chain
        protein_selection = f"protein and not (chainid == {ligand_chain_idx})"
        # Ligand: just the specified chain (works for both protein and nucleic)
        ligand_selection = f"chainid == {ligand_chain_idx}"
        
        print(f"   Receptor (chain not {args.ligand_id}): '{protein_selection}'")
        print(f"   Ligand (chain {args.ligand_id}): '{ligand_selection}'")
        
        # Compute Jaccard contact matrix
        clusterer.compute_jaccard_contact_matrix(
            distance_cutoff=args.distance_cutoff,
            protein_selection=protein_selection,
            nucleic_selection=ligand_selection,  # Despite name, works for any chain
            apply_offset=args.apply_offset,
            n_jobs=args.n_jobs
        )
        
        # Perform HDBSCAN clustering
        clusterer.cluster_hdbscan(
            min_cluster_size=args.min_cluster_size,
            min_samples=args.min_samples,
            filter_duplicates=args.filter_duplicates
        )
        
        analysis_time = time.time() - start_time
        print(f"\n✅ Analysis complete in {analysis_time:.1f} seconds")
        
        # Extract results
        results = {
            'labels': clusterer.labels,
            'probabilities': getattr(clusterer, 'probabilities', None),
            'jaccard_matrix': clusterer.jaccard_distance_matrix
        }
        
    except Exception as e:
        print(f"\n❌ ERROR during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Get cluster summary
    cluster_summary = clusterer.get_cluster_summary()
    
    print(f"\n📊 Clustering Results:")
    print(f"   Total clusters: {len(cluster_summary[cluster_summary['cluster'] != -1])}")
    noise_count = cluster_summary[cluster_summary['cluster'] == -1]['n_structures'].sum() if -1 in cluster_summary['cluster'].values else 0
    print(f"   Noise points: {noise_count}")
    print(f"   Largest cluster: {cluster_summary[cluster_summary['cluster'] != -1]['n_structures'].max()} structures")
    
    # Save cluster summary CSV
    csv_path = output_dirs['csv'] / 'cluster_summary.csv'
    cluster_summary.to_csv(csv_path, index=False)
    print(f"\n💾 Saved cluster summary: {csv_path.name}")
    
    # Save structure assignments
    structures = [os.path.basename(str(f)) for f in pdb_files]
    labels = results['labels']
    probs = results.get('probabilities')

    print("pdb_files:", len(structures))
    print("labels:", len(labels))
    probs = results.get('probabilities')
    print("probabilities:", len(probs) if probs is not None else 0)

    if probs is None or len(probs) != len(labels):
       probs = [1.0] * len(labels)

    structures = structures[:len(labels)]

    assignments_df = pd.DataFrame({
        'structure': [os.path.basename(str(f)) for f in pdb_files],
        'cluster': results['labels'],
        'probability': results.get('probabilities', [1.0] * len(pdb_files))
    })
    assignments_path = output_dirs['csv'] / 'structure_assignments.csv'
    assignments_df.to_csv(assignments_path, index=False)
    print(f"💾 Saved structure assignments: {assignments_path.name}")
    
    # Save Jaccard matrix
    if 'jaccard_matrix' in results:
        jaccard_path = output_dirs['csv'] / 'jaccard_matrix.csv'
        np.savetxt(jaccard_path, results['jaccard_matrix'], delimiter=',')
        print(f"💾 Saved Jaccard matrix: {jaccard_path.name}")
    
    # Generate visualizations
    print(f"\n🎨 Generating interactive visualizations...")
    
    # 1. Hierarchical dendrogram
    print("   • Dendrogram...")
    try:
        fig_dend = interactive.create_distance_heatmap_with_dendrogram(
            clusterer.jaccard_distance_matrix,
            results['labels'],
            metric_name='Jaccard Distance',
            pdb_files=[os.path.basename(str(f)) for f in pdb_files]
        )
        save_html_report(fig_dend, output_dirs['html'] / 'dendrogram.html', 'Dendrogram')
    except Exception as e:
        print(f"     ⚠️ Warning: Could not generate dendrogram: {e}")
        import traceback
        traceback.print_exc()
    
    # 2. Dimensionality Reduction - All methods in one page with dropdown
    print("   • Dimensionality Reduction Views (MDS, t-SNE, PCA, UMAP)...")
    try:
        df_for_scatter = pd.DataFrame({
            'structure': [os.path.basename(str(f)) for f in pdb_files],
            'cluster': results['labels'],
            'n_protein_residues': [0] * len(pdb_files),
            'protein_residues': [''] * len(pdb_files)
        })
        
        dimred_figures = {}
        for method in ['mds', 'tsne', 'pca', 'umap']:
            try:
                fig = interactive.create_scatter_multimethod(
                    df_for_scatter,
                    clusterer.jaccard_distance_matrix,
                    method=method
                )
                dimred_figures[method.upper()] = fig
            except Exception as e:
                if method != 'umap':  # UMAP optional
                    print(f"       ⚠️ {method.upper()} failed: {e}")
        
        if dimred_figures:
            create_multi_view_html(
                dimred_figures,
                output_dirs['html'] / 'dimensionality_reduction.html',
                '🎯 Cluster Visualization - Dimensionality Reduction',
                default_view='MDS'
            )
    except Exception as e:
        print(f"     ⚠️ Warning: Could not generate dim-red views: {e}")
        import traceback
        traceback.print_exc()
    
    # 3. Jaccard heatmap
    print("   • Jaccard Heatmap...")
    try:
        fig_jaccard = interactive.create_distance_heatmap(
            clusterer.jaccard_distance_matrix,
            results['labels'],
            metric_name='Jaccard Distance'
        )
        save_html_report(fig_jaccard, output_dirs['html'] / 'jaccard_heatmap.html', 'Jaccard Heatmap')
    except Exception as e:
        print(f"     ⚠️ Warning: Could not generate heatmap: {e}")
        import traceback
        traceback.print_exc()
    
    # 4. Condensed tree
    print("   • HDBSCAN Condensed Tree...")
    try:
        html_tree = interactive.create_hdbscan_condensed_tree(
            clusterer,
            results['labels'],
            pdb_files=[os.path.basename(str(f)) for f in pdb_files]
        )
        tree_path = output_dirs['html'] / 'condensed_tree.html'
        with open(tree_path, 'w') as f:
            f.write(f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>HDBSCAN Condensed Tree</title>
    <style>
        body {{ margin: 0; padding: 20px; font-family: Arial, sans-serif; }}
    </style>
</head>
<body>
    <h1>HDBSCAN Condensed Tree</h1>
    {html_tree}
</body>
</html>''')
        print(f"  ✓ Saved: condensed_tree.html")
    except Exception as e:
        print(f"     ⚠️ Warning: Could not generate condensed tree: {e}")
        import traceback
        traceback.print_exc()
    
    # 5. Intra-cluster heatmaps - All in one page with dropdown
    print("   • Intra-Cluster Distance Heatmaps...")
    try:
        intra_matrices, cluster_members = interactive.compute_intra_cluster_matrices(
            clusterer.jaccard_distance_matrix,
            results['labels']
        )
        
        intra_figures = {}
        for cluster_id, members in cluster_members.items():
            if cluster_id == -1 or len(members) <= 1:
                continue
            
            fig = interactive.create_single_intra_cluster_heatmap(
                clusterer.jaccard_distance_matrix,
                results['labels'],
                cluster_id,
                pdb_files=[os.path.basename(str(f)) for f in pdb_files],
                metric_name='Jaccard Distance'
            )
            intra_figures[f"Cluster {cluster_id} (n={len(members)})"] = fig
        
        if intra_figures:
            create_multi_view_html(
                intra_figures,
                output_dirs['html'] / 'intra_cluster_heatmaps.html',
                '📊 Intra-Cluster Distance Heatmaps'
            )
            print(f"     ✓ Generated heatmaps for {len(intra_figures)} clusters")
    except Exception as e:
        print(f"     ⚠️ Warning: Could not generate intra-cluster heatmaps: {e}")
        import traceback
        traceback.print_exc()
    
    # 6. Contact/Hotspot Analysis
    print("   • Contact Residue Hotspots...")
    try:
        # First, generate OVERALL hotspots (all structures combined, like Shiny app)
        overall_hotspots = clusterer.get_binding_hotspots()
        
        if not overall_hotspots.empty:
            fig_overall = interactive.create_hotspot_histogram(
                overall_hotspots,
                top_n=None,
                smooth=True,
                split_by_chain=False
            )
            fig_overall.update_layout(
                title=f'Overall Binding Hotspots - All Structures (n={len(pdb_files)})'
            )
            
            # Save as standalone HTML
            overall_html_path = output_dirs['html'] / 'hotspots_overall.html'
            with open(overall_html_path, 'w') as f:
                f.write(fig_overall.to_html(
                    full_html=True,
                    include_plotlyjs='cdn'
                ))
            print(f"     ✓ Generated overall hotspots across {len(pdb_files)} structures")
        
        # Then, generate per-cluster hotspots + contact heatmaps (as separate blocks)
        hotspot_figures = {}
        unique_clusters = sorted([c for c in set(results['labels']) if c != -1])
        
        for cluster_id in unique_clusters:
            # Get hotspots for this cluster using the proper method
            cluster_hotspots = clusterer.get_binding_hotspots(cluster_id=cluster_id)
            
            if not cluster_hotspots.empty:
                n_structures = sum(1 for l in results['labels'] if l == cluster_id)
                
                # Create histogram
                fig_histogram = interactive.create_hotspot_histogram(
                    cluster_hotspots,
                    top_n=None,
                    smooth=True,
                    split_by_chain=False
                )
                fig_histogram.update_layout(
                    title=f'Cluster {cluster_id} Binding Frequency Pattern ({n_structures} structures)',
                    height=400
                )
                
                # Create contact distance heatmap
                fig_heatmap = interactive.create_contact_residue_heatmap(
                    clusterer.contact_residue_pairs,
                    results['labels'],
                    pdb_files,
                    cluster_id=cluster_id
                )
                
                # Generate separate HTML blocks for each figure
                histogram_html = fig_histogram.to_html(
                    full_html=False,
                    include_plotlyjs='cdn',
                    div_id=f'histogram-cluster-{cluster_id}'
                )
                
                heatmap_html = fig_heatmap.to_html(
                    full_html=False,
                    include_plotlyjs='cdn',
                    div_id=f'heatmap-cluster-{cluster_id}'
                )
                
                # Combine as separate blocks with spacing
                combined_html = f'''
                <div style="margin-bottom: 30px;">
                    {histogram_html}
                </div>
                <div style="margin-bottom: 30px;">
                    {heatmap_html}
                </div>
                '''
                
                # Store the combined HTML string instead of a figure
                hotspot_figures[f"Cluster {cluster_id} (n={n_structures})"] = combined_html
        
        if hotspot_figures:
            # Create custom multi-view HTML that handles HTML strings instead of figures
            view_names = list(hotspot_figures.keys())
            default_view = view_names[0]
            
            dropdown_options = ''.join([
                f'<option value="{name}"{" selected" if name == default_view else ""}>{name}</option>'
                for name in view_names
            ])
            
            plot_divs_html = ''.join([
                f'<div id="view-{name.replace(" ", "-")}" class="plot-view" style="display: {"block" if name == default_view else "none"};">{html}</div>'
                for name, html in hotspot_figures.items()
            ])
            
            full_html = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>🔥 Contact Residue Hotspots by Cluster</title>
    <script src="https://cdn.plot.ly/plotly-2.18.0.min.js"></script>
    <style>
        body {{
            margin: 0;
            padding: 20px;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            background: #f5f5f5;
        }}
        .control-panel {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .control-panel label {{
            font-weight: bold;
            margin-right: 10px;
            color: #333;
        }}
        .control-panel select {{
            padding: 8px 12px;
            font-size: 14px;
            border: 1px solid #ddd;
            border-radius: 4px;
            background: white;
            cursor: pointer;
        }}
        #plot-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .plot-view {{
            width: 100%;
        }}
    </style>
</head>
<body>
    <h1>🔥 Contact Residue Hotspots by Cluster</h1>
    
    <div class="control-panel">
        <label for="view-selector">Select View:</label>
        <select id="view-selector" onchange="updatePlot()">
            {dropdown_options}
        </select>
    </div>
    
    <div id="plot-container">
        {plot_divs_html}
    </div>
    
    <script>
        function updatePlot() {{
            const selector = document.getElementById('view-selector');
            const selectedView = selector.value;
            
            // Hide all views
            const views = document.querySelectorAll('.plot-view');
            views.forEach(view => view.style.display = 'none');
            
            // Show selected view
            const targetId = 'view-' + selectedView.replace(/ /g, '-');
            const targetView = document.getElementById(targetId);
            if (targetView) {{
                targetView.style.display = 'block';
            }}
        }}
    </script>
</body>
</html>'''
            
            output_path = output_dirs['html'] / 'contact_hotspots_by_cluster.html'
            with open(output_path, 'w') as f:
                f.write(full_html)
            
            print(f"     ✓ Generated per-cluster hotspots + heatmaps for {len(hotspot_figures)} clusters")
    except Exception as e:
        print(f"     ⚠️ Warning: Could not generate hotspots: {e}")
        import traceback
        traceback.print_exc()
    
    # 7. QC Summary Page
    print("   • QC Summary Report...")
    try:
        # Generate QC plots
        qc_figures = {}
        
        # Cluster size distribution
        fig_sizes = interactive.create_cluster_size_distribution(results['labels'])
        qc_figures["Cluster Sizes"] = fig_sizes
        
        # Cluster stability (if HDBSCAN)
        if hasattr(clusterer, 'hdbscan_clusterer'):
            fig_stability = interactive.create_cluster_stability_plot(
                clusterer.hdbscan_clusterer,
                results['labels'],
                distance_matrix=clusterer.jaccard_distance_matrix
            )
            qc_figures["Cluster Stability"] = fig_stability
        
        # Inter-cluster distance matrix
        cluster_matrix, cluster_ids, cluster_sizes = interactive.compute_cluster_distance_matrix(
            clusterer.jaccard_distance_matrix,
            results['labels']
        )
        fig_inter = interactive.create_cluster_distance_heatmap(
            cluster_matrix,
            cluster_ids,
            cluster_sizes,
            metric_name='Jaccard Distance'
        )
        qc_figures["Inter-Cluster Distances"] = fig_inter
        
        if qc_figures:
            create_multi_view_html(
                qc_figures,
                output_dirs['html'] / 'qc_summary.html',
                '📊 Clustering QC Summary'
            )
            print(f"     ✓ Generated QC summary with {len(qc_figures)} plots")
    except Exception as e:
        print(f"     ⚠️ Warning: Could not generate QC summary: {e}")
        import traceback
        traceback.print_exc()
    
    # 8. Sequence alignment QC
    if reference_sequence:
        print("   • Sequence Alignment QC...")
        try:
            from visualization.sequence_alignment import create_alignment_visualization_medoids
            fig_align = create_alignment_visualization_medoids(
                pdb_files=[str(f) for f in pdb_files],
                labels=results['labels'],
                reference_pdb=reference_pdb,
                reference_sequence=reference_sequence,
                cluster_summary=cluster_summary
            )
            save_html_report(fig_align, output_dirs['html'] / 'sequence_alignment.html', 'Sequence Alignment')
        except Exception as e:
            print(f"     ⚠️ Warning: Could not generate sequence alignment: {e}")
    
    # Export contact lists
    print(f"\n💾 Exporting contact lists...")
    contact_export_dir = output_dirs['csv'] / 'contact_exports'
    contact_export_dir.mkdir(exist_ok=True)
    
    if hasattr(clusterer, 'contact_sets_'):
        for i, (pdb_file, contacts) in enumerate(zip(pdb_files, clusterer.contact_lists_)):
            if len(contacts) > 0:
                basename = os.path.basename(str(pdb_file)).replace('.pdb', '')
                contact_df = pd.DataFrame(contacts)
                contact_path = contact_export_dir / f'{basename}_contacts.csv'
                contact_df.to_csv(contact_path, index=False)
        print(f"   Exported {len(clusterer.contact_lists_)} contact lists")
    
    # Copy representative structures
    print(f"\n📦 Copying cluster representatives...")
    import shutil
    for _, row in cluster_summary.iterrows():
        if row['cluster'] != -1:
            cluster_id = row['cluster']
            cluster_indices = [i for i, lbl in enumerate(results['labels']) if lbl == cluster_id]
            if cluster_indices:
                medoid_idx = cluster_indices[0]
                source_pdb = pdb_files[medoid_idx]
                dest_pdb = output_dirs['pdb'] / f"cluster_{cluster_id}_medoid.pdb"
                shutil.copy2(source_pdb, dest_pdb)
    print(f"   Copied {len(cluster_summary[cluster_summary['cluster'] != -1])} representative structures")
    
    # Generate main summary HTML
    print(f"\n📄 Generating summary report...")
    summary_path = generate_summary_html(
        output_dirs,
        {
            'models': args.models,
            'n_structures': len(pdb_files),
            'ligand_id': args.ligand_id,
            'distance_cutoff': args.distance_cutoff,
            'min_cluster_size': args.min_cluster_size,
            'min_samples': args.min_samples,
            'apply_offset': args.apply_offset,
            'filter_duplicates': args.filter_duplicates
        },
        cluster_summary,
        analysis_time
    )
    
    print("\n" + "=" * 80)
    print("✅ ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\n📊 Open the main report in your browser:")
    print(f"   file://{summary_path.absolute()}")
    print(f"\n📁 All outputs saved to: {args.output}/")
    print("=" * 80 + "\n")


if __name__ == '__main__':
    main()
