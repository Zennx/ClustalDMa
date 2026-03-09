"""
Interactive visualization functions using Plotly and Molstar
"""

import os
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import json
import base64
import re
from pathlib import Path
import pandas as pd
import plotly.figure_factory as ff
#import py3Dmol
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE

def create_interactive_scatter(df, distance_matrix):
    """
    Create interactive scatter plot with structure preview on hover
    
    Parameters:
    -----------
    df : pd.DataFrame
        Interface statistics DataFrame
    distance_matrix : np.ndarray
        Distance matrix for PCA projection
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    # Reduce to 2D using PCA (MDS on distance matrix)
    from sklearn.manifold import MDS
    
    if distance_matrix.shape[0] > 2:
        # Use MDS for distance matrix visualization (better than PCA for distances)
        # print(f"Computing MDS projection for {distance_matrix.shape[0]} structures...")
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, n_init=1, max_iter=300)
        coords_2d = mds.fit_transform(distance_matrix)
        # print(f"✓ MDS complete. Coords range: X=[{coords_2d[:, 0].min():.2f}, {coords_2d[:, 0].max():.2f}], Y=[{coords_2d[:, 1].min():.2f}, {coords_2d[:, 1].max():.2f}]")
    else:
        # Not enough points, just use first two dimensions
        coords_2d = np.column_stack([np.arange(distance_matrix.shape[0]), np.zeros(distance_matrix.shape[0])])
        # print(f"✓ Too few points for MDS, using simple layout")
    
    # Create hover text with interface info
    hover_text = []
    for idx, row in df.iterrows():
        text = f"<b>{row['structure']}</b><br>"
        text += f"Cluster: {row['cluster']}<br>"
        text += f"Protein residues: {row['n_protein_residues']}<br>"
        if len(row['protein_residues']) > 50:
            text += f"Interface: {row['protein_residues'][:47]}..."
        else:
            text += f"Interface: {row['protein_residues']}"
        hover_text.append(text)
    
    # Create figure
    fig = go.Figure()
    
    # Add scatter points colored by cluster
    colors = df['cluster'].values
    # print(f"Creating scatter plot with {len(colors)} points, clusters: {set(colors)}")
    
    fig.add_trace(go.Scatter(
        x=coords_2d[:, 0],
        y=coords_2d[:, 1],
        mode='markers',
        marker=dict(
            size=12,
            color=colors,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Cluster ID"),
            line=dict(width=1, color='white')
        ),
        text=hover_text,
        hoverinfo='text',
        customdata=df['structure'].values,
    ))
    
    # print(f"✓ Scatter plot created successfully")
    
    fig.update_layout(
        title="Interactive Structure Clustering Map (MDS Projection)",
        xaxis_title="Dimension 1",
        yaxis_title="Dimension 2",
        height=600,
        hovermode='closest',
        template='plotly', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
        yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0')
    )
    
    return fig


def create_hotspot_histogram(hotspots_df, top_n=None, smooth=True, split_by_chain=False):
    """
    Create binding hotspot histogram - genome browser style
    
    Parameters:
    -----------
    hotspots_df : pd.DataFrame
        Hotspot data from InterfaceAnalyzer
    top_n : int
        Number of top hotspots to display (None = all)
    smooth : bool
        If True, use bar chart; if False, use gradient bars
    split_by_chain : bool
        If True, return dict of figures by chain; if False, return single combined figure
    
    Returns:
    --------
    plotly.graph_objects.Figure or dict of Figures (if split_by_chain=True)
    """
    import re
    
    # Extract residue numbers for sorting
    def extract_resnum(res_str):
        # Handle both "ARG123" and "A:ARG123" formats
        match = re.search(r'(\d+)', res_str)
        return int(match.group(1)) if match else 0
    
    # Extract chain ID for grouping
    def extract_chain(res_str):
        if ':' in res_str:
            return res_str.split(':')[0]
        return 'A'  # Default chain if not specified
    
    # Add residue number and chain columns
    hotspots_df = hotspots_df.copy()
    hotspots_df['resnum'] = hotspots_df['residue'].apply(extract_resnum)
    hotspots_df['chain'] = hotspots_df['residue'].apply(extract_chain)
    
    # Sort by chain, then residue number - CRITICAL: reset_index to avoid misalignment
    hotspots_df = hotspots_df.sort_values(['chain', 'resnum']).reset_index(drop=True)
    
    # Optionally limit to top N
    if top_n is not None:
        top_residues = hotspots_df.nlargest(top_n, 'percentage')
        hotspots_df = top_residues.sort_values(['chain', 'resnum']).reset_index(drop=True)
    
    # Color scheme for chains
    chain_colors = {
        'A': 'rgb(220, 53, 69)',   # Red
        'B': 'rgb(13, 110, 253)',   # Blue
        'C': 'rgb(25, 135, 84)',    # Green
        'D': 'rgb(255, 193, 7)',    # Yellow/Orange
        'E': 'rgb(111, 66, 193)',   # Purple
        'F': 'rgb(13, 202, 240)',   # Cyan
    }
    
    if split_by_chain:
        # Create separate figures for each chain
        figures = {}
        chains = sorted(hotspots_df['chain'].unique())
        
        for chain in chains:
            chain_data = hotspots_df[hotspots_df['chain'] == chain]
            
            if len(chain_data) == 0:
                continue
            
            color = chain_colors.get(chain, 'rgb(108, 117, 125)')  # Gray for others
            
            fig = go.Figure()
            
            fig.add_trace(go.Bar(
                x=chain_data['resnum'],
                y=chain_data['percentage'],
                marker=dict(
                    color=color,
                    line=dict(width=0)
                ),
                hovertemplate='<b>%{text}</b><br>Position: %{x}<br>Count: %{customdata[0]}<br>Frequency: %{y:.1f}%<extra></extra>',
                text=chain_data['residue'],
                customdata=chain_data[['frequency']].values,
                name=f'Chain {chain}'
            ))
            
            fig.update_layout(
                title=f'Chain {chain} Binding Hotspots',
                xaxis_title='Residue Position',
                yaxis_title='Binding Frequency (%)',
                height=350,
                template='plotly', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
                hovermode='closest',
                bargap=0.1,
                xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
                yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0')
            )
            
            figures[chain] = fig
        
        return figures
    
    else:
        # Single combined figure with color by chain
        fig = go.Figure()
        
        # Color each bar by its chain
        colors = [chain_colors.get(chain, 'rgb(108, 117, 125)') for chain in hotspots_df['chain']]
        
        fig.add_trace(go.Bar(
            x=hotspots_df['resnum'],
            y=hotspots_df['percentage'],
            marker=dict(
                color=colors,
                line=dict(width=0)
            ),
            hovertemplate='<b>%{text}</b><br>Position: %{x}<br>Count: %{customdata[0]}<br>Frequency: %{y:.1f}%<extra></extra>',
            text=hotspots_df['residue'],
            customdata=hotspots_df[['frequency']].values,
            name='Binding Frequency'
        ))
        
        fig.update_layout(
            title='Protein Binding Hotspots (All Chains)',
            xaxis_title='Residue Position',
            yaxis_title='Binding Frequency (%)',
            height=400,
            template='plotly', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
            hovermode='closest',
            bargap=0.1,
            xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
            yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0')
        )
        
        return fig


def create_molstar_viewer_html(pdb_file, width=800, height=600, 
                               protein_style='surface', protein_color='cyan',
                               dna_style='stick', show_surface=True, 
                               surface_opacity=0.5, zoom_level=1.0, 
                               protein_chains=None, nucleic_chains=None):
    """
    Create Molstar HTML for structure viewing
    Uses molviewspec Python library for proper state generation
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file
    width : int
        Viewer width in pixels
    height : int
        Viewer height in pixels
    protein_style : str
        Style for protein: 'cartoon', 'stick', 'sphere', 'line'
    protein_color : str
        Color for protein representation
    dna_style : str
        Style for DNA: 'stick', 'sphere', 'line', 'cartoon'
    show_surface : bool
        Show protein surface
    surface_opacity : float
        Surface opacity (0-1)
    zoom_level : float
        Zoom multiplier (1.0 = default)
    protein_chains : list, optional
        List of chain IDs that are protein (e.g., ['A', 'B'])
    nucleic_chains : list, optional
        List of chain IDs that are nucleic acid (e.g., ['C', 'D'])
    
    Returns:
    --------
    str : HTML string for embedding
    """
    import base64
    import os
    import json
    
    # Read PDB file and encode
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()
    pdb_base64 = base64.b64encode(pdb_data.encode()).decode()
    
    print(f"Building Molstar viewer: cartoon with pLDDT coloring")
    
    
    # Build HTML with Molstar and custom preset
    html_template = f'''
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8" />
    <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css" />
    <style>
        /* Dark mode override for Molstar */
        @media (prefers-color-scheme: dark) {{
            #molstar-container {{
                background-color: #1a1a1a !important;
            }}
            .msp-plugin {{
                background-color: #1a1a1a !important;
                color: #e0e0e0 !important;
            }}
        }}
    </style>
    <style>
        #molstar-container {{
            position: relative;
            width: {width}px;
            height: {height}px;
            margin: 0 auto;
        }}
    </style>
</head>
<body>
    <div id="molstar-container"></div>
    
    <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js"></script>
    <script type="text/javascript">
        // Wait for both DOM and Molstar to load
        function initMolstar() {{
            // Check DOM is ready
            if (!document.body) {{
                console.log('Waiting for DOM...');
                setTimeout(initMolstar, 50);
                return;
            }}
            
            // Check Molstar is loaded
            if (typeof molstar === 'undefined') {{
                console.log('Waiting for Molstar to load...');
                setTimeout(initMolstar, 100);
                return;
            }}
            
            startVisualization();
        }}
        
        async function startVisualization() {{
            const container = document.getElementById('molstar-container');
            
            console.log('=== Initializing Molstar ===');
            console.log('Molstar object:', molstar);
            console.log('Available methods:', Object.keys(molstar));
            
            // Initialize Molstar plugin using PluginContext
            const plugin = await molstar.Viewer.create(container, {{
                layoutIsExpanded: false,
                layoutShowControls: false,
                layoutShowRemoteState: false,
                layoutShowSequence: true,
                layoutShowLog: false,
                layoutShowLeftPanel: false,
                viewportShowExpand: true,
                viewportShowSelectionMode: false,
                viewportShowAnimation: false
            }});
            
            console.log('✓ Molstar viewer initialized');
            console.log('Plugin wrapper:', plugin);
            console.log('Plugin wrapper properties:', Object.keys(plugin));
            
            // Access the actual plugin from the wrapper
            const actualPlugin = plugin.plugin;
            console.log('Actual plugin:', actualPlugin);
            console.log('Actual plugin properties:', Object.keys(actualPlugin));
            
            if (!actualPlugin) {{
                console.error('✗ Actual plugin is null/undefined');
                return;
            }}
            
            try {{
                // Load structure from base64 using proper Molstar API
                const pdbData = atob('{pdb_base64}');
                
                console.log('Loading PDB structure...');
                
                // Use the plugin's data manager to load PDB data
                const data = await actualPlugin.builders.data.rawData({{
                    data: pdbData,
                    label: 'Structure'
                }});
                
                const trajectory = await actualPlugin.builders.structure.parseTrajectory(data, 'pdb');
                const model = await actualPlugin.builders.structure.createModel(trajectory);
                const structure = await actualPlugin.builders.structure.createStructure(model);
                
                console.log('✓ Structure loaded');
                console.log('Structure object:', structure);
                
                // Clear any default representations
                const state = actualPlugin.state.data;
                const reprs = state.selectQ(q => q.ofType('representation'));
                console.log(`Clearing ${{reprs.length}} default representations...`);
                for (const repr of reprs) {{
                    await actualPlugin.build().delete(repr).commit();
                }}
                
                // Create protein component
                console.log('Creating protein component...');
                const proteinComp = await actualPlugin.builders.structure.tryCreateComponentStatic(
                    structure, 'protein'
                );
                
                if (proteinComp) {{
                    console.log('✓ Protein component created');
                    
                    // Add cartoon with uncertainty coloring (pLDDT)
                    // Following MolStar's pLDDT colouring, same as AlphaFold server/ AlphaFoldDB
                    await actualPlugin.builders.structure.representation.addRepresentation(
                        proteinComp,
                        {{
                            type: 'cartoon',
                            typeParams: {{}},
                            color: 'plddt-confidence'
                        }}  

                    );
                    console.log('✓ Cartoon representation added with (pLDDT) coloring');
                }} else {{
                    console.error('✗ Failed to create protein component');
                }}
                
                // Create nucleic acid component
                console.log('Creating nucleic acid component...');
                const nucleicComp = await actualPlugin.builders.structure.tryCreateComponentStatic(
                    structure, 'nucleic'
                );
                
                if (nucleicComp) {{
                    console.log('✓ Nucleic acid component created');
                    
                    console.log('Adding nucleic acid ball-and-stick representation');
                    await actualPlugin.builders.structure.representation.addRepresentation(
                        nucleicComp,
                        {{
                            type: 'ball-and-stick',
                            typeParams: {{ sizeFactor: 0.3 }},
                            color: 'element-symbol'
                        }}
                    );
                    console.log('✓ Nucleic acid representation added');
                }} else {{
                    console.warn('No nucleic acid found (this is OK for protein-only structures)');
                }}
                
                // Reset camera and apply zoom
                console.log('Resetting camera...');
                actualPlugin.canvas3d.requestCameraReset();
                setTimeout(() => {{
                    if (actualPlugin.canvas3d?.camera?.state) {{
                        const currentRadius = actualPlugin.canvas3d.camera.state.radius;
                        actualPlugin.canvas3d.camera.setState({{
                            radius: currentRadius / {zoom_level}
                        }}, 200);
                    }}
                    console.log('✓ Camera adjusted');
                }}, 500);
                
                console.log('=== Visualization complete ===');
                
            }} catch (error) {{
                console.error('✗ Error:', error.message);
                console.error('Stack:', error.stack);
            }}
        }}
        
        // Start when page loads
        if (document.readyState === 'loading') {{
            document.addEventListener('DOMContentLoaded', initMolstar);
        }} else {{
            initMolstar();
        }}
    </script>
</body>
</html>
'''
    
    return html_template


def create_scatter_multimethod(df, distance_matrix, method='mds'):
    """
    Create interactive scatter plot with multiple dimensionality reduction methods
    
    Parameters:
    -----------
    df : pd.DataFrame
        Interface statistics DataFrame
    distance_matrix : np.ndarray
        Distance matrix for projection
    method : str
        Dimensionality reduction method: 'mds', 'tsne', or 'pca'
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    if distance_matrix.shape[0] <= 2:
        coords_2d = np.column_stack([np.arange(distance_matrix.shape[0]), np.zeros(distance_matrix.shape[0])])
    elif method == 'mds':
        # print(f"Computing MDS projection...")
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, n_init=1, max_iter=300)
        coords_2d = mds.fit_transform(distance_matrix)
    elif method == 'tsne':
        # print(f"Computing t-SNE projection...")
        # t-SNE with precomputed metric requires init='random' (not 'pca')
        tsne = TSNE(n_components=2, metric='precomputed', random_state=42, 
                    perplexity=min(30, distance_matrix.shape[0]-1),
                    init='random')  # FIX: Must use 'random' init with precomputed metric
        coords_2d = tsne.fit_transform(distance_matrix)
    elif method == 'pca':
        # print(f"Computing PCA projection (on distance matrix)...")
        pca = PCA(n_components=2, random_state=42)
        coords_2d = pca.fit_transform(distance_matrix)
    elif method == 'umap':
        try:
            from umap import UMAP
            # print(f"Computing UMAP projection...")
            reducer = UMAP(n_components=2, metric='precomputed', random_state=42, n_neighbors=min(15, distance_matrix.shape[0]-1))
            coords_2d = reducer.fit_transform(distance_matrix)
        except ImportError:
            raise ImportError("UMAP not installed. Install with: pip install umap-learn")
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # print(f"✓ {method.upper()} complete. Coords range: X=[{coords_2d[:, 0].min():.2f}, {coords_2d[:, 0].max():.2f}], Y=[{coords_2d[:, 1].min():.2f}, {coords_2d[:, 1].max():.2f}]")
    
    # Create hover text
    hover_text = []
    for idx, row in df.iterrows():
        text = f"<b>{row['structure']}</b><br>"
        text += f"Cluster: {row['cluster']}<br>"
        text += f"Protein residues: {row['n_protein_residues']}<br>"
        if len(row['protein_residues']) > 50:
            text += f"Interface: {row['protein_residues'][:47]}..."
        else:
            text += f"Interface: {row['protein_residues']}"
        hover_text.append(text)
    
    # Color map for clusters
    cluster_colors = {}
    unique_clusters = sorted(df['cluster'].unique())
    color_palette = px.colors.qualitative.Plotly
    for i, cluster_id in enumerate(unique_clusters):
        if cluster_id == -1:
            cluster_colors[cluster_id] = '#888888'  # Gray for noise
        else:
            cluster_colors[cluster_id] = color_palette[i % len(color_palette)]
    
    df['color'] = df['cluster'].map(cluster_colors)
    
    # Create figure
    fig = go.Figure()
    
    for cluster_id in unique_clusters:
        mask = df['cluster'] == cluster_id
        cluster_df = df[mask]
        cluster_coords = coords_2d[mask]
        cluster_hover = [h for h, m in zip(hover_text, mask) if m]
        
        label = f"Cluster {cluster_id}" if cluster_id != -1 else "Noise"
        
        fig.add_trace(go.Scatter(
            x=cluster_coords[:, 0],
            y=cluster_coords[:, 1],
            mode='markers',
            name=label,
            marker=dict(
                size=10,
                color=cluster_colors[cluster_id],
                line=dict(width=1, color='white')
            ),
            text=cluster_hover,
            hovertemplate='%{text}<extra></extra>'
        ))
    
    method_name = {'mds': 'MDS', 'tsne': 't-SNE', 'pca': 'PCA', 'umap': 'UMAP'}[method]
    fig.update_layout(
        title=f'Cluster Visualization ({method_name} Projection)',
        xaxis_title=f'{method_name} Component 1',
        yaxis_title=f'{method_name} Component 2',
        hovermode='closest',
        showlegend=True,
        height=600,
        template='plotly', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
        yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0')
    )
    
    return fig


def compute_cluster_distance_matrix(distance_matrix, labels):
    """
    Compute inter-cluster distance matrix by averaging distances between all pairs of structures
    in different clusters.
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Full pose-level distance matrix (N x N)
    labels : array-like
        Cluster labels for each structure
    
    Returns:
    --------
    cluster_matrix : np.ndarray
        Cluster-level distance matrix (n_clusters x n_clusters)
    cluster_ids : list
        List of cluster IDs corresponding to matrix rows/cols
    cluster_sizes : dict
        Dictionary mapping cluster_id -> number of members
    """
    unique_clusters = sorted(set(labels))
    n_clusters = len(unique_clusters)
    cluster_matrix = np.zeros((n_clusters, n_clusters))
    cluster_sizes = {}
    
    # Build cluster membership lists
    cluster_members = {}
    for cluster_id in unique_clusters:
        members = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
        cluster_members[cluster_id] = members
        cluster_sizes[cluster_id] = len(members)
    
    # Compute average inter-cluster distances
    for i, cluster_i in enumerate(unique_clusters):
        for j, cluster_j in enumerate(unique_clusters):
            if i == j:
                # Intra-cluster: average pairwise distance within cluster
                members = cluster_members[cluster_i]
                if len(members) > 1:
                    distances = []
                    for m1 in members:
                        for m2 in members:
                            if m1 < m2:
                                distances.append(distance_matrix[m1, m2])
                    cluster_matrix[i, j] = np.mean(distances) if distances else 0.0
                else:
                    cluster_matrix[i, j] = 0.0
            else:
                # Inter-cluster: average distance between all pairs
                members_i = cluster_members[cluster_i]
                members_j = cluster_members[cluster_j]
                distances = []
                for m1 in members_i:
                    for m2 in members_j:
                        distances.append(distance_matrix[m1, m2])
                cluster_matrix[i, j] = np.mean(distances) if distances else 0.0
    
    return cluster_matrix, unique_clusters, cluster_sizes


def compute_intra_cluster_matrices(distance_matrix, labels):
    """
    Compute separate distance matrices for each cluster (intra-cluster distances only).
    Useful for validating cluster cohesiveness.
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Full pose-level distance matrix (N x N)
    labels : array-like
        Cluster labels for each structure
    
    Returns:
    --------
    intra_matrices : dict
        Dictionary mapping cluster_id -> distance matrix for that cluster
    cluster_members : dict
        Dictionary mapping cluster_id -> list of member indices
    """
    unique_clusters = sorted(set(labels))
    intra_matrices = {}
    cluster_members = {}
    
    for cluster_id in unique_clusters:
        members = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
        cluster_members[cluster_id] = members
        
        if len(members) > 1:
            # Extract submatrix for this cluster
            member_array = np.array(members)
            submatrix = distance_matrix[np.ix_(member_array, member_array)]
            intra_matrices[cluster_id] = submatrix
        else:
            # Single member - create 1x1 matrix with zero
            intra_matrices[cluster_id] = np.array([[0.0]])
    
    return intra_matrices, cluster_members


def create_distance_heatmap(distance_matrix, labels, metric_name='Distance'):
    """
    Create interactive heatmap of distance matrix
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Distance matrix
    labels : array-like
        Cluster labels for each structure
    metric_name : str
        Name of the distance metric
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    # Sort by cluster for better visualization
    sorted_indices = np.argsort(labels)
    sorted_matrix = distance_matrix[sorted_indices][:, sorted_indices]
    sorted_labels = np.array(labels)[sorted_indices]
    
    # Create hover text with cluster info
    hover_text = []
    for i in range(len(sorted_labels)):
        row_text = []
        for j in range(len(sorted_labels)):
            row_text.append(
                f"Structure {i} (Cluster {sorted_labels[i]}) vs<br>"
                f"Structure {j} (Cluster {sorted_labels[j]})<br>"
                f"{metric_name}: {sorted_matrix[i,j]:.3f}"
            )
        hover_text.append(row_text)
    
    fig = go.Figure(data=go.Heatmap(
        z=sorted_matrix,
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        colorscale='Viridis',
        colorbar=dict(title=metric_name)
    ))
    
    fig.update_layout(
        title=f'{metric_name} Matrix (sorted by cluster)',
        xaxis_title='Structure Index',
        yaxis_title='Structure Index',
        width=700,
        height=700,
        paper_bgcolor='rgba(0,0,0,0)', 
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(scaleanchor='y', scaleratio=1, constrain='domain'),
        yaxis=dict(scaleanchor='x', constrain='domain'),
        margin=dict(l=60, r=60, t=80, b=60)
    )
    
    return fig


def create_cluster_distance_heatmap(cluster_matrix, cluster_ids, cluster_sizes, metric_name='Distance'):
    """
    Create interactive heatmap of cluster-level distance matrix
    
    Parameters:
    -----------
    cluster_matrix : np.ndarray
        Cluster-level distance matrix (n_clusters x n_clusters)
    cluster_ids : list
        List of cluster IDs
    cluster_sizes : dict
        Dictionary mapping cluster_id -> number of members
    metric_name : str
        Name of the distance metric
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    # Create hover text with cluster info
    hover_text = []
    for i, cluster_i in enumerate(cluster_ids):
        row_text = []
        for j, cluster_j in enumerate(cluster_ids):
            if i == j:
                row_text.append(
                    f"Cluster {cluster_i} (n={cluster_sizes[cluster_i]})<br>"
                    f"Intra-cluster avg {metric_name}: {cluster_matrix[i,j]:.3f}"
                )
            else:
                row_text.append(
                    f"Cluster {cluster_i} (n={cluster_sizes[cluster_i]}) vs<br>"
                    f"Cluster {cluster_j} (n={cluster_sizes[cluster_j]})<br>"
                    f"Inter-cluster avg {metric_name}: {cluster_matrix[i,j]:.3f}"
                )
        hover_text.append(row_text)
    
    # Create custom tick labels with cluster sizes
    tick_labels = [f"C{cid}<br>(n={cluster_sizes[cid]})" for cid in cluster_ids]
    n_clusters = len(cluster_ids)
    
    fig = go.Figure(data=go.Heatmap(
        z=cluster_matrix,
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        colorscale='Viridis',
        colorbar=dict(title=f"Avg {metric_name}"),
        x=tick_labels,
        y=tick_labels
    ))
    
    # Hide tick labels when there are too many clusters (unreadable)
    if n_clusters > 50:
        xaxis_cfg = dict(scaleanchor='y', scaleratio=1, constrain='domain', showticklabels=False)
        yaxis_cfg = dict(scaleanchor='x', constrain='domain', showticklabels=False)
    else:
        xaxis_cfg = dict(scaleanchor='y', scaleratio=1, constrain='domain', tickangle=-45)
        yaxis_cfg = dict(scaleanchor='x', constrain='domain')

    fig.update_layout(
        title=f'Cluster-Level {metric_name} Matrix<br><sub>Diagonal = intra-cluster, Off-diagonal = inter-cluster</sub>',
        xaxis_title='Cluster',
        yaxis_title='Cluster',
        width=700,
        height=700,
        xaxis=xaxis_cfg,
        yaxis=yaxis_cfg,
        margin=dict(l=100, r=60, t=100, b=100)
    )
    
    return fig


def create_intra_cluster_heatmaps(intra_matrices, cluster_members, pdb_files, metric_name='Distance'):
    """
    Create combined figure showing intra-cluster distance matrices for all clusters.
    
    Parameters:
    -----------
    intra_matrices : dict
        Dictionary mapping cluster_id -> distance matrix
    cluster_members : dict
        Dictionary mapping cluster_id -> list of member indices
    pdb_files : list
        List of PDB filenames
    metric_name : str
        Name of the distance metric
    
    Returns:
    --------
    plotly.graph_objects.Figure or list of figures
    """
    from plotly.subplots import make_subplots
    
    cluster_ids = sorted(intra_matrices.keys())
    n_clusters = len(cluster_ids)
    
    # If only a few clusters, create subplots
    if n_clusters <= 6:
        # Create grid layout
        cols = min(3, n_clusters)
        rows = (n_clusters + cols - 1) // cols
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[f"Cluster {cid} (n={len(cluster_members[cid])})" for cid in cluster_ids],
            vertical_spacing=0.12,
            horizontal_spacing=0.08
        )
        
        for idx, cluster_id in enumerate(cluster_ids):
            row = idx // cols + 1
            col = idx % cols + 1
            
            matrix = intra_matrices[cluster_id]
            members = cluster_members[cluster_id]
            
            # Create hover text
            hover_text = []
            for i in range(len(members)):
                row_text = []
                for j in range(len(members)):
                    pose_i = pdb_files[members[i]] if pdb_files else f"Pose {members[i]}"
                    pose_j = pdb_files[members[j]] if pdb_files else f"Pose {members[j]}"
                    row_text.append(
                        f"{pose_i} vs<br>{pose_j}<br>{metric_name}: {matrix[i,j]:.3f}"
                    )
                hover_text.append(row_text)
            
            fig.add_trace(
                go.Heatmap(
                    z=matrix,
                    text=hover_text,
                    hovertemplate='%{text}<extra></extra>',
                    colorscale='Viridis',
                    showscale=(idx == 0),  # Only show colorbar for first plot
                    colorbar=dict(title=metric_name, x=1.02) if idx == 0 else None
                ),
                row=row, col=col
            )
        
        fig.update_layout(
            title=f'Intra-Cluster {metric_name} Matrices - Cohesiveness Validation<br><sub>Lower values = more cohesive cluster</sub>',
            height=300 * rows,
            showlegend=False
        )
        
        return fig
    
    else:
        # Too many clusters - return message to select specific cluster
        fig = go.Figure()
        cluster_stats = "<br>".join([f"Cluster {cid}: {len(cluster_members[cid])} members, "
                                     f"avg intra-distance = {intra_matrices[cid].mean():.3f}"
                                     for cid in cluster_ids[:20]])  # Show first 20
        
        if n_clusters > 20:
            cluster_stats += f"<br>... and {n_clusters - 20} more clusters"
        
        fig.add_annotation(
            text=f"<b>{n_clusters} clusters found - too many to display all matrices</b><br><br>"
                 f"Showing cluster statistics:<br>{cluster_stats}<br><br>"
                 f"Tip: Use cluster-level heatmap to identify cohesiveness",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=12),
            align="left"
        )
        fig.update_layout(height=600)
        
        return fig


def create_distance_heatmap_with_dendrogram(distance_matrix, labels=None, metric_name='Distance', pdb_files=None):
    """
    Create heatmap with hierarchical clustering dendrogram.
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        N×N symmetric distance matrix
    labels : array-like, optional
        Cluster labels for annotation
    metric_name : str
        Name of the distance metric
    pdb_files : list, optional
        List of PDB filenames for labels
        
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform
    import plotly.figure_factory as ff
    import os
    
    n = distance_matrix.shape[0]
    
    # Handle small matrices
    if n < 3:
        # For very small matrices, skip dendrogram and just show heatmap
        hover_labels = [os.path.basename(pdb_files[i]) if pdb_files else f"Structure {i}" for i in range(n)]
        fig = go.Figure(data=go.Heatmap(
            z=distance_matrix,
            x=hover_labels,
            y=hover_labels,
            colorscale='Viridis',
            hovertemplate='%{y} vs %{x}<br>' + metric_name + ': %{z:.3f}<extra></extra>',
            colorbar=dict(title=metric_name)
        ))
        fig.update_layout(
            title=f'{metric_name} (too few samples for dendrogram)',
            height=600,
            template='plotly',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)'
        )
        return fig
    
    # Convert distance matrix to condensed form for scipy
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Perform hierarchical clustering
    Z = linkage(condensed_dist, method='average')
    
    # Create dendrogram to get ordering
    dend = dendrogram(Z, no_plot=True)
    idx = dend['leaves']
    
    # Reorder matrix according to dendrogram
    ordered_matrix = distance_matrix[idx, :][:, idx]
    
    # Create hover labels
    if pdb_files is not None:
        # Check if pdb_files are already strings (e.g., cluster labels) or file paths
        if isinstance(pdb_files[0], str) and not os.path.exists(pdb_files[0]):
            # Already labels, use directly
            hover_labels = [pdb_files[i] for i in idx]
        else:
            # File paths, extract basename
            hover_labels = [os.path.basename(pdb_files[i]) for i in idx]
    else:
        hover_labels = [f"Structure {i}" for i in idx]
    
    # Use plotly figure_factory for cleaner dendrogram
    # Create the clustered heatmap using Plotly's built-in function
    try:
        fig = ff.create_dendrogram(
            distance_matrix,
            distfun=lambda x: condensed_dist,
            linkagefun=lambda x: Z,
            labels=hover_labels if len(hover_labels) <= 50 else None,  # Only show labels if not too many
            orientation='bottom'
        )
        
        # Add the heatmap
        dendro_leaves = fig['layout']['xaxis']['ticktext']
        dendro_side = fig['data'][0]['y']
        
        # Create new figure with heatmap below dendrogram
        from plotly.subplots import make_subplots
        
        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[0.15, 0.85],
            vertical_spacing=0.01,
            subplot_titles=('Dendrogram', f'{metric_name} Heatmap'),
            specs=[[{"type": "scatter"}],
                   [{"type": "heatmap"}]]
        )
        
        # Get dendrogram data
        dendro_fig = ff.create_dendrogram(
            distance_matrix,
            distfun=lambda x: condensed_dist,
            linkagefun=lambda x: Z,
            orientation='bottom'
        )
        
        for trace in dendro_fig['data']:
            fig.add_trace(trace, row=1, col=1)
        
        # Add heatmap
        heatmap = go.Heatmap(
            z=ordered_matrix,
            x=hover_labels,
            y=hover_labels,
            colorscale='Viridis',
            hovertemplate='%{y} vs %{x}<br>' + metric_name + ': %{z:.3f}<extra></extra>',
            colorbar=dict(title=metric_name)
        )
        fig.add_trace(heatmap, row=2, col=1)
        
        # Update layout
        fig.update_layout(
            height=800,
            showlegend=False,
            template='plotly',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            title=f'{metric_name} with Hierarchical Clustering'
        )
        
        # Update axes
        fig.update_xaxes(showticklabels=False, row=1, col=1)
        fig.update_yaxes(showticklabels=False, row=1, col=1)
        
        # Show labels on heatmap only if not too many
        if len(hover_labels) > 50:
            fig.update_xaxes(showticklabels=False, row=2, col=1)
            fig.update_yaxes(showticklabels=False, row=2, col=1)
        
    except Exception as e:
        # Fallback to simple heatmap if dendrogram fails
        print(f"Dendrogram creation failed: {e}, falling back to simple heatmap")
        fig = go.Figure(data=go.Heatmap(
            z=ordered_matrix,
            x=hover_labels,
            y=hover_labels,
            colorscale='Viridis',
            hovertemplate='%{y} vs %{x}<br>' + metric_name + ': %{z:.3f}<extra></extra>',
            colorbar=dict(title=metric_name)
        ))
        fig.update_layout(
            title=f'{metric_name} (dendrogram unavailable)',
            height=600,
            template='plotly',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)'
        )
    
    return fig


def create_cluster_size_distribution(labels, duplicate_groups=None, rescued_cluster_ids=None):
    """
    Create bar chart of cluster size distribution
    
    Parameters:
    -----------
    labels : array-like
        Cluster labels
    duplicate_groups : dict, optional
        Dictionary mapping representative index -> list of duplicate indices
    rescued_cluster_ids : set, optional
        Set of cluster IDs that were rescued from noise
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    import pandas as pd
    
    cluster_counts = pd.Series(labels).value_counts().sort_index()
    
    # Exclude noise cluster (-1) from visualization
    cluster_counts = cluster_counts[cluster_counts.index != -1]
    
    # Identify clusters with duplicates
    clusters_with_duplicates = set()
    if duplicate_groups:
        for rep_idx, duplicates in duplicate_groups.items():
            if len(duplicates) > 1 and rep_idx < len(labels):
                cluster_label = labels[rep_idx]
                if cluster_label != -1:
                    clusters_with_duplicates.add(cluster_label)
    
    # Create labels with asterisks
    labels_text = []
    for idx in cluster_counts.index:
        label = f"Cluster {idx}"
        if rescued_cluster_ids and idx in rescued_cluster_ids:
            label += "**"  # Rescued (all duplicates)
        elif idx in clusters_with_duplicates:
            label += "*"   # Contains some duplicates
        labels_text.append(label)
    
    colors = ['#1f77b4' for _ in cluster_counts.index]
    
    fig = go.Figure(data=[
        go.Bar(
            x=labels_text,
            y=cluster_counts.values,
            marker_color=colors,
            text=cluster_counts.values,
            textposition='auto',
        )
    ])
    
    # Add note about asterisks
    title_note = ""
    if rescued_cluster_ids:
        n_rescued = len([c for c in cluster_counts.index if c in rescued_cluster_ids])
        n_with_dups = len([c for c in cluster_counts.index if c in clusters_with_duplicates and c not in rescued_cluster_ids])
        if n_rescued > 0 or n_with_dups > 0:
            title_note = f"<br><sub>* = has duplicates ({n_with_dups}), ** = all duplicates ({n_rescued})</sub>"
    elif clusters_with_duplicates:
        title_note = f"<br><sub>* = contains near-duplicates ({len(clusters_with_duplicates)} clusters)</sub>"
    
    fig.update_layout(
        title=f'Cluster Size Distribution{title_note}',
        xaxis_title='Cluster',
        yaxis_title='Number of Structures',
        height=400,
        template='plotly', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False,
        xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
        yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0')
    )
    
    return fig


def create_cluster_stability_plot(clusterer, labels, distance_matrix=None, duplicate_groups=None, rescued_cluster_ids=None):
    """
    Create cluster stability plot showing HDBSCAN persistence scores with three-tier fallback
    
    Fallback hierarchy:
    1. HDBSCAN cluster_persistence_ (λ-integral, most accurate)
    2. Hybrid: (λ_range × size) / median_intra_cluster_distance
    3. Pure distance: 1 / median_intra_cluster_distance
    
    Parameters:
    -----------
    clusterer : hdbscan.HDBSCAN
        Fitted HDBSCAN clusterer (must have cluster_persistence_ attribute)
    labels : array-like
        Cluster labels
    distance_matrix : np.ndarray, optional
        Distance matrix for fallback calculations (required for fallback levels 2-3)
    duplicate_groups : dict, optional
        Dictionary mapping representative index -> list of duplicate indices
        Used to mark clusters containing duplicates with asterisks
    rescued_cluster_ids : set, optional
        Set of cluster IDs that were rescued from noise (formed from duplicate groups)
        These clusters are marked with ** to indicate they consist entirely of near-duplicates
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    import pandas as pd
    
    if not hasattr(clusterer, 'cluster_persistence_'):
        # Fallback: return empty plot with message
        fig = go.Figure()
        fig.add_annotation(
            text="Cluster persistence scores not available<br>(requires HDBSCAN clustering)",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Get cluster sizes
    cluster_counts = pd.Series(labels).value_counts().sort_index()
    cluster_ids = [idx for idx in cluster_counts.index if idx != -1]  # Exclude noise
    
    if len(cluster_ids) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No clusters found (all noise)",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Get persistence scores for each cluster
    persistence_scores = clusterer.cluster_persistence_
    cluster_sizes = [cluster_counts[cid] for cid in cluster_ids]
    
    # Initialize metric tracking
    metric_name = "Persistence Score (λ-integral)"
    warning_msg = ""
    fallback_level = 0  # 0 = primary, 1 = hybrid, 2 = distance-only
    
    # Check if persistence scores are valid (not all 1.0 or inf)
    if np.all(persistence_scores == 1.0) or np.any(np.isinf(persistence_scores)):
        # FALLBACK LEVEL 1: Try hybrid approach (λ_range × size) / median_distance
        try:
            if distance_matrix is None:
                raise ValueError("Distance matrix required for hybrid fallback")
            
            tree_df = clusterer.condensed_tree_.to_pandas()
            tree_df_finite = tree_df[np.isfinite(tree_df['lambda_val'])].copy()
            
            if len(tree_df_finite) == 0:
                raise ValueError("No finite lambda values in condensed tree")
            
            # Compute median intra-cluster distances and lambda-based scores
            median_distances = []
            lambda_range_scores = []
            
            print(f"[Stability Fallback] Computing hybrid metric for {len(cluster_ids)} clusters")
            print(f"[Stability Fallback] Tree has {len(tree_df_finite)} nodes with finite lambda")
            
            for cid in cluster_ids:
                # Get cluster members (actual point indices)
                cluster_mask = np.array(labels) == cid
                cluster_indices = np.where(cluster_mask)[0]
                
                # Compute median intra-cluster distance (optimized for large clusters)
                if len(cluster_indices) > 1:
                    # Use numpy indexing for fast extraction of submatrix
                    cluster_submatrix = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
                    # Get upper triangle (excluding diagonal) for pairwise distances
                    upper_tri_indices = np.triu_indices_from(cluster_submatrix, k=1)
                    cluster_dists = cluster_submatrix[upper_tri_indices]
                    
                    # Compute median and add small noise to prevent exact zeros
                    median_dist = np.median(cluster_dists)
                    # Add tiny epsilon to prevent division by zero (1e-6 is ~0.000001 Angstrom)
                    median_dist = max(median_dist, 1e-6)
                else:
                    median_dist = 1.0  # Single member, use neutral value
                
                median_distances.append(median_dist)
                
                # Get lambda range by looking at all splits involving this cluster's points
                # The condensed tree records when each point leaves the cluster hierarchy
                # We look for rows where the child is one of our cluster's points
                cluster_point_rows = tree_df_finite[tree_df_finite['child'].isin(cluster_indices)]
                
                if len(cluster_point_rows) > 0:
                    lambda_vals = cluster_point_rows['lambda_val'].values
                    lambda_vals_finite = lambda_vals[np.isfinite(lambda_vals)]
                    
                    if len(lambda_vals_finite) > 1:
                        # Lambda range: max (when first point splits) - min (when last point splits)
                        # Larger range = cluster persists across more density levels = more stable
                        lambda_range = lambda_vals_finite.max() - lambda_vals_finite.min()
                        # Add small epsilon to lambda_range to prevent exact zeros from identical structures
                        lambda_range = max(lambda_range, 1e-6)
                        # Multiply by cluster size for stability score
                        lambda_range_scores.append(lambda_range * len(cluster_indices))
                    else:
                        # All points split at same lambda - no persistence
                        # Fall back to cluster size
                        lambda_range_scores.append(float(len(cluster_indices)))
                else:
                    # No tree nodes found - use cluster size as proxy
                    lambda_range_scores.append(float(len(cluster_indices)))
                
                print(f"  Cluster {cid}: median_dist={median_dist:.3f}, lambda_range_score={lambda_range_scores[-1]:.3f}, n_points_in_tree={len(cluster_point_rows)}")
            
            # Check if we got valid lambda scores
            if all(lrs == 0.0 for lrs in lambda_range_scores):
                print("[Stability Fallback] All lambda scores are zero, falling back to distance-only")
                raise ValueError("Lambda range scores all zero")
            
            # Hybrid stability: (λ_range × size) / median_distance
            # Higher lambda range = more stable across densities
            # Lower median distance = more cohesive
            persistence_scores = np.array([
                (lrs / md) if md > 0 and lrs > 0 else 0.0
                for lrs, md in zip(lambda_range_scores, median_distances)
            ])
            
            print(f"[Stability Fallback] Hybrid scores: {persistence_scores}")
            
            metric_name = "Hybrid Stability"
            warning_msg = (
                "<br><sub>ℹ️ <b>Fallback Level 1 (Hybrid):</b> HDBSCAN persistence invalid (infinite λ values).<br>"
                "Using: (λ_range × avg_size) / median_intra_cluster_distance<br>"
                "Higher = more stable across densities + more cohesive</sub>"
            )
            fallback_level = 1
                
        except Exception as e:
            # FALLBACK LEVEL 2: Pure distance-based (1 / median_intra_cluster_distance)
            if distance_matrix is not None:
                median_distances = []
                
                print(f"[Stability Fallback] Lambda extraction failed: {str(e)[:80]}")
                print(f"[Stability Fallback] Using distance-only metric for {len(cluster_ids)} clusters")
                
                for cid in cluster_ids:
                    # Get cluster members
                    cluster_mask = np.array(labels) == cid
                    cluster_indices = np.where(cluster_mask)[0]
                    
                    # Compute median intra-cluster distance (optimized)
                    if len(cluster_indices) > 1:
                        # Use numpy indexing for fast extraction
                        cluster_submatrix = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
                        upper_tri_indices = np.triu_indices_from(cluster_submatrix, k=1)
                        cluster_dists = cluster_submatrix[upper_tri_indices]
                        
                        median_dist = np.median(cluster_dists)
                        # Add epsilon to prevent division by zero
                        median_dist = max(median_dist, 1e-6)
                    else:
                        median_dist = 1.0  # Single member, use neutral value
                    
                    median_distances.append(median_dist)
                
                # Stability = 1 / median_distance (lower distance = higher stability)
                persistence_scores = np.array([
                    1.0 / md if md > 0 else 0.0
                    for md in median_distances
                ])
                
                print(f"[Stability Fallback] Distance-based scores: {persistence_scores}")
                
                metric_name = "Distance-Based Stability"
                warning_msg = (
                    "<br><sub>⚠️ <b>Fallback Level 2 (Distance-Only):</b> HDBSCAN persistence and λ-values unavailable.<br>"
                    f"Using: 1 / median_intra_cluster_distance (Reason: {str(e)[:80]})<br>"
                    "Higher = more cohesive cluster (lower internal distances)</sub>"
                )
                fallback_level = 2
            else:
                # No distance matrix - use cluster sizes as last resort
                persistence_scores = np.array(cluster_sizes, dtype=float)
                metric_name = "Cluster Size"
                warning_msg = (
                    f"<br><sub>⚠️ <b>Fallback Level 3 (Size-Only):</b> No distance matrix available.<br>"
                    f"Using cluster size as stability proxy (Reason: {str(e)[:80]})</sub>"
                )
                fallback_level = 3
    else:
        # Primary metric - HDBSCAN persistence is valid
        metric_name = "Persistence Score (λ-integral)"
        warning_msg = "<br><sub>✓ <b>Primary Metric:</b> HDBSCAN cluster_persistence_ (λ-integral)</sub>"
        fallback_level = 0
    
    # Sort by persistence (stability)
    sorted_indices = np.argsort(persistence_scores)[::-1]
    sorted_clusters = [cluster_ids[i] for i in sorted_indices]
    sorted_persistence = [persistence_scores[i] for i in sorted_indices]
    sorted_sizes = [cluster_sizes[i] for i in sorted_indices]
    
    # Check which clusters contain duplicates
    clusters_with_duplicates = set()
    if duplicate_groups:
        # Find which cluster each representative belongs to
        for rep_idx, duplicates in duplicate_groups.items():
            if len(duplicates) > 1:  # Has duplicates
                # Find the cluster label for this representative
                if rep_idx < len(labels):
                    cluster_label = labels[rep_idx]
                    if cluster_label != -1:
                        clusters_with_duplicates.add(cluster_label)
    
    # Use log scale for better visualization of wide-ranging values
    # Add 1 to handle zero values: log10(x + 1)
    log_persistence = [np.log10(p + 1) for p in sorted_persistence]
    
    # Compute visual cap at Q3 + 3×IQR on log scale (gentler cap for log)
    # Using 3×IQR instead of 1.5×IQR since log scale already compresses outliers
    q1_log = np.percentile(log_persistence, 25)
    q3_log = np.percentile(log_persistence, 75)
    iqr_log = q3_log - q1_log
    visual_cap_log = q3_log + 3 * iqr_log
    
    # Identify extreme outliers that exceed even the log cap
    outlier_mask = np.array(log_persistence) > visual_cap_log
    n_outliers = outlier_mask.sum()
    
    # Create display values (capped on log scale) but keep real values for hover
    display_log_persistence = [min(lp, visual_cap_log) for lp in log_persistence]
    
    # Color by log(size)
    log_sizes = np.log10(sorted_sizes)
    
    # Create x-axis labels with asterisks for clusters containing duplicates
    x_labels = []
    for cid in sorted_clusters:
        label = f"Cluster {cid}"
        if rescued_cluster_ids and cid in rescued_cluster_ids:
            label += "**"  # Double asterisk for rescued (all duplicates)
        elif cid in clusters_with_duplicates:
            label += "*"   # Single asterisk for contains some duplicates
        x_labels.append(label)
    
    fig = go.Figure(data=[
        go.Bar(
            x=x_labels,
            y=display_log_persistence,  # Use log scale with cap
            marker=dict(
                color=log_sizes,
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(
                    title="log₁₀(Size)",
                )
            ),
            text=[
                f"{p:.1e}" if p > 1000 else f"{p:.1f}"  # Scientific notation for large values
                for p in sorted_persistence
            ],
            textposition='outside',
            hovertemplate=(
                '%{x}<br>' + 
                metric_name + ': %{customdata[0]:.3e}<br>' +  # Scientific notation in hover
                'log₁₀(' + metric_name + ' + 1): %{y:.2f}<br>' +  # Show log value too
                'Size: %{customdata[1]}<br>' +
                '<extra></extra>'
            ),
            customdata=[[true_p, size] for true_p, size in zip(sorted_persistence, sorted_sizes)]
        )
    ])
    
    # Add warning message as annotation below plot to prevent overlap with bars
    # Update warning to include log scale info and outlier info if present
    log_scale_info = "📊 Using log₁₀ scale to handle wide range of values (1 to 7M+)"
    
    if n_outliers > 0:
        # Convert log cap back to original scale for display
        original_cap = 10**visual_cap_log - 1
        outlier_info = f"\n⚠️ {n_outliers} extreme outlier(s) capped at log={visual_cap_log:.1f} (≈{original_cap:.1e}) - hover to see true values"
        final_warning = f"{log_scale_info}\n{outlier_info}"
        if warning_msg:
            final_warning = warning_msg + f"\n{final_warning}"
    else:
        final_warning = f"{log_scale_info}"
        if warning_msg:
            final_warning = warning_msg + f"\n{final_warning}"
    
    # Add note about asterisks if we have duplicates
    duplicate_note = ""
    if rescued_cluster_ids:
        n_rescued = len([c for c in sorted_clusters if c in rescued_cluster_ids])
        n_with_dups = len([c for c in sorted_clusters if c in clusters_with_duplicates and c not in rescued_cluster_ids])
        duplicate_note = f" | * = has duplicates ({n_with_dups}), ** = all duplicates ({n_rescued})"
    elif clusters_with_duplicates:
        duplicate_note = f" | * = contains near-duplicates ({len(clusters_with_duplicates)} clusters)"
    
    fig.update_layout(
        title=f'Cluster Stability ({metric_name}){duplicate_note}',
        xaxis_title='Cluster (sorted by stability)',
        yaxis_title=f'log₁₀({metric_name} + 1)' + (' [capped]' if n_outliers > 0 else ''),
        height=450,  # Standard height without annotation space
        template='plotly',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False,
        xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
        yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
        margin=dict(b=80)  # Normal bottom margin
    )
    
    # Store warning message in figure metadata for UI to display
    # Extract clean text without HTML tags
    import re
    clean_msg = re.sub('<[^<]+?>', '', final_warning) if final_warning else ""
    clean_msg = clean_msg.replace('&times;', '×')
    fig.update_layout(meta={'warning_text': clean_msg})
    
    return fig


def create_cluster_diversity_plot(labels, rmsd_matrix, pdb_files, sasa_data=None):
    """
    Create pose diversity stripplot showing RMSD from medoid for each cluster
    
    Parameters:
    -----------
    labels : array-like
        Cluster labels
    rmsd_matrix : np.ndarray
        RMSD distance matrix
    pdb_files : list
        List of PDB filenames
    sasa_data : pd.DataFrame, optional
        DataFrame with ΔSASA data (must have 'structure' and 'total_delta_sasa' columns)
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    import pandas as pd
    import os
    
    # Get unique clusters (excluding noise)
    unique_clusters = sorted([c for c in set(labels) if c != -1])
    
    if len(unique_clusters) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No clusters found (all noise)",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Sort clusters by size (descending)
    cluster_sizes = {c: sum(1 for l in labels if l == c) for c in unique_clusters}
    sorted_clusters = sorted(unique_clusters, key=lambda c: cluster_sizes[c], reverse=True)
    
    # Prepare data
    plot_data = []
    
    for cluster_id in sorted_clusters:
        # Get indices for this cluster
        cluster_indices = [i for i, l in enumerate(labels) if l == cluster_id]
        
        if len(cluster_indices) < 2:
            continue
        
        # Find medoid (structure with minimum average RMSD to others in cluster)
        cluster_rmsd = rmsd_matrix[np.ix_(cluster_indices, cluster_indices)]
        medoid_idx_local = np.argmin(cluster_rmsd.sum(axis=0))
        medoid_idx_global = cluster_indices[medoid_idx_local]
        
        # Compute RMSD from medoid for each structure
        for i, global_idx in enumerate(cluster_indices):
            rmsd_from_medoid = rmsd_matrix[global_idx, medoid_idx_global]
            
            # Get ΔSASA if available
            color_value = None
            if sasa_data is not None:
                structure_name = os.path.basename(pdb_files[global_idx])
                sasa_row = sasa_data[sasa_data['structure'] == structure_name]
                if len(sasa_row) > 0 and 'total_delta_sasa' in sasa_row.columns:
                    color_value = sasa_row['total_delta_sasa'].values[0]
            
            plot_data.append({
                'cluster': cluster_id,
                'cluster_label': f"Cluster {cluster_id} (n={cluster_sizes[cluster_id]})",
                'rmsd': rmsd_from_medoid,
                'structure': os.path.basename(pdb_files[global_idx]),
                'delta_sasa': color_value,
                'is_medoid': (global_idx == medoid_idx_global)
            })
    
    df = pd.DataFrame(plot_data)
    
    if len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No data to plot",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Compute global SASA range for consistent colorbar across all clusters
    use_sasa_coloring = sasa_data is not None and df['delta_sasa'].notna().any()
    if use_sasa_coloring:
        sasa_min = df['delta_sasa'].min()
        sasa_max = df['delta_sasa'].max()
        sasa_mean = df['delta_sasa'].mean()
    
    # Create figure
    fig = go.Figure()
    
    # Add traces for each cluster
    first_sasa_trace = True  # Track first trace with SASA coloring for colorbar
    
    for cluster_id in sorted_clusters:
        cluster_df = df[df['cluster'] == cluster_id]
        
        if len(cluster_df) == 0:
            continue
        
        # Separate medoid and non-medoid
        medoid_df = cluster_df[cluster_df['is_medoid']]
        non_medoid_df = cluster_df[~cluster_df['is_medoid']]
        
        # Non-medoid points
        if len(non_medoid_df) > 0:
            if use_sasa_coloring:
                # Color by ΔSASA - use global range and show colorbar only once
                fig.add_trace(go.Scatter(
                    x=non_medoid_df['rmsd'],
                    y=non_medoid_df['cluster_label'],
                    mode='markers',
                    marker=dict(
                        color=non_medoid_df['delta_sasa'],
                        colorscale='RdBu_r',
                        cmin=sasa_min,
                        cmax=sasa_max,
                        showscale=first_sasa_trace,  # Only show colorbar for first trace
                        colorbar=dict(
                            title="ΔSASA (Å²)",
                            tickmode="array",
                            tickvals=[sasa_min, sasa_mean, sasa_max],
                            ticktext=[f"{sasa_min:.0f}", f"{sasa_mean:.0f}", f"{sasa_max:.0f}"]
                        ) if first_sasa_trace else None,
                        size=8
                    ),
                    text=non_medoid_df['structure'],
                    hovertemplate='%{text}<br>RMSD: %{x:.2f} Å<br>ΔSASA: %{marker.color:.1f} Å²<extra></extra>',
                    name=f'Cluster {cluster_id}',
                    showlegend=False
                ))
                first_sasa_trace = False  # Colorbar shown, don't show again
            else:
                # Uniform color
                fig.add_trace(go.Scatter(
                    x=non_medoid_df['rmsd'],
                    y=non_medoid_df['cluster_label'],
                    mode='markers',
                    marker=dict(
                        size=8,
                        color='steelblue',
                        line=dict(width=0.5, color='white')
                    ),
                    text=non_medoid_df['structure'],
                    hovertemplate='%{text}<br>RMSD: %{x:.2f} Å<extra></extra>',
                    name=f'Cluster {cluster_id}',
                    showlegend=False
                ))
        
        # Medoid point (highlighted)
        if len(medoid_df) > 0:
            fig.add_trace(go.Scatter(
                x=medoid_df['rmsd'],
                y=medoid_df['cluster_label'],
                mode='markers',
                marker=dict(
                    size=12,
                    color='gold',
                    symbol='star',
                    line=dict(width=1, color='black')
                ),
                text=medoid_df['structure'],
                hovertemplate='<b>MEDOID</b><br>%{text}<br>RMSD: %{x:.2f} Å<extra></extra>',
                name='Medoid',
                showlegend=False  # Show legend only once
            ))
    
    fig.update_layout(
        title='Pose Diversity within Clusters (RMSD from Medoid)',
        xaxis_title='RMSD from Cluster Medoid (Å)',
        yaxis_title='Cluster (sorted by size)',
        height=500,  # Fixed height similar to motif match chart, users can scroll if needed
        template='plotly',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
        yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0', categoryorder='array', categoryarray=[f"Cluster {cid} (n={cluster_sizes[cid]})" for cid in sorted_clusters]),
        hovermode='closest',
        margin=dict(l=120, r=20, t=60, b=60)  # Tighter margins
    )
    
    return fig


def create_single_intra_cluster_heatmap(distance_matrix, labels, cluster_id, pdb_files=None, metric_name='Distance'):
    """
    Create heatmap for a single cluster's intra-cluster distances.
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Full distance matrix
    labels : array-like
        Cluster labels
    cluster_id : int
        Which cluster to display
    pdb_files : list, optional
        PDB filenames for labels
    metric_name : str
        Name of the distance metric
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    # Get members of this cluster
    members = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
    
    if len(members) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Cluster {cluster_id} has no members",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    if len(members) == 1:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Cluster {cluster_id} has only 1 member<br>Cannot compute intra-cluster distances",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Extract submatrix
    member_array = np.array(members)
    submatrix = distance_matrix[np.ix_(member_array, member_array)]
    
    # Create hover text with structure names
    hover_text = []
    for i in range(len(members)):
        row_text = []
        for j in range(len(members)):
            if pdb_files:
                struct_i = pdb_files[members[i]]
                struct_j = pdb_files[members[j]]
            else:
                struct_i = f"Pose {members[i]}"
                struct_j = f"Pose {members[j]}"
            
            row_text.append(
                f"{struct_i} vs<br>{struct_j}<br>{metric_name}: {submatrix[i,j]:.3f}"
            )
        hover_text.append(row_text)
    
    # Create labels for axes (use just filename, not full path)
    if pdb_files:
        import os
        tick_labels = [os.path.basename(pdb_files[idx]) for idx in members]
    else:
        tick_labels = [f"Pose {idx}" for idx in members]
    
    fig = go.Figure(data=go.Heatmap(
        z=submatrix,
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        colorscale='Viridis',
        colorbar=dict(title=metric_name),
        x=tick_labels,
        y=tick_labels
    ))
    
    avg_distance = submatrix[np.triu_indices_from(submatrix, k=1)].mean()
    
    fig.update_layout(
        title=f'Cluster {cluster_id} Intra-Cluster {metric_name}<br>'
              f'<sub>{len(members)} members | Avg {metric_name}: {avg_distance:.3f}</sub>',
        xaxis_title='Structure',
        yaxis_title='Structure',
        width=700,
        height=700,
        xaxis=dict(scaleanchor='y', scaleratio=1, constrain='domain', tickangle=-45),
        yaxis=dict(scaleanchor='x', constrain='domain'),
        margin=dict(l=100, r=60, t=100, b=150)
    )
    
    return fig


def create_motif_match_score_plot(labels, motif_match_scores, pdb_files):
    """
    Create bar chart showing motif match scores per structure and cluster averages
    
    Parameters:
    -----------
    labels : array-like
        Cluster labels for each structure
    motif_match_scores : list
        Motif match scores (% of motif residues making contact) for each structure
    pdb_files : list
        List of PDB filenames
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    import pandas as pd
    
    # Check if we have valid motif scores
    if motif_match_scores is None or len(motif_match_scores) == 0 or all(s is None for s in motif_match_scores):
        fig = go.Figure()
        fig.add_annotation(
            text="No motif match scores available<br>(motif residues not specified)",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Get unique clusters (excluding noise)
    unique_clusters = sorted([c for c in set(labels) if c != -1])
    
    if len(unique_clusters) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No clusters found (all noise)",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Compute cluster average match scores
    cluster_avg_scores = {}
    cluster_sizes = {}
    
    for cluster_id in unique_clusters:
        cluster_indices = [i for i, l in enumerate(labels) if l == cluster_id]
        cluster_scores = [motif_match_scores[i] for i in cluster_indices if motif_match_scores[i] is not None]
        
        if len(cluster_scores) > 0:
            cluster_avg_scores[cluster_id] = np.mean(cluster_scores)
            cluster_sizes[cluster_id] = len(cluster_indices)
        else:
            cluster_avg_scores[cluster_id] = 0.0
            cluster_sizes[cluster_id] = len(cluster_indices)
    
    # Sort clusters by average match score (descending)
    sorted_clusters = sorted(unique_clusters, key=lambda c: cluster_avg_scores[c], reverse=True)
    
    # Create bar chart
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=[f"Cluster {cid}" for cid in sorted_clusters],
        y=[cluster_avg_scores[cid] for cid in sorted_clusters],
        marker=dict(
            color=[cluster_avg_scores[cid] for cid in sorted_clusters],
            colorscale='RdYlGn',  # Red (low) to Yellow to Green (high)
            cmin=0,
            cmax=100,
            showscale=True,
            colorbar=dict(title="Match %")
        ),
        text=[f"{cluster_avg_scores[cid]:.1f}%" for cid in sorted_clusters],
        textposition='outside',
        hovertemplate=(
            '%{x}<br>' +
            'Avg Match Score: %{y:.1f}%<br>' +
            'Cluster Size: %{customdata}<br>' +
            '<extra></extra>'
        ),
        customdata=[cluster_sizes[cid] for cid in sorted_clusters]
    ))
    
    fig.update_layout(
        title='Motif Match Score by Cluster<br><sub>% of motif residues making contact (higher = better match)</sub>',
        xaxis_title='Cluster (sorted by match score)',
        yaxis_title='Average Motif Match Score (%)',
        height=450,
        template='plotly',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False,
        xaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0'),
        yaxis=dict(zerolinecolor='#808080', gridcolor='#e0e0e0', range=[0, 105]),
        margin=dict(l=80, r=20, t=80, b=80)
    )
    
    return fig


def create_contact_residue_heatmap(contact_residue_pairs, labels, pdb_files, cluster_id=None):
    """
    Create intra-cluster contact heatmap showing which residues each structure contacts.
    X-axis: protein residues, Y-axis: structures in cluster
    Color: actual contact distances in Ångströms (red = closer, blue = farther)
    
    Parameters:
    -----------
    contact_residue_pairs : list of sets
        Contact residue pairs for each structure (from clusterer)
        Each element is a set of (protein_res, nucleic_res, distance) tuples or dicts
    labels : array-like
        Cluster labels for each structure
    pdb_files : list
        List of PDB filenames
    cluster_id : int, optional
        Which cluster to display. If None, shows first cluster.
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    import pandas as pd
    import numpy as np
    
    # Get unique clusters (exclude noise)
    unique_clusters = sorted([c for c in set(labels) if c != -1])
    
    if len(unique_clusters) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No clusters found",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Select cluster to display
    if cluster_id is None or cluster_id not in unique_clusters:
        cluster_id = unique_clusters[0]
    
    # Get structures in this cluster
    cluster_mask = np.array(labels) == cluster_id
    cluster_indices = np.where(cluster_mask)[0]
    
    if len(cluster_indices) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Cluster {cluster_id} has no structures",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        fig.update_layout(height=400)
        return fig
    
    # Extract all unique protein residues contacted by this cluster
    cluster_protein_residues = set()
    for idx in cluster_indices:
        contact_pairs = contact_residue_pairs[idx]
        for pair in contact_pairs:
            # pair can be: dict, tuple, or set element
            if isinstance(pair, dict):
                # Dictionary format: {'protein_residue': 'A:ARG10', 'nucleic_residue': ...}
                protein_res = pair.get('protein_residue', pair.get('protein_res', ''))
            elif isinstance(pair, tuple):
                # Tuple format: ('A:ARG10', 'B:DG5') or ('A:ARG10', 'B:DG5', distance)
                protein_res = pair[0]
            else:
                # Unknown format, skip
                continue
            
            if protein_res:
                cluster_protein_residues.add(protein_res)
    
    # Sort residues by chain and number
    def parse_residue(res_str):
        # Format: "CHAIN:RESNAME123"
        if ':' in res_str:
            chain, res_part = res_str.split(':', 1)
            res_num = int(''.join(filter(str.isdigit, res_part)))
            return (chain, res_num, res_str)
        return ('', 0, res_str)
    
    sorted_residues = sorted(cluster_protein_residues, key=parse_residue)
    
    # Limit to most frequent residues if too many
    MAX_RESIDUES = 80
    if len(sorted_residues) > MAX_RESIDUES:
        # Count frequency within cluster
        residue_counts = {}
        for res in sorted_residues:
            count = 0
            for idx in cluster_indices:
                for pair in contact_residue_pairs[idx]:
                    # Extract protein residue based on format
                    if isinstance(pair, dict):
                        protein_res = pair.get('protein_residue', pair.get('protein_res', ''))
                    elif isinstance(pair, tuple):
                        protein_res = pair[0]
                    else:
                        continue
                    
                    if protein_res == res:
                        count += 1
                        break  # Count each structure only once per residue
            residue_counts[res] = count
        
        # Keep top MAX_RESIDUES most frequent
        top_residues = sorted(residue_counts.keys(), key=lambda r: residue_counts[r], reverse=True)[:MAX_RESIDUES]
        sorted_residues = sorted(top_residues, key=parse_residue)  # Re-sort by position
    
    # Build contact distance matrix: structures x residues
    # Binary mode: 1.0 = contact, NaN = no contact
    # Distance mode: actual minimum distance in Å
    contact_matrix = np.full((len(cluster_indices), len(sorted_residues)), np.nan)
    
    for struct_idx_in_cluster, global_idx in enumerate(cluster_indices):
        contact_pairs = contact_residue_pairs[global_idx]
        
        for res_idx, residue in enumerate(sorted_residues):
            # Find minimum distance to this residue
            min_dist = None
            for pair in contact_pairs:
                # Extract protein residue identifier and distance based on format
                if isinstance(pair, dict):
                    protein_res = pair.get('protein_residue', pair.get('protein_res', ''))
                    # Extract distance (prefer min_residue_distance, fallback to distance)
                    if 'min_residue_distance' in pair:
                        pair_distance = pair['min_residue_distance']
                    elif 'distance' in pair:
                        pair_distance = pair['distance']
                    else:
                        pair_distance = 1.0  # Fallback for old data
                elif isinstance(pair, tuple):
                    protein_res = pair[0]
                    # Check if distance is available (3-element tuple)
                    if len(pair) >= 3:
                        pair_distance = pair[2]
                    else:
                        pair_distance = 1.0  # Fallback
                else:
                    continue
                
                if protein_res == residue:
                    # Use minimum distance
                    if min_dist is None or pair_distance < min_dist:
                        min_dist = pair_distance
            
            if min_dist is not None:
                contact_matrix[struct_idx_in_cluster, res_idx] = min_dist
    
    # Cluster structures by contact pattern similarity using hierarchical clustering
    # This orders structures with similar contact patterns together
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist
    
    if len(cluster_indices) > 2:
        # Compute distance between structures based on contact patterns
        # Replace NaN with large value for distance computation (only for clustering)
        contact_matrix_for_clustering = np.nan_to_num(contact_matrix.copy(), nan=10.0)
        
        # Check if we have variation in the data (not all identical rows)
        row_stds = np.std(contact_matrix_for_clustering, axis=1)
        has_variation = np.any(row_stds > 0)
        
        # Compute pairwise distances between structures (rows)
        if has_variation:
            try:
                row_distances = pdist(contact_matrix_for_clustering, metric='euclidean')
                # Check for valid distances (not all zero)
                if np.max(row_distances) > 1e-10:
                    row_linkage = linkage(row_distances, method='average')
                    row_dendrogram = dendrogram(row_linkage, no_plot=True)
                    row_order = row_dendrogram['leaves']
                    
                    # Reorder matrix and labels according to clustering (keep original NaN values)
                    contact_matrix = contact_matrix[row_order, :]
                    cluster_indices_ordered = [cluster_indices[i] for i in row_order]
                else:
                    # All distances are zero (identical structures), keep original order
                    cluster_indices_ordered = cluster_indices
            except Exception as e:
                # If clustering fails, keep original order
                print(f"[Contact Heatmap] Clustering failed: {e}, using original order")
                cluster_indices_ordered = cluster_indices
        else:
            # All rows identical, no need to cluster
            cluster_indices_ordered = cluster_indices
    else:
        # Too few structures to cluster
        cluster_indices_ordered = cluster_indices
    
    # Create structure labels (now in clustered order)
    structure_labels = [os.path.basename(pdb_files[idx]) for idx in cluster_indices_ordered]
    
    # Create hover text
    hover_text = []
    for struct_idx_in_cluster, global_idx in enumerate(cluster_indices_ordered):
        row_text = []
        for res_idx, residue in enumerate(sorted_residues):
            value = contact_matrix[struct_idx_in_cluster, res_idx]
            struct_name = structure_labels[struct_idx_in_cluster]
            if np.isnan(value):
                row_text.append(f"{struct_name}<br>{residue}<br>No contact")
            else:
                row_text.append(f"{struct_name}<br>{residue}<br>Distance: {value:.2f} Å")
        hover_text.append(row_text)
    
    # Create figure
    fig = go.Figure(data=go.Heatmap(
        z=contact_matrix,
        x=sorted_residues,
        y=structure_labels,
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        colorscale='RdYlBu_r',  # Red (close) -> Yellow -> Blue (far)
        reversescale=False,
        colorbar=dict(title="Contact<br>Distance (Å)"),
        zmin=0,
        zmax=6  # Typical contact cutoff
    ))
    
    # Update layout with tight margins
    # Calculate optimal height: base height + rows * pixels per row
    n_structures = len(cluster_indices)
    row_height = 25 if n_structures < 20 else 18 if n_structures < 40 else 15
    plot_height = max(300, n_structures * row_height)
    
    fig.update_layout(
        title=f'Cluster {cluster_id} Contact Distance Pattern<br><sub>Structures ordered by contact similarity (hierarchical clustering). Redder = closer contact.</sub>',
        xaxis_title='Protein Residue',
        yaxis_title='Structure',
        height=plot_height,
        template='plotly',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(
            tickangle=-45,
            tickfont=dict(size=9 if len(sorted_residues) > 50 else 10),
            side='bottom'
        ),
        yaxis=dict(
            tickfont=dict(size=8 if n_structures > 30 else 9),
            automargin=True
        ),
        margin=dict(l=100, r=80, t=80, b=100),
        autosize=True
    )
    
    return fig


def create_hdbscan_condensed_tree(clusterer, labels=None, labels_for_tree=None, pdb_files=None):
    """
    Create HDBSCAN condensed tree plot showing cluster hierarchy.
    
    Uses tree-only visualization (no ellipses) which works reliably.
    
    Parameters:
    -----------
    clusterer : hdbscan.HDBSCAN or PDBClusterer
        Fitted HDBSCAN clusterer object with condensed_tree_ attribute
    labels : np.ndarray, optional
        Cluster labels for statistics display
    labels_for_tree, pdb_files : optional
        Unused (kept for backwards compatibility)
    
    Returns:
    --------
    str : HTML string with base64-encoded PNG image
    """
    # Get the actual HDBSCAN clusterer
    actual_clusterer = clusterer
    if hasattr(clusterer, 'hdbscan_clusterer'):
        actual_clusterer = clusterer.hdbscan_clusterer
        print(f"[Tree Viz] Extracted hdbscan_clusterer from PDBClusterer", flush=True)
        print(f"[Tree Viz] Clusterer type: {type(actual_clusterer)}", flush=True)
        print(f"[Tree Viz] Labels shape: {len(actual_clusterer.labels_)}", flush=True)
    else:
        print(f"[Tree Viz] Using clusterer directly (no hdbscan_clusterer attr)", flush=True)
    
    if not hasattr(actual_clusterer, 'condensed_tree_'):
        return '<div style="padding: 40px; text-align: center; color: #999;"><p>HDBSCAN condensed tree not available.<br>Re-run clustering to generate.</p></div>'
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        #from matplotlib import axes
        from io import BytesIO
        import base64
        import seaborn as sns
        
        print("[Tree Viz] Rendering tree-only visualization (ellipses disabled due to matplotlib bug)...", flush=True)
        
        fig_mpl, ax = plt.subplots(1, 1, figsize=(14, 8))
        actual_clusterer.condensed_tree_.plot(
            select_clusters=True,
            selection_palette=sns.color_palette('deep'),
            axis=ax,
            log_size=True,  # Use log scale for cluster size
            cmap='viridis'
        )
        
        # Add title with statistics
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0) if labels is not None else "N/A"
        n_total = len(labels) if labels is not None else "N/A"
        n_noise = sum(1 for l in labels if l == -1) if labels is not None else "N/A"
        
        plt.title(
            f'HDBSCAN Condensed Tree | {n_clusters} Clusters, {n_total} Poses, {n_noise} Noise\n'
            f'Tree structure shows cluster hierarchy (color gradient indicates density)',
            fontsize=14, pad=20
        )
        ax.set_xlabel('Cluster Size (log scale)', fontsize=12)
        ax.set_ylabel('Lambda (λ = 1/distance)', fontsize=12)
        
        # Convert to base64 image
        buf = BytesIO()
        plt.tight_layout()
        plt.savefig(buf, dpi=300, bbox_inches='tight', facecolor='white')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode()
        plt.close()
        
        print("[Tree Viz] ✅ Tree rendered successfully", flush=True)
        
        return f'''
        <div style="text-align: center; background: white; padding: 10px;">
            <img src="data:image/png;base64,{img_base64}" 
                 style="max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px;"
                 alt="HDBSCAN Condensed Tree" />
        </div>'''
        
    except Exception as e:
        import traceback
        print(f"[Tree Viz] ❌ Failed: {e}", flush=True)
        print(f"[Tree Viz] Traceback:\n{traceback.format_exc()}", flush=True)
        return f'''<div style="padding: 40px; text-align: center; color: #ff9800;">
            <p><strong>⚠️ Tree visualization failed</strong></p>
            <pre style="text-align: left; font-size: 11px; color: #666;">{e}</pre>
        </div>'''


def create_hdbscan_tree(clusterer, labels=None, pdb_files=None):
    """
    Create HDBSCAN single linkage tree plot (matplotlib version).
    
    The single linkage tree shows the minimum spanning tree used by HDBSCAN.
    This is the HDBSCAN-recommended visualization for understanding the hierarchical structure.
    Unlike traditional dendrograms, this scales well to large datasets (1000+ points).
    
    Parameters:
    -----------
    clusterer : hdbscan.HDBSCAN
        Fitted HDBSCAN clusterer object (must have single_linkage_tree_ attribute)
    labels : np.ndarray, optional
        Cluster labels for statistics
    pdb_files : list, optional
        PDB filenames for statistics
    
    Returns:
    --------
    str : HTML string with base64-encoded PNG image
    """
    if not hasattr(clusterer, 'single_linkage_tree_'):
        return '<div style="padding: 40px; text-align: center; color: #999;"><p>HDBSCAN linkage tree not available.<br>Re-run clustering to generate.</p></div>'
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from io import BytesIO
        import base64
        
        # Create matplotlib figure using HDBSCAN's built-in plot
        fig_mpl, ax = plt.subplots(1, 1, figsize=(16, 8))
        clusterer.single_linkage_tree_.plot(
            cmap='viridis',
            colorbar=True,
            axis=ax
        )
        
        # Add cluster statistics to title
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0) if labels is not None else "N/A"
        n_total = len(labels) if labels is not None else "N/A"
        n_noise = sum(1 for l in labels if l == -1) if labels is not None else "N/A"
        
        ax.set_title(
            f'HDBSCAN Single Linkage Tree | {n_clusters} Clusters, {n_total} Poses, {n_noise} Noise\\n'
            f'Minimum Spanning Tree - Color gradient shows merge distance',
            fontsize=14, pad=20
        )
        ax.set_xlabel('Data Point Index', fontsize=12)
        ax.set_ylabel('Distance', fontsize=12)
        
        # Convert to base64 image
        buf = BytesIO()
        plt.tight_layout()
        plt.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='white')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode()
        plt.close(fig_mpl)
        
        # Return raw image HTML (can be right-clicked and saved)
        html = f'''
        <div style="text-align: center; background: white; padding: 10px;">
            <img src="data:image/png;base64,{img_base64}" 
                 style="max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px;"
                 alt="HDBSCAN Single Linkage Tree" />
        </div>
        '''
        return html
        
    except Exception as e:
        import traceback
        error_html = f'''
        <div style="padding: 40px; text-align: center; color: red;">
            <p><strong>Error generating linkage tree:</strong></p>
            <pre style="text-align: left; background: #f5f5f5; padding: 10px; border-radius: 4px;">{str(e)}

{traceback.format_exc()}</pre>
        </div>
        '''
        return error_html
