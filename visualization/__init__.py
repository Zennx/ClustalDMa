"""
Visualization Module
Interactive plots using Plotly and Molstar 3D viewers
"""

from .interactive import (
    create_interactive_scatter,
    create_hotspot_histogram,
    create_mol_viewer_html,  # Legacy py3Dmol (deprecated)
    create_molstar_viewer_html,  # Modern Molstar viewer
    create_scatter_multimethod,
    create_distance_heatmap,
    create_cluster_size_distribution,
)

__all__ = [
    'create_interactive_scatter',
    'create_hotspot_histogram', 
    'create_mol_viewer_html',
    'create_molstar_viewer_html',
    'create_scatter_multimethod',
    'create_distance_heatmap',
    'create_cluster_size_distribution',
]
