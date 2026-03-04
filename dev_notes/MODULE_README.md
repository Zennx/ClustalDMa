# ClustalDM - Modular Structure

## 📁 Project Organization

```
ClustalDM/
├── core/                      # Core clustering functionality
│   ├── __init__.py
│   ├── clusterer.py          # Main PDBClusterer class
│   ├── metrics.py            # Distance metrics (RMSD, Jaccard, etc.)
│   ├── analysis.py           # Interface analysis & hotspots
│   └── io_utils.py           # File finding & I/O utilities
│
├── visualization/             # Visualization components
│   ├── __init__.py
│   └── interactive.py        # Plotly & py3Dmol interactive plots
│
├── app.py                    # Main Shiny web app
├── app_interactive.py        # Interactive features app (NEW!)
├── jaccard_prototype.py      # Standalone Jaccard test script
├── clustaldemo.py           # Legacy file (for backwards compatibility)
└── requirements.txt          # Python dependencies
```

## 🚀 Quick Start

### Using the Modular API

```python
from core import PDBClusterer, find_pdb_files
from visualization import create_interactive_scatter, create_hotspot_histogram

# Find PDB files
pdb_files = find_pdb_files("/path/to/hdock/results", filter_pattern='model')

# Create clusterer
clust = PDBClusterer(pdb_files)

# Compute Jaccard distance matrix
clust.compute_jaccard_contact_matrix(distance_cutoff=4.5)

# Cluster
clust.cluster_dbscan(eps=0.3, min_samples=2)

# Get statistics
stats_df = clust.get_interface_stats()
hotspots_df = clust.get_binding_hotspots()

# Create interactive plots
fig_scatter = create_interactive_scatter(stats_df, clust.distance_matrix)
fig_hotspots = create_hotspot_histogram(hotspots_df)
```

### Running the Interactive App

```bash
conda activate structural
shiny run app_interactive.py
```

Then open your browser and:
1. Enter PDB directory path
2. Click "Run Jaccard Clustering"
3. Explore interactive visualizations!

## 🎯 Key Features

### Core Modules

**clusterer.py** - Main clustering class
- `PDBClusterer`: Main clustering interface
- Methods: `compute_distance_matrix()`, `compute_jaccard_contact_matrix()`, `cluster_dbscan()`, `cluster_hierarchical()`

**metrics.py** - Distance computation
- `DistanceMetrics.compute_rmsd_matrix()`: RMSD-based distances
- `DistanceMetrics.compute_jaccard_matrix()`: Jaccard contact similarity
  - **Key Innovation**: Compares protein residues only, not residue pairs
  - This makes binding modes comparable even with different DNA sequences!

**analysis.py** - Interface characterization
- `InterfaceAnalyzer.get_interface_stats()`: Per-structure statistics
- `InterfaceAnalyzer.get_binding_hotspots()`: Identify hotspot residues
- `InterfaceAnalyzer.get_cluster_signature()`: Consensus interface residues

**io_utils.py** - File management
- `find_pdb_files()`: Smart PDB file discovery with filtering

### Visualization Modules

**interactive.py** - Interactive plots
- `create_interactive_scatter()`: PCA projection with hover details
- `create_hotspot_histogram()`: Binding hotspot bar chart
- `create_mol_viewer_html()`: 3D structure viewer (centered on DNA!)

## 💡 Design Benefits

1. **Separation of Concerns**: Each module has a single, well-defined responsibility
2. **Testability**: Easy to test individual components
3. **No Redundancy**: Shared code in one place
4. **Flexibility**: Import only what you need
5. **Backwards Compatible**: Old imports still work via `clustaldemo.py`
6. **Maintainability**: Clear structure for future development

## 🔬 Technical Details

### Jaccard Contact Distance

The Jaccard metric compares **protein residues in contact**, not protein-nucleic pairs:

**Why?** Because DNA shape is similar regardless of sequence:
- `ARG123-DG5` and `ARG123-DA5` → Same binding mode (both use ARG123)
- Different nucleotide sequence shouldn't affect binding mode clustering

**Formula:**
```
distance = 1 - (|A ∩ B| / |A ∪ B|)
where A, B = sets of protein residues in contact
```

### Distance Matrix Convention

- **0.0 = identical** (purple in heatmap, same location)
- **1.0 = completely different** (yellow in heatmap, far apart)
- This follows standard distance matrix conventions

## 📊 Example Workflow

```python
# 1. Setup
from core import PDBClusterer, find_pdb_files

files = find_pdb_files("./hdock_results", filter_pattern='model')
clust = PDBClusterer(files)

# 2. Compute distances
clust.compute_jaccard_contact_matrix(distance_cutoff=4.5)

# 3. Cluster
clust.cluster_dbscan(eps=0.3, min_samples=2)

# 4. Analyze
stats = clust.get_interface_stats()
hotspots = clust.get_binding_hotspots()
signature = clust.get_cluster_signature(cluster_id=0)

print(f"Cluster 0 signature residues: {signature}")
```

## 🐛 Troubleshooting

**Import errors?**
- Make sure you're in the `ClustalDM` directory
- Python needs to find the `core/` and `visualization/` packages

**No structures loaded?**
- Use `find_pdb_files()` helper to locate files
- Check filter_pattern matches your filenames

**Clustering fails?**
- Make sure to run a distance metric first
- For Jaccard: `compute_jaccard_contact_matrix()`
- For RMSD: `load_structures()` then `compute_distance_matrix()`

## 📝 Notes

- The old `clustaldemo.py` file remains for backwards compatibility
- New code should use the modular imports from `core/` and `visualization/`
- All distance matrices use the "precomputed" format for clustering

##✨ New in This Version

- **Modular architecture**: Clean separation into focused modules
- **Interactive visualizations**: Plotly scatter plots, histograms, 3D viewer
- **Protein-only Jaccard**: Fixed to compare protein residues, not pairs
- **Smart file finding**: Automatic filtering for model files
- **Enhanced analysis**: Interface stats, hotspots, cluster signatures
