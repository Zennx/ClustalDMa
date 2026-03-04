# 🧬 ClustalDM - Interactive HDOCK Clustering Analysis

An interactive web application for clustering and analyzing HDOCK protein-nucleic acid docking results using DBSCAN and structural RMSD.

## ✨ Features

- **Interactive Web Interface**: User-friendly Shiny app for parameter tuning
- **Real-time Clustering**: Adjust DBSCAN parameters and see results instantly
- **Multiple Visualizations**:
  - Distance matrix heatmaps
  - 2D projections (MDS, PCA, t-SNE)
  - Hierarchical clustering dendrograms
  - Cluster size distributions
- **Quality Metrics**: Silhouette scores and intra-cluster statistics
- **PyMOL Integration**: Generate PyMOL scripts for visual inspection
- **Efficient Export**: Save results with optional symbolic links

## 🚀 Quick Start

### Installation

```bash
# Clone or navigate to the repository
cd ClustalDM

# Install dependencies
pip install -r requirements.txt
```

### Launch the App

```bash
# Option 1: Using the launcher script
python run_app.py

# Option 2: Direct shiny command
shiny run app.py
```

The app will open in your default browser at `http://127.0.0.1:8000`

### Using the CLI (Original)

```bash
python clustaldemo.py -d /path/to/HDOCK/output \
    -s "nucleic and name P" \
    --no-align \
    -e 12.0 \
    -m 2 \
    -o output_dir \
    --use-symlinks
```

## 📖 Usage Guide

### Web App Workflow

1. **Input Settings**
   - Enter your HDOCK output directory path
   - Specify atom selection (e.g., `nucleic and name P` for DNA phosphate backbone)
   - Set filename filter (default: `model`)

2. **Clustering Parameters**
   - Adjust DBSCAN `eps` (distance threshold in Ångströms)
   - Set `min_samples` (minimum cluster size)
   - Toggle alignment (disable for docking results)

3. **Run Analysis**
   - Click "🚀 Run Clustering Analysis"
   - Wait for processing (progress shown in status panel)

4. **Explore Results**
   - **Overview**: Cluster summary and quality metrics
   - **Distance Matrix**: RMSD heatmap sorted by clusters
   - **Scatter Plots**: 2D projections with different methods
   - **Dendrogram**: Hierarchical clustering tree
   - **Cluster Details**: Examine individual clusters
   - **Export**: Save all results to disk

### Recommended Settings for HDOCK

- **Selection**: `nucleic and name P` (for DNA/RNA ligand)
- **No Alignment**: ✅ Checked (preserves docking poses)
- **eps**: 10-15 Å (adjust based on ligand flexibility)
- **min_samples**: 2-3 structures

## 📊 Output Files

When you export results, the following structure is created:

```
output_dir/
├── distance_matrix.png           # RMSD heatmap
├── distance_matrix_clustered.png # Sorted by cluster
├── dbscan_clusters.png           # Cluster sizes
├── cluster_dendrogram.png        # Hierarchical tree
├── dbscan_scatter.png            # MDS, t-SNE, PCA projections
├── dendrogram_labeled.png        # Model numbers labeled
├── cluster_info.txt              # Text summary
├── representatives/              # Centroid structures
│   ├── cluster_0_representative.pdb
│   ├── cluster_1_representative.pdb
│   └── ...
└── cluster_pdbs/                 # Organized by cluster
    ├── load_cluster_0.pml        # PyMOL script
    ├── load_cluster_1.pml
    ├── cluster_0/
    │   ├── model_1.pdb (or symlink)
    │   ├── model_2.pdb
    │   └── ...
    └── ...
```

## 🔬 PyMOL Visualization

After exporting, visualize clusters in PyMOL:

```bash
cd output_dir/cluster_pdbs
pymol load_cluster_0.pml
```

The script will:
- Load all structures with unique names (e.g., `cluster_0_model_1`)
- Color protein (cyan) and DNA (orange)
- Set up cartoon + stick representation
- Provide alignment commands (commented out)

In PyMOL, you can:
- Toggle structures: `enable cluster_0_model_1`
- Align structures: `alignto cluster_0_model_1`
- Compare poses side-by-side

## 🛠️ Advanced Features

### Command-Line Options

```
-d, --directories    : Input directories with PDB files
-s, --selection      : MDAnalysis atom selection string
-e, --eps            : DBSCAN epsilon (distance threshold)
-m, --min-samples    : DBSCAN minimum samples
-c, --chain          : Auto-select chain (alternative to -s)
-o, --output-dir     : Output directory
-f, --filter         : Filename filter pattern
--no-align           : Skip alignment before RMSD
--align-selection    : Custom selection for alignment
--use-symlinks       : Use symbolic links instead of copying
```

### MDAnalysis Selection Syntax

Examples:
- `nucleic and name P` - DNA/RNA phosphate atoms
- `protein and name CA` - Protein alpha carbons
- `resname LIG` - Specific ligand residue
- `protein and segid A` - Specific protein chain

## 📝 Citation

If you use ClustalDM in your research, please cite:
- HDOCK: Yan et al., Nature Protocols (2020)
- MDAnalysis: Michaud-Agrawal et al., J. Comput. Chem. (2011)
- DBSCAN: Ester et al., KDD (1996)

## 🐛 Troubleshooting

**App won't start:**
```bash
pip install --upgrade shiny
```

**Import errors:**
```bash
pip install -r requirements.txt --upgrade
```

**No clusters found:**
- Increase `eps` value
- Decrease `min_samples`
- Check your selection string is valid

**PyMOL can't find files:**
- Use full paths in PyMOL scripts
- Or `cd` to the cluster directory before loading

## 📄 License

MIT License - See LICENSE file for details

## 🤝 Contributing

Contributions welcome! Please open an issue or pull request.

---

**Developed for interactive analysis of HDOCK molecular docking results** 🧬
