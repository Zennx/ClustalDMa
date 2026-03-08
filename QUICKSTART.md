# Quick Start Guide - ClustalDM for AlphaFold2

## Installation

```bash
pip install -r requirements.txt
```

## Launch the App

```bash
./launch_app.sh
# or
shiny run app_main.py
```

## Basic Usage

1. **Specify your AlphaFold2 models directory**
   - Enter the full path to a directory containing PDB or mmCIF files
   - Example: `/path/to/alphafold/predictions`

2. **(Optional) Upload a reference model**
   - Click "Choose File" under "Reference Model"
   - Select a PDB or CIF file to use for superposition/RMSD

3. **(Optional) Define motif residues**
   - Format: `A:10-20,B:30-35`
   - Leave blank for global interface analysis

4. **Adjust clustering parameters**
   - `min_cluster_size`: Minimum structures per cluster (default: 5)
   - `min_samples`: Controls noise sensitivity (default: 2)
   - `Filter Near-Duplicates`: Remove very similar structures

5. **Run Analysis**
   - Click "🚀 Run Analysis"
   - Monitor progress in the status log

6. **Explore Results**
   - Navigate through tabs: Overview, Distance Matrices, Interactive Maps, Cluster Trees, 3D Viewer
   - Download results using "📥 Download Results" button

## Key Features

- **Jaccard Contact Similarity**: Cluster by interface residue contacts (robust to conformational changes)
- **RMSD Validation**: Quality control using Ca-atom RMSD
- **HDBSCAN Clustering**: Density-based clustering that finds optimal number of clusters
- **Interactive 3D Viewer**: Molstar integration with pLDDT coloring support
- **Reference Model Support**: Align all structures to a reference for consistent RMSD calculation

## Tips

- For large datasets (100+ models), enable "Filter Near-Duplicates"
- Use motif residues if you want to focus on specific binding sites
- Check the HDBSCAN condensed tree to evaluate cluster stability
- Export cluster PDBs for further analysis in external tools
