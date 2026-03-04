# PDB Clustering Script - Usage Guide

## Overview
This script performs DBSCAN clustering on PDB files using MDAnalysis and generates distance maps and visualization plots.

## Features
1. **Accepts multiple directories** containing PDB files
2. **DBSCAN clustering** based on RMSD distances
3. **Generates multiple plots**:
   - Distance matrix heatmap
   - Clustered distance matrix (sorted by clusters)
   - Cluster assignments bar plot
   - Cluster size distribution
   - Dendrogram-style visualization
4. **Saves detailed cluster information** to text file

## Requirements
Install required packages:
```bash
pip install MDAnalysis numpy scikit-learn matplotlib scipy
```

Optional (for enhanced plots):
```bash
pip install seaborn
```

## Usage

### Basic Usage
```bash
python clustaldemo.py -d /path/to/pdb/files
```

### Specify DBSCAN Parameters
```bash
python clustaldemo.py -d /path/to/pdb/files -e 2.5 -m 3
```

### Multiple Directories
```bash
python clustaldemo.py -d /path/to/dir1 /path/to/dir2 /path/to/dir3
```

### Custom Selection String
```bash
python clustaldemo.py -d /path/to/pdb/files -s "protein and name CA" -e 2.0
```

### Specify Output Directory
```bash
python clustaldemo.py -d /path/to/pdb/files -o my_results
```

### Filter for Specific Files (e.g., HDOCK models)
```bash
# Only process files with "model" in the name (e.g., model_1.pdb, model_2.pdb)
python clustaldemo.py -d /path/to/hdock/output -f model

# Process all PDB files (no filter)
python clustaldemo.py -d /path/to/pdb/files -f ""
```

### Specify Chain for Clustering (Important for HDOCK!)
```bash
# Cluster only chain B (protein ligand chain) - for protein-protein docking
python clustaldemo.py -d /path/to/hdock/output -c B -e 3.0

# Cluster DNA/RNA ligand - for protein-nucleic acid docking
python clustaldemo.py -d /path/to/hdock/output -s "nucleic" -e 3.0

# Cluster only phosphate backbone of nucleic acids (faster)
python clustaldemo.py -d /path/to/hdock/output -s "nucleic and name P" -e 3.0
```

## Command Line Arguments

| Argument | Short | Description | Default |
|----------|-------|-------------|---------|
| `--directories` | `-d` | Directory(ies) containing PDB files (required) | - |
| `--eps` | `-e` | DBSCAN eps parameter (max distance in Å) | 2.0 |
| `--min-samples` | `-m` | DBSCAN min_samples parameter | 2 |
| `--selection` | `-s` | MDAnalysis selection string | Auto-detect based on --chain |
| `--chain` | `-c` | Specific chain ID to analyze (e.g., "B" for ligand) | None (all chains) |
| `--output-dir` | `-o` | Output directory for results | clustering_output |
| `--filter` | `-f` | Filter pattern for PDB filenames | "model" |
| `--no-align` | - | Skip structure alignment (recommended for docking results) | False (aligns by default) |
| `--align-selection` | - | Custom selection for alignment (e.g., "protein" to align on receptor) | Same as --selection |

## DBSCAN Parameters

### eps (epsilon)
- The maximum distance between two samples for them to be in the same neighborhood
- For protein structures, typical values range from 1.5 to 5.0 Å
- **Smaller eps** → more clusters, stricter grouping
- **Larger eps** → fewer clusters, looser grouping

### min_samples
- The minimum number of samples in a neighborhood for a point to be a core point
- Typical values: 2-5
- **Smaller min_samples** → more points classified as core points
- **Larger min_samples** → stricter cluster membership

## Selection Strings

Common MDAnalysis selection strings:

- `"protein and name CA"` - Only C-alpha atoms (default, fast)
- `"protein"` - All protein atoms (slower, more accurate)
- `"backbone"` - Backbone atoms only
- `"name CA C N O"` - Specific atom types
- `"resid 1-100"` - Specific residue range
- `"nucleic"` - All nucleic acid (DNA/RNA) atoms
- `"nucleic and name P"` - Only phosphate atoms of nucleic acids (very fast)
- `"nucleic and backbone"` - Nucleic acid backbone atoms

## Output Files

The script generates the following files in the output directory:

### Plots
1. **distance_matrix.png** - Heatmap of pairwise RMSD distances
2. **distance_matrix_clustered.png** - Distance matrix sorted by cluster with boundaries
3. **dbscan_clusters.png** - Cluster assignments and size distribution
4. **cluster_dendrogram.png** - Hierarchical dendrogram and intra-cluster RMSD
5. **dbscan_scatter.png** - DBSCAN visualization with MDS/PCA projections, RMSD distributions, and cluster statistics
6. **dendrogram_labeled.png** - Hierarchical dendrogram with model numbers as labels (color-coded by DBSCAN cluster)

### Text Output
7. **cluster_info.txt** - Detailed text file listing cluster assignments with model numbers

### Validation & Inspection
8. **representatives/** - Directory containing representative PDB file for each cluster (centroid structure)
9. **cluster_models/** - Directory containing multi-model PDB files (one per cluster) for visualization in PyMOL/ChimeraX

## Example Workflow

```bash
# 1. Navigate to your working directory
cd /Users/Zen/Documents/ClustalDM

# 2. For HDOCK with DNA/RNA ligand (e.g., ssDNA, dsRNA)
python clustaldemo.py \
  -d /path/to/hdock/output \
  -s "nucleic and name P" \
  -f model \
  --no-align \
  -e 3.0 \
  -m 2 \
  -o results_dna_eps3

# 3. For HDOCK with protein ligand
python clustaldemo.py \
  -d /path/to/hdock/output \
  -c B \
  -f model \
  --no-align \
  -e 3.0 \
  -m 2 \
  -o results_protein_eps3

# 4. Check the results
ls results_dna_eps3/
cat results_dna_eps3/cluster_info.txt

# 5. If you want to try different eps values, run again
python clustaldemo.py -d /path/to/hdock/output -s "nucleic and name P" --no-align -e 5.0 -o results_eps5
```

## HDOCK-Specific Usage

When working with HDOCK output, which typically contains:
- `model_1.pdb`, `model_2.pdb`, ..., `model_N.pdb` (docking predictions - complete complexes)
- `receptor.pdb` (receptor structure alone)
- `ligand.pdb` or `lig_*.pdb` (ligand structure alone)

### Important: Selecting the Right Molecule for Clustering

HDOCK keeps the **receptor fixed** and varies the **ligand** position across different models.

**⚠️ You MUST select only the ligand for clustering!**

#### Case 1: Protein-Nucleic Acid Docking (DNA/RNA ligand)

If your ligand is any nucleic acid (ssDNA, dsDNA, sRNA, dsRNA, mixed RNA/DNA, etc.):

```bash
# CORRECT: Cluster nucleic acid ligand positions
python clustaldemo.py -d /path/to/hdock/output -s "nucleic" --no-align -e 3.0

# Faster (using only phosphate backbone atoms) - RECOMMENDED for large datasets:
python clustaldemo.py -d /path/to/hdock/output -s "nucleic and name P" --no-align -e 3.0

# Full recommended command for nucleic acid ligands:
python clustaldemo.py \
  -d /path/to/hdock/output \
  -s "nucleic and name P" \
  -f model \
  --no-align \
  -e 3.0 \
  -m 2 \
  -o hdock_nucleic_clusters
```

**Important**: Use `--no-align` flag for docking results! This measures the actual spatial differences in ligand positions. Without it, structures are aligned which defeats the purpose of clustering docking poses.

**Note**: Using `"nucleic and name P"` (phosphate atoms only) is recommended because:
- Much faster computation
- Captures the overall backbone conformation
- Works for any nucleic acid type (DNA, RNA, hybrids)
- Reduces noise from base flexibility

#### Case 2: Protein-Protein Docking

If your ligand is a protein (usually Chain B):

```bash
# CORRECT: Cluster protein ligand (Chain B)
python clustaldemo.py -d /path/to/hdock/output -c B --no-align -e 3.0

# Full recommended command for protein-protein:
python clustaldemo.py \
  -d /path/to/hdock/output \
  -c B \
  -f model \
  --no-align \
  -e 3.0 \
  -m 2 \
  -o hdock_protein_clusters
```

**❌ WRONG: This will cluster the receptor (always same) → RMSD = 0**
```bash
python clustaldemo.py -d /path/to/hdock/output  # Uses default "protein and name CA" - includes receptor!
```

### How to identify your ligand type:

Check the ligand file or model file structure:

```bash
# Check what's in the ligand file
head -10 lig_*.pdb

# If you see residues like "A", "G", "C", "T", "U", "DA", "DG", "DC", "DT":
#   → It's DNA/RNA, use: -s "nucleic"

# If you see residues like "ALA", "GLY", "SER", "VAL", etc.:
#   → It's protein, use: -c B (or whichever chain is the ligand)

# Count atoms in each chain of model files
grep "^ATOM" model_1.pdb | awk '{print $5}' | sort | uniq -c
#  1472 A  <- Usually receptor (larger) + DNA (if present)
#   680 B  <- Usually protein ligand (if protein-protein docking)

# Check for nucleic acids in the file
grep "^ATOM" model_1.pdb | grep " P " | head -5
# If you see phosphate (P) atoms → DNA/RNA is present
```

The script now **defaults to filtering for "model"** in filenames, so it will automatically skip the receptor and ligand files.

### Quick Start for HDOCK:

```bash
# For any nucleic acid ligand (ssDNA, dsDNA, RNA, etc.) - RECOMMENDED:
python clustaldemo.py -d /path/to/hdock/output -s "nucleic and name P" --no-align -e 3.0

# For protein ligand:
python clustaldemo.py -d /path/to/hdock/output -c B --no-align -e 3.0

# For mixed analysis (align on receptor, measure ligand RMSD):
python clustaldemo.py -d /path/to/hdock/output -s "nucleic and name P" --align-selection "protein" -e 3.0
```

**Tip**: If all structures end up in one cluster, your `eps` value is too large. Check the `dbscan_scatter.png` plot to see the RMSD distribution and adjust `eps` accordingly (try using a value close to the mean or median RMSD).

## Validating and Inspecting Clusters

After clustering, you can validate the results using several methods:

### 1. Check Model Numbers in cluster_info.txt
Open `cluster_info.txt` to see which model numbers belong to each cluster:
```
Cluster 0 (25 structures):
  Model numbers: 1, 5, 12, 23, 34, 45, ...
```
You can then manually inspect specific models in PyMOL or your preferred viewer.

### 2. View Labeled Dendrogram
The `dendrogram_labeled.png` shows model numbers on the x-axis, color-coded by DBSCAN cluster. This helps you:
- See which models cluster together
- Identify outliers
- Validate clustering decisions

### 3. Inspect Representative Structures
The `representatives/` folder contains one PDB file per cluster (the centroid structure):
```bash
ls results/representatives/
# cluster_0_representative.pdb
# cluster_1_representative.pdb
# cluster_2_representative.pdb
```
These are the "most typical" structures from each cluster - great for quick inspection!

### 4. Visualize All Cluster Members in PyMOL/ChimeraX
The `cluster_models/` folder contains multi-model PDB files (one per cluster):
```bash
# In PyMOL
load results/cluster_models/cluster_0_multimodel.pdb
# All structures in cluster 0 will load as different states
# Use the state slider to browse through them

# Or align and visualize all at once
intra_fit cluster_0_multimodel
show cartoon
spectrum count, rainbow

# In ChimeraX
open results/cluster_models/cluster_0_multimodel.pdb
tile
```

This lets you visually confirm that structures in the same cluster are indeed similar!

## Interpreting Results

### Silhouette Score
- Range: -1 to 1
- **> 0.5**: Good clustering
- **0.2 - 0.5**: Reasonable clustering
- **< 0.2**: Poor clustering, consider adjusting parameters

### Cluster Labels
- **Cluster 0, 1, 2...**: Valid clusters
- **Cluster -1**: Noise points (outliers)

### Distance Matrix
- **Dark colors**: Low RMSD (similar structures)
- **Light colors**: High RMSD (different structures)
- **Block patterns**: Indicate distinct structural clusters

## Tips for Parameter Selection

1. **Start with default parameters** and examine the results
2. **Check the RMSD distribution** in `dbscan_scatter.png` - the histogram shows you the range of RMSD values
3. **Set eps based on your data**:
   - Look at the mean/median RMSD in the histogram
   - For meaningful clusters, set `eps` to 0.5-1.0× the median RMSD
   - If all points are in one cluster → decrease `eps`
   - If all points are noise → increase `eps`
4. **If too many noise points**: Increase `eps` or decrease `min_samples`
5. **If too few clusters**: Decrease `eps`
6. **If too many small clusters**: Increase `eps` or increase `min_samples`
7. **Use the MDS/PCA scatter plots** to visualize structure similarity in 2D

## Programmatic Usage

You can also use the `PDBClusterer` class in your own scripts:

```python
from clustaldemo import PDBClusterer, find_pdb_files

# Find PDB files
pdb_files = find_pdb_files(['/path/to/directory'])

# Initialize clusterer
clusterer = PDBClusterer(pdb_files, selection='protein and name CA')

# Load and compute distances
clusterer.load_structures()
clusterer.compute_distance_matrix()

# Cluster
labels = clusterer.cluster_dbscan(eps=2.5, min_samples=3)

# Generate plots
clusterer.plot_distance_matrix('my_distance_matrix.png')
clusterer.plot_clusters('my_clusters.png')
clusterer.save_cluster_info('my_cluster_info.txt')

# Access results
print(f"Cluster labels: {clusterer.labels}")
print(f"Distance matrix shape: {clusterer.distance_matrix.shape}")
```

## Troubleshooting

### "No module named 'MDAnalysis'"
```bash
pip install MDAnalysis
```

### "Could not load PDB file"
- Check that PDB files are valid and not corrupted
- Ensure files have .pdb extension
- Check file permissions

### Memory Issues
- Use C-alpha atoms only (`-s "protein and name CA"`)
- Process smaller batches of structures
- Increase system swap space

### All points are noise
- Increase the `eps` parameter
- Check RMSD range in output and set eps accordingly
- Decrease `min_samples` parameter

## Contact & Support

For issues or questions, check:
- MDAnalysis documentation: https://docs.mdanalysis.org/
- scikit-learn DBSCAN: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
