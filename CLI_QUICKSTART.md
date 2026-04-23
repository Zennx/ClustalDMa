# ClustalDMα CLI - Quick Start Guide

## Overview
Command-line interface for large-scale protein-nucleic acid docking analysis (handles 1,000+ models efficiently).

**Key advantages:**
- ✅ Scalable to 10,000+ models
- ✅ No Shiny app overhead
- ✅ Generates standalone HTML reports
- ✅ Batch processing friendly
- ✅ All the same analyses as GUI version

---

## Installation

No additional dependencies needed - uses same environment as Shiny app:
```bash
conda activate structural  # or your environment name
```

---

## Basic Usage

### Minimal Example
```bash
python clustal_cli.py -m models/ -o results/
```

### With Reference Sequence
```bash
python clustal_cli.py \
  -m models/ \
  -r reference.pdb \
  -o results/
```

### Full Parameter Control
```bash
python clustal_cli.py \
  -m models/ \
  -r reference.pdb \
  -o results/ \
  --protein-chains A \
  --apply-offset \
  --min-cluster-size 10 \
  --min-samples 2 \
  --hdbscan-epsilon 0.05 \
  --hdbscan-selection-method mds \
  --filter-duplicates \
  --distance-cutoff 4.5 \
  --n-jobs -1
```

---

## Required Arguments

| Argument | Description |
|----------|-------------|
| `-m, --models` | Directory containing AlphaFold-Multimer PDB models (searches recursively in subdirectories) |
| `-o, --output` | Output directory for all results |

**Note:** The `--models` argument will **recursively search** all subdirectories for `.pdb` and `.cif` files. This means you can organize your models in subfolders and the CLI will find them all automatically.

Example directory structure that works:
```
models/
├── run1/
│   ├── model_1.pdb
│   └── model_2.pdb
├── run2/
│   ├── model_1.pdb
│   └── model_2.pdb
└── other_batch/
    └── predictions/
        ├── model_1.pdb
        └── model_2.pdb
```

All 6 models will be found and analyzed!

---

## Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-r, --reference` | None | Reference PDB file for sequence alignment |
| `--protein-chains` | `A` | Protein chain ID(s), comma/space separated (GUI-aligned semantics) |
| `--apply-offset` | False | Apply residue number corrections for chopped models |
| `--min-cluster-size` | 5 | HDBSCAN minimum cluster size |
| `--min-samples` | 2 | HDBSCAN minimum samples parameter |
| `--hdbscan-epsilon` | 0.0 | HDBSCAN `cluster_selection_epsilon` (post-selection cluster merging) |
| `--hdbscan-selection-method` | `eom` | One of: `eom`, `leaf`, `mds`, `tsne`, `umap`, `pca` |
| `--filter-duplicates` | False | Filter near-duplicate structures (RMSD < 0.5 Å) |
| `--distance-cutoff` | 4.5 | Contact distance cutoff in Ångströms |
| `--n-jobs` | -1 | Number of parallel jobs (-1 = all CPUs) |

**Compatibility note:** `--ligand-id` is still accepted as a deprecated fallback, but new workflows should use `--protein-chains`.

### HDBSCAN Selection Modes

- `eom` / `leaf`: cluster directly on the precomputed distance matrix
- `mds` / `tsne` / `umap` / `pca`: first project distances into 2D, then cluster those coordinates with HDBSCAN `eom`

Example:
```bash
python clustal_cli.py \
  -m models/ \
  -o results_mds/ \
  --hdbscan-selection-method mds \
  --hdbscan-epsilon 0.02
```

---

## Output Structure

After running, your output directory will contain:

```
results/
├── index.html                      # ← Main summary report (open this first!)
├── html_reports/
│   ├── dendrogram.html            # Hierarchical clustering
│   ├── 3d_clusters.html           # Interactive 3D visualization
│   ├── jaccard_heatmap.html       # Similarity matrix
│   ├── condensed_tree.html        # HDBSCAN cluster tree
│   ├── sequence_alignment.html    # QC visualization
│   └── contact_distance_*.html    # Contact pattern analysis
├── csv_exports/
│   ├── cluster_summary.csv        # Cluster statistics
│   ├── structure_assignments.csv   # Structure → cluster mapping
│   ├── jaccard_matrix.csv         # Full similarity matrix
│   └── contact_exports/           # Per-structure contact lists
│       ├── model1_contacts.csv
│       ├── model2_contacts.csv
│       └── ...
├── cluster_representatives/
│   ├── cluster_0_medoid.pdb
│   ├── cluster_1_medoid.pdb
│   └── ...
└── plots/                          # (reserved for future exports)
```

---

## Workflow Examples

### Example 1: Quick Analysis (Default Parameters)
```bash
# Analyze 10,000 models with default settings
python clustal_cli.py \
  -m /path/to/alphafold_outputs/ \
  -o results_10k/

# Open the report
open results_10k/index.html
```

**Expected time:** ~5-10 minutes for 10,000 models (with optimizations)

---

### Example 2: High-Quality Clustering (Strict Parameters)
```bash
# Stricter clustering for high-confidence results
python clustal_cli.py \
  -m models/ \
  -r reference_structure.pdb \
  -o results_strict/ \
  --min-cluster-size 20 \
  --min-samples 5 \
  --filter-duplicates \
  --distance-cutoff 3.5
```

Use for: Finding high-confidence, reproducible binding modes

---

### Example 3: Chopped AlphaFold Models
```bash
# When using fragments of a large protein (e.g., residues 500-700)
python clustal_cli.py \
  -m chopped_models/ \
  -r full_length.pdb \
  -o results_chopped/ \
  --apply-offset
```

The `--apply-offset` flag automatically detects and corrects residue numbering from filenames.

---

### Example 4: Specify Protein Chains (GUI-Aligned)
```bash
# Use chains A and B as the protein target; all non-protein chains become partner selection
python clustal_cli.py \
  -m models/ \
  -o results/ \
  --protein-chains A,B
```

---

### Example 5: Permissive Clustering (Exploratory)
```bash
# Find even small/weak clusters
python clustal_cli.py \
  -m models/ \
  -o results_exploratory/ \
  --min-cluster-size 3 \
  --min-samples 1 \
  --distance-cutoff 5.0
```

Use for: Exploratory analysis, finding rare binding modes

---

## Performance Tips

### For 10,000+ Models:

1. **Use all CPUs:**
   ```bash
   --n-jobs -1
   ```

2. **Filter duplicates** to reduce dataset size:
   ```bash
   --filter-duplicates
   ```

3. **Run on HPC/cluster** for very large datasets:
   ```bash
   sbatch run_clustal.sh  # Example SLURM script
   ```

4. **Split analysis** by subdirectories if needed:
   ```bash
   python clustal_cli.py -m models_batch1/ -o results1/
   python clustal_cli.py -m models_batch2/ -o results2/
   ```

### Expected Runtimes (approximate):

| Models | Time (optimized) | Peak RAM |
|--------|------------------|----------|
| 100    | ~10 seconds      | ~500 MB  |
| 1,000  | ~1 minute        | ~2 GB    |
| 10,000 | ~10 minutes      | ~10 GB   |
| 50,000 | ~1 hour          | ~30 GB   |

*With KDTree optimization and parallelization*

---

## Viewing Results

### Main Report
```bash
# macOS
open results/index.html

# Linux
xdg-open results/index.html

# Windows
start results/index.html
```

The main report (`index.html`) contains:
- Analysis parameters summary
- Cluster statistics table
- Links to all interactive visualizations
- List of exported files

### Interactive Visualizations
All HTML reports are **fully interactive** with:
- Zoom, pan, hover tooltips
- Downloadable as PNG/SVG
- No server needed - works offline
- Shareable with collaborators

---

## Comparison: CLI vs Shiny App

| Feature | Shiny App | CLI |
|---------|-----------|-----|
| **Max models** | ~1,000 | 50,000+ |
| **Interface** | GUI (browser) | Command-line |
| **Real-time updates** | ✅ Yes | ❌ No |
| **Batch processing** | ❌ No | ✅ Yes |
| **HPC compatible** | ❌ No | ✅ Yes |
| **Output** | Temporary | Persistent HTML |
| **Resource usage** | High (server running) | Low (one-time) |

**Use Shiny app for:**
- Interactive exploration
- Small-medium datasets (<1,000 models)
- Real-time parameter tuning

**Use CLI for:**
- Large datasets (1,000+ models)
- Production pipelines
- HPC/cluster environments
- Reproducible analyses

---

## Troubleshooting

### "No PDB files found"
- Check that `-m` points to directory with `.pdb` files
- Verify file extensions are `.pdb` (not `.pdb.gz`)

### "Out of memory"
- Reduce dataset size with `--filter-duplicates`
- Use smaller `--min-cluster-size` to reduce matrix operations
- Run on machine with more RAM

### "Analysis taking too long"
- Verify KDTree optimization is active (check code in `core/metrics.py`)
- Reduce `--n-jobs` if causing overhead
- Consider splitting dataset

### "Visualizations not generating"
- Check console output for specific errors
- Some visualizations require reference sequence
- Contact distance maps need valid clusters

---

## Advanced: SLURM Script Template

```bash
#!/bin/bash
#SBATCH --job-name=clustal_analysis
#SBATCH --output=clustal_%j.log
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4:00:00

module load anaconda3
conda activate structural

python clustal_cli.py \
  -m /scratch/models/ \
  -r /scratch/reference.pdb \
  -o /scratch/results_${SLURM_JOB_ID}/ \
  --min-cluster-size 10 \
  --filter-duplicates \
  --n-jobs ${SLURM_CPUS_PER_TASK}

echo "Results saved to: /scratch/results_${SLURM_JOB_ID}/"
```

---

## Getting Help

```bash
# Show all options
python clustal_cli.py --help
```

For issues or questions:
- Check the main [README.md](README.md)
- Review [USAGE.md](USAGE.md) for concepts
- Examine log output for error messages

---

## Next Steps

After analysis completes:

1. **Review main report:** `results/index.html`
2. **Explore clusters:** Check 3D visualization and dendrogram
3. **Validate contacts:** Review sequence alignment QC
4. **Export data:** Use CSV files for downstream analysis
5. **Share results:** Entire `results/` folder is self-contained

Happy clustering! 🧬
