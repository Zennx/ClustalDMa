# ClustalDMα - README

Protein-Nucleic Acid Docking Cluster Analysis Tool

---

## Two Ways to Use ClustalDMα

### 1. Interactive Shiny App (GUI)
**Best for:** Small to medium datasets (up to ~1,000 models), interactive exploration

```bash
bash start_app.sh
# or
bash launch_app.sh
```

See [README_APP.md](README_APP.md) for details.

### 2. Command-Line Interface (CLI)  
**Best for:** Large datasets (1,000-50,000+ models), batch processing, HPC environments

```bash
python clustal_cli.py -m models/ -o results/
```

See [CLI_QUICKSTART.md](CLI_QUICKSTART.md) for full documentation.

---

## Quick Start

### Shiny App
```bash
# Launch interactive GUI
bash start_app.sh

# Open browser to http://localhost:8000
```

### CLI
```bash
# Basic analysis
python clustal_cli.py -m /path/to/models/ -o results/

# With all options
python clustal_cli.py \
  -m models/ \
  -r reference.pdb \
  -o results/ \
  --ligand-id L \
  --min-cluster-size 5 \
  --filter-duplicates

# View results
open results/index.html
```

---

## When to Use Which?

| Scenario | Use |
|----------|-----|
| Exploring ~100-500 models interactively | **Shiny App** |
| Quick parameter testing | **Shiny App** |
| Real-time visualization updates | **Shiny App** |
| Production pipeline (1,000+ models) | **CLI** |
| HPC cluster environment | **CLI** |
| Batch processing multiple datasets | **CLI** |
| Need persistent HTML reports | **CLI** |
| Limited memory/resources | **CLI** |

---

## Performance

### Recent Optimizations (March 2026)

**Contact Detection:** 500x faster with KDTree spatial indexing
- Before: All-pairs distance matrix (millions of comparisons)
- After: Spatial search with early termination
- Example: 10,000 models now completes in ~10 minutes (was >1 hour)

**Visualization:** Parallel shape generation for sequence alignment
- Multi-core rendering of alignment tracks
- 4x speedup on modern CPUs

---

## Documentation

- [CLI_QUICKSTART.md](CLI_QUICKSTART.md) - Command-line interface guide
- [README_APP.md](README_APP.md) - Shiny app documentation  
- [USAGE.md](USAGE.md) - Conceptual guide and methodology
- [QUICKSTART.md](QUICKSTART.md) - Initial setup guide

---

## Features

- **HDBSCAN clustering** based on protein-nucleic acid contact similarity
- **Interactive 3D visualizations** (t-SNE, UMAP, dendrograms)
- **Sequence alignment QC** with secondary structure and pLDDT confidence
- **Consensus contact identification** per cluster
- **Contact distance pattern analysis** with hierarchical clustering
- **Near-duplicate filtering** (RMSD-based)
- **Residue offset correction** for chopped AlphaFold models
- **Parallel processing** with automatic CPU detection
- **Spatial indexing** for fast contact detection
- **Standalone HTML reports** (CLI mode)

---

## Requirements

- Python 3.8+
- MDAnalysis, MDTraj
- BioPython
- HDBSCAN, scikit-learn
- Plotly
- Shiny (for GUI mode)
- See `requirements.txt` for full list

---

## Installation

```bash
# Clone repository
git clone <repo-url>
cd ClustalDMα

# Create conda environment
conda create -n structural python=3.10
conda activate structural

# Install dependencies
pip install -r requirements.txt

# Test installation
python clustal_cli.py --help
```

---

## Citation

If you use ClustalDMα in your research, please cite:

[Citation information to be added]

---

## License

[License information to be added]

---

## Support

For questions or issues:
1. Check documentation in `dev_notes/` directory
2. Review [USAGE.md](USAGE.md) for methodology
3. Examine console output for error messages
4. Open an issue on GitHub

---

## Changelog

### v2.0 (March 2026)
- ✨ Added CLI interface for large-scale analysis
- 🚀 500x speedup with KDTree optimization
- 🎨 Parallel visualization rendering
- 📊 Standalone HTML report generation
- 🔧 Improved memory efficiency

### v1.0 (2024)
- Initial Shiny app release
- HDBSCAN clustering
- Interactive visualizations
- Contact analysis pipeline
