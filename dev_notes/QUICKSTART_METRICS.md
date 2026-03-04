# Quick Start Guide: Advanced Distance Metrics

## Phase 1 Features Now Available! 🎉

ClustalDM now supports **5 distance metrics** for clustering HDOCK protein-nucleic acid docking results.

---

## Using the Web App

### Launch the app:
```bash
./launch_app.sh
```

### Select Your Distance Metric:

In the sidebar, you'll see a **"Distance Metric"** dropdown with these options:

---

### 1️⃣ RMSD (Root Mean Square Deviation)
**When to use**: Standard structural alignment comparison  
**Best for**: Overall backbone similarity  
**Speed**: ⚡⚡⚡ Very Fast  
**Tip**: Default and most commonly used

**Example Settings**:
- eps: 8-15 Å (depending on your structures)
- No alignment: Checked (for docking results)

---

### 2️⃣ Jaccard Contact Score
**When to use**: Comparing binding interfaces  
**Best for**: Finding similar contact patterns  
**Speed**: ⚡⚡ Fast  
**Tip**: Great for identifying similar binding modes even with different poses

**What it does**: 
- Identifies protein-nucleic contacts (residue pairs < 4.5 Å)
- Compares contact sets between structures
- Distance = 1 - (shared contacts / total contacts)

**Example Settings**:
- eps: 0.3-0.7 (normalized scale)
- min_samples: 2-5

---

### 3️⃣ Contact Area Distance (CAD)
**When to use**: Comparing binding surface areas  
**Best for**: Finding poses with similar interface sizes  
**Speed**: ⚡⚡⚡ Very Fast  
**Tip**: Good for filtering by "buried surface area"

**What it does**:
- Calculates SASA (Solvent Accessible Surface Area)
- Compares deltaSASA between structures
- Identifies similar contact areas

**Example Settings**:
- eps: 0.2-0.5 (normalized scale)
- min_samples: 2-5

---

### 4️⃣ H-bond Pattern Similarity
**When to use**: Focus on electrostatic interactions  
**Best for**: Grouping by hydrogen bonding networks  
**Speed**: ⚡⚡ Fast  
**Tip**: Excellent for RNA/DNA binding analysis

**What it does**:
- Detects H-bonds using Baker-Hubbard criterion
- Compares donor-acceptor residue pairs
- Clusters structures with similar H-bond patterns

**Example Settings**:
- eps: 0.4-0.8 (normalized scale)
- min_samples: 2-5

---

### 5️⃣ Combined Weighted Metrics
**When to use**: Comprehensive multi-criteria analysis  
**Best for**: Publication-quality clustering  
**Speed**: ⚡ Slower (computes all metrics)  
**Tip**: Most thorough but takes longer

**What it does**:
- Combines RMSD + Jaccard + CAD + H-bond
- Default weights: RMSD (40%), Jaccard (30%), CAD (20%), H-bond (10%)
- Normalized weighted sum

**Example Settings**:
- eps: 0.3-0.6 (normalized scale)
- min_samples: 2-5

---

## Workflow Example

### Scenario: Analyzing 100 HDOCK docking poses

1. **Start with RMSD** (fast overview)
   ```
   Distance: RMSD
   eps: 12.0
   min_samples: 2
   ```
   Result: Get overall structural clusters

2. **Refine with Jaccard** (contact specificity)
   ```
   Distance: Jaccard Contact Score
   eps: 0.5
   min_samples: 3
   ```
   Result: Find poses with similar binding interfaces

3. **Check H-bonds** (electrostatic verification)
   ```
   Distance: H-bond Pattern
   eps: 0.6
   min_samples: 2
   ```
   Result: Verify electrostatic complementarity

4. **Final analysis with Combined** (comprehensive)
   ```
   Distance: Combined Weighted Metrics
   eps: 0.4
   min_samples: 3
   ```
   Result: Best overall clustering

---

## Tips & Tricks

### 🎯 Choosing eps (distance threshold)

| Metric | Typical Range | Notes |
|--------|--------------|-------|
| RMSD | 5-20 Å | Depends on structure size |
| Jaccard | 0.3-0.7 | Normalized [0,1] |
| CAD | 0.2-0.5 | Normalized [0,1] |
| H-bond | 0.4-0.8 | Normalized [0,1] |
| Combined | 0.3-0.6 | Normalized [0,1] |

### 🎯 Choosing min_samples

- **2**: Maximum sensitivity (many small clusters)
- **3-5**: Balanced (moderate cluster sizes)
- **>5**: Conservative (only strong clusters)

### 🎯 Comparing Results

Run clustering with different metrics and compare:
- **Silhouette Score**: Higher is better (0.5+ is good)
- **Number of Clusters**: Too many? Increase eps or min_samples
- **Noise Points**: Too many? Decrease eps

---

## Common Issues & Solutions

### ❌ "ImportError: mdtraj is required"
**Solution**: Dependencies already installed! This shouldn't happen. Try restarting the app.

### ❌ "No PDB files found"
**Solution**: 
- Check input directory path
- Verify filter pattern (default: `model`)
- Make sure files end in `.pdb`

### ❌ "All structures are noise"
**Solution**:
- Decrease eps value
- Decrease min_samples to 2
- Try different distance metric

### ❌ Clustering takes too long
**Solution**:
- Use RMSD or CAD (fastest metrics)
- Reduce number of input structures
- Avoid Combined metric for >200 structures

---

## Example Results Interpretation

### Scenario: DNA binding protein analysis

**Results with different metrics**:

- **RMSD clustering**: 5 clusters
  - Groups by overall protein position
  - Cluster 0: Protein on major groove
  - Cluster 1: Protein on minor groove

- **Jaccard clustering**: 8 clusters  
  - Groups by contact residues
  - More specific than RMSD
  - Separates different binding residues

- **H-bond clustering**: 6 clusters
  - Groups by H-bond partners
  - Identifies specific base interactions
  - Good for mechanism studies

- **Combined clustering**: 4 clusters
  - Most conservative
  - Combines all criteria
  - Best for final selection

---

## Advanced Usage (CLI)

You can also use the new metrics via command line:

```python
from clustaldemo import PDBClusterer

# Load structures
clust = PDBClusterer(pdb_files, selection='nucleic and name P')
clust.load_structures()

# Option 1: Jaccard
jaccard_matrix = clust.compute_jaccard_contact_matrix()
clust.distance_matrix = jaccard_matrix
clust.cluster_dbscan(eps=0.5, min_samples=2)

# Option 2: Combined with custom weights
weights = {'rmsd': 0.6, 'jaccard': 0.4}
combined_matrix = clust.compute_combined_distance_matrix(weights=weights)
clust.distance_matrix = combined_matrix
clust.cluster_dbscan(eps=0.4, min_samples=3)

# Generate all plots
clust.plot_distance_matrix('distance.png')
clust.plot_dendrogram('dendrogram.png')
clust.save_cluster_separate_pdbs('output/')
```

---

## Next Phase Preview

**Coming in Phase 2**:
- 📊 Per-structure contact statistics
- 🎚️ Interactive weight sliders for combined metrics
- 📈 Metric comparison visualizations
- 🔬 Per-residue contribution analysis
- ⚡ PropKa integration for charge analysis

---

## Need Help?

Check these files:
- `ROADMAP.md` - Future development plans
- `PHASE1_SUMMARY.md` - Technical implementation details
- `README_APP.md` - General app documentation
- `USAGE.md` - Command-line usage guide

---

**Enjoy exploring your docking results with advanced metrics! 🧬**
