# Understanding "All Yellow" Jaccard Distance Matrix

## The Question

**Why is my Jaccard distance matrix all yellow (high distance ~1.0) with no clustering?**

## Short Answer

**This might be CORRECT behavior!** 

If your HDOCK poses represent **different binding modes** with **different interface contacts**, Jaccard distance SHOULD be high (yellow).

## Detailed Explanation

### What Jaccard Measures

Jaccard contact score measures **contact pattern similarity**:

```
Jaccard Similarity = (Shared Contacts) / (All Unique Contacts)

Example:
Pose A contacts: ARG123-DG5, LYS45-DA3, SER67-DG5
Pose B contacts: ARG123-DG5, THR89-DC2, ILE34-DA4

Shared: 1 (ARG123-DG5)
Total unique: 5
Jaccard similarity = 1/5 = 0.2
Jaccard distance = 1 - 0.2 = 0.8 (YELLOW)
```

### What "All Yellow" Means

**All yellow = High Jaccard distance = Low contact overlap**

This happens when:
1. ✅ **Different binding modes** - HDOCK explored different orientations
2. ✅ **Different interface regions** - Poses contact different parts of protein/nucleic acid
3. ✅ **Low scoring poses** - Varied positions, not clustered around optimum

This is **EXPECTED** for:
- Unfiltered HDOCK output (top 100-1000 poses)
- Protein-nucleic systems with multiple binding sites
- Flexible ligands with many possible orientations

### When to Use Jaccard vs RMSD

| Metric | Best For | What it Measures |
|--------|----------|------------------|
| **RMSD** | Docking pose clustering | Overall 3D structural similarity |
| **Jaccard** | Interface analysis | Contact pattern similarity |

**For HDOCK clustering**: **RMSD is usually better**

### Expected Behavior

#### HDOCK Top 100 Poses (Typical)
```
RMSD matrix:
- Some clustering (similar orientations)
- Blue blocks = pose families
- Yellow = different orientations
→ Good for grouping similar poses

Jaccard matrix:
- Mostly yellow (different contacts)
- Few blue spots (rarely identical contacts)
- Expected: Distance 0.7-1.0 for most pairs
→ Shows poses explore different interfaces
```

#### HDOCK Top 10 Poses (High scoring)
```
RMSD matrix:
- More clustering
- Similar poses score well
→ Identifies best binding mode

Jaccard matrix:
- More blue if top poses agree on interface
- Still yellow if multiple binding modes
→ Confirms if top poses use same interface
```

### Diagnostic Questions

**Q: Are contacts being detected?**
Check contact summary CSV:
- If all structures have 30-60 contacts → Working ✓
- If all have 0-5 contacts → Problem ✗

**Q: What distance range do you see?**
- 0.8-1.0 (yellow) = Different contacts → Use RMSD instead
- 0.0-0.5 (purple-blue) = Similar contacts → Jaccard useful
- All exactly 1.0 = Bug (no overlap at all)

**Q: What does RMSD show?**
- RMSD shows clustering → Poses structurally similar but contact different residues
- RMSD all yellow too → Poses completely different (expected for low-ranked HDOCK)

### What to Check

1. **Run debug script:**
```bash
python debug_jaccard.py your_directory/
```

Look for:
- Contact counts per structure (should be 20-80 for protein-DNA)
- Shared contacts between poses (might be 0-5 for different binding modes)
- Jaccard distance values (0.8-1.0 is normal for diverse poses)

2. **Check contact summary CSV:**
- Do structures have contacts?
- Are contact residues different across structures?
- This tells you if poses really are different

3. **Compare with RMSD:**
- If RMSD clusters but Jaccard doesn't → Use RMSD
- If both are yellow → Poses are truly diverse

### Recommendation

**For HDOCK docking pose clustering:**

1. **Primary method: RMSD**
   - Clusters similar 3D orientations
   - Proven for docking analysis
   - Referenced in papers

2. **Secondary analysis: Jaccard**
   - After RMSD clustering, check if clusters have consistent interfaces
   - Identify if different RMSD clusters contact same residues
   - Useful for interface analysis, not primary clustering

3. **Workflow:**
```python
# Step 1: Cluster by RMSD
clust.compute_distance_matrix()
clust.cluster(eps=5.0, min_samples=2, metric='precomputed')

# Step 2: Analyze interface contacts
clust.compute_jaccard_contact_matrix()
clust.export_contacts()  # Check which residues each cluster contacts
```

### When Jaccard IS Useful

- **Protein-protein docking**: Different poses might have similar RMSD but different interfaces
- **After RMSD filtering**: Cluster RMSD-similar poses by interface similarity
- **Interface hot-spot identification**: Find consensus contacts across pose families
- **Mutation design**: Identify residues consistently contacted

### When Jaccard Is NOT Useful

- **Primary clustering of diverse HDOCK output** ← YOUR CASE
- **Very different binding modes** (no shared contacts)
- **Low-scoring poses** (random-ish orientations)

## Summary

**"All yellow" Jaccard matrix is likely CORRECT if:**
- ✅ HDOCK poses have different binding orientations
- ✅ Contact summary shows varied residue contacts
- ✅ You're clustering unfiltered top N poses

**Action:**
- Use **RMSD for clustering**
- Use **Jaccard for interface analysis** after clustering
- Export contacts CSV to see which residues are involved

## New Exports Available

After running Jaccard, you now get:
1. `contact_residues.csv` - All contacts, all structures
2. `contact_residues_summary.csv` - Contact counts and residue lists per structure
3. `distance_matrix.csv` - Pairwise distance values (the actual numbers)

Use these to verify if poses really are different!
