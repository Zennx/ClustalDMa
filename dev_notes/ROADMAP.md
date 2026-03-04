# ClustalDM Development Roadmap

## 🎯 Current Status: v1.0 - Core RMSD Clustering

### ✅ Implemented Features
- RMSD-based structural clustering with DBSCAN
- Interactive web interface with Shiny
- Multiple visualization methods (MDS, PCA, t-SNE)
- Hierarchical clustering dendrogram
- PyMOL script generation
- Representative structure extraction
- Cluster organization and export

---

## 🚀 Planned Features

### Phase 1: Additional Distance Metrics (High Priority)

#### 1.1 Jaccard Contact Score
**Purpose**: Measure similarity based on interface contacts

**Implementation**:
```python
def compute_jaccard_contact_matrix(pdb_files, distance_cutoff=4.5):
    """
    Compute Jaccard similarity for protein-nucleic acid contacts
    
    Parameters:
    - distance_cutoff: Å threshold for contact definition
    
    Returns:
    - Jaccard distance matrix (1 - Jaccard similarity)
    """
    # For each structure:
    # 1. Identify protein-nucleic acid contacts
    # 2. Create contact set (residue pairs < cutoff)
    # 3. Compute Jaccard similarity between all pairs
    # 4. Convert to distance: 1 - similarity
```

**Dependencies**: 
- `mdtraj` for contact calculation
- Custom contact definition logic

**Use Case**: Identify poses with similar binding interfaces regardless of backbone RMSD

---

#### 1.2 Contact Area Distance (CAD)
**Purpose**: Quantify differences in binding interface area

**Implementation**:
```python
def compute_cad_matrix(pdb_files):
    """
    Compute Contact Area Distance between structures
    
    Steps:
    1. Calculate interface SASA for each structure
    2. Identify contacting residue pairs
    3. Compute area contribution per contact
    4. Calculate pairwise CAD scores
    """
    # CAD = sum of differences in contact areas
```

**Dependencies**:
- `freesasa` for SASA calculation
- Contact identification algorithm

**Use Case**: Find poses with similar contact surface areas

---

#### 1.3 Hydrogen Bond Distance
**Purpose**: Cluster by H-bond pattern similarity

**Implementation**:
```python
def compute_hbond_distance_matrix(pdb_files):
    """
    Compute distance based on H-bond patterns
    
    Steps:
    1. Identify H-bonds using MDTraj
    2. Create H-bond fingerprint (donor-acceptor pairs)
    3. Compare fingerprints between structures
    4. Distance = fraction of non-matching H-bonds
    """
```

**Dependencies**:
- `mdtraj.baker_hubbard()` for H-bond detection

**Use Case**: Group poses with similar hydrogen bonding networks

---

#### 1.4 Combined Weighted Metrics
**Purpose**: Multi-criteria clustering with user-defined weights

**Implementation**:
```python
def compute_combined_distance_matrix(pdb_files, weights):
    """
    Combine multiple distance metrics
    
    Parameters:
    weights: dict like {'rmsd': 0.4, 'jaccard': 0.3, 'cad': 0.2, 'hbond': 0.1}
    
    Returns:
    Weighted combined distance matrix
    """
    # Normalize each metric to [0, 1]
    # Apply weights
    # Sum to create combined distance
```

**UI Addition**:
- Sliders for weight adjustment
- Real-time weight normalization
- Save/load weight presets

---

### Phase 2: Extended Structural Metrics

#### 2.1 deltaSASA (ΔSASA)
**Purpose**: Measure burial of surface area upon binding

```python
def compute_delta_sasa(protein_pdb, ligand_pdb, complex_pdb):
    """
    ΔSASA = (SASA_protein + SASA_ligand) - SASA_complex
    
    Indicates binding interface size
    """
```

**Display**: 
- Per-structure ΔSASA values in cluster tables
- ΔSASA distribution plots
- Correlation with RMSD/clustering

---

#### 2.2 Contact Number Analysis
**Purpose**: Count and categorize interface contacts

```python
class ContactAnalyzer:
    def count_contacts(self, pdb_file, cutoff=4.5):
        """Count different contact types"""
        return {
            'total_contacts': int,
            'polar_contacts': int,
            'nonpolar_contacts': int,
            'charged_contacts': int,
            'aromatic_contacts': int,
            'h_bonds': int
        }
```

**Display**:
- Contact count distribution per cluster
- Contact type breakdown pie charts
- Identify "hot spot" residues

---

#### 2.3 Charge Ratio and Electrostatics
**Purpose**: Analyze electrostatic complementarity

```python
def analyze_electrostatics(pdb_file):
    """
    Using PropKa for pKa prediction and charge analysis
    
    Returns:
    - Charge distribution at interface
    - Favorable vs unfavorable charge pairs
    - Electrostatic energy estimates
    """
```

**Dependencies**: `propka`

**Display**:
- Charge complementarity scores
- Salt bridge identification
- Favorable interaction counts

---

### Phase 3: Advanced Analysis Tools

#### 3.1 Per-Residue Contribution Analysis
```python
def compute_residue_contributions(cluster_structures):
    """
    Identify which residues contribute most to:
    - Binding energy (if scores available)
    - Interface area
    - H-bond formation
    - Cluster differentiation
    """
```

**Visualization**:
- Heatmap of residue importance
- Export for PyMOL coloring scripts
- Residue conservation across clusters

---

#### 3.2 Cluster Consensus Structure
```python
def generate_consensus_structure(cluster_members):
    """
    Create representative structure showing:
    - Average ligand position
    - Conserved contacts
    - Variable regions
    """
```

**Output**:
- Consensus PDB file per cluster
- Confidence/variability annotations
- PyMOL scripts to visualize uncertainty

---

#### 3.3 Dynamic Features (Future)
**If trajectories available:**
```python
def analyze_dynamics(trajectory_files):
    """
    - RMSD fluctuations
    - Contact persistence
    - Conformational sampling
    """
```

---

### Phase 4: Enhanced User Interface

#### 4.1 Interactive Parameter Tuning
- [ ] Real-time clustering updates
- [ ] Parameter sensitivity analysis
- [ ] Automatic optimal eps suggestion
- [ ] Cluster stability visualization

#### 4.2 Advanced Visualizations
- [ ] 3D scatter plots (plotly)
- [ ] Interactive heatmaps
- [ ] Contact network graphs
- [ ] Time-series if multiple runs

#### 4.3 Comparison Mode
- [ ] Compare clustering from different metrics
- [ ] Sankey diagrams showing cluster membership
- [ ] Metric correlation analysis

---

## 📊 Implementation Priority

### High Priority (Next Sprint)
1. ✅ Fix table display issue (Model/Index/RMSD columns)
2. ✅ Center PyMOL download button
3. ✅ Add export directory browser
4. 🔄 Implement Jaccard contact score
5. 🔄 Add contact counting (basic)
6. 🔄 Integrate FreeSASA for ΔSASA

### Medium Priority
7. CAD (Contact Area Distance) implementation
8. H-bond pattern analysis
9. PropKa integration for charge analysis
10. Combined weighted metrics

### Lower Priority (Advanced)
11. Per-residue contribution analysis
12. Consensus structure generation
13. 3D interactive visualizations
14. Batch processing multiple HDOCK runs

---

## 🛠️ Technical Requirements

### New Dependencies
```python
# requirements_extended.txt
mdtraj>=1.9.7        # Contact analysis, H-bonds
freesasa>=2.1.0      # SASA calculation
propka>=3.5.0        # pKa prediction, charge analysis
plotly>=5.0.0        # Interactive 3D plots (optional)
networkx>=2.8.0      # Contact network analysis (optional)
```

### Data Structure Extensions
```python
class ClusterMetrics:
    """Store extended metrics per structure"""
    rmsd_matrix: np.ndarray
    jaccard_matrix: np.ndarray  # Coming
    cad_matrix: np.ndarray      # Coming
    hbond_matrix: np.ndarray    # Coming
    
    # Per-structure metrics
    contact_counts: List[Dict]
    delta_sasa: List[float]
    charge_ratios: List[Dict]
    hbond_lists: List[List[Tuple]]
```

---

## 🎯 Success Criteria

### Phase 1 Complete When:
- [x] RMSD clustering works reliably
- [ ] Jaccard contact metric implemented
- [ ] Users can choose distance metric in UI
- [ ] Basic contact statistics displayed
- [ ] FreeSASA integration complete

### Phase 2 Complete When:
- [ ] All individual metrics calculate correctly
- [ ] Combined weighted metrics functional
- [ ] Comprehensive metric comparison available
- [ ] Export includes all computed metrics

### Phase 3 Complete When:
- [ ] Advanced analysis tools ready
- [ ] Publication-quality figures generated
- [ ] Batch processing supported
- [ ] Full documentation complete

---

## 📝 Notes

- **Validation**: Each new metric needs validation against known test cases
- **Performance**: Large datasets (>1000 structures) may need optimization
- **Memory**: Consider lazy loading for large distance matrices
- **Compatibility**: Ensure all tools work with various PDB formats from different docking programs

---

## 🤝 Contributing

Areas needing development:
1. Algorithm optimization for large datasets
2. Additional distance metrics
3. Visualization improvements
4. Documentation and tutorials
5. Test cases and benchmarks

---

**Last Updated**: November 10, 2025
**Version**: 1.0 (RMSD baseline)
**Next Milestone**: v1.1 - Multi-metric clustering
