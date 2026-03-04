# Cluster Summary Feature

## Date: 2025-11-13

## New Features Added

### 1. Cluster Summary Table
**Purpose**: Show consensus binding residues for each cluster

**New Function**: `get_cluster_summary(threshold=50)`
- Returns DataFrame with columns:
  - `cluster`: Cluster ID
  - `n_structures`: Number of structures in cluster
  - `n_consensus`: Number of consensus residues (>50% frequency)
  - `consensus_residues`: Comma-separated list of consensus residues
  - `top_residues`: Top 10 most frequent residues
  - `binding_mode`: Description (e.g., "Mode 0" or "Noise")

**Location**: 
- `core/analysis.py`: `InterfaceAnalyzer.get_cluster_summary()`
- `core/clusterer.py`: `PDBClusterer.get_cluster_summary()`

### 2. Per-Cluster Hotspot Analysis
**Purpose**: View binding pattern for individual clusters

**Features**:
- Cluster selector dropdown
- Smooth genome-browser style plot for selected cluster
- Shows only residues from structures in that cluster
- Title shows number of structures in cluster

### 3. Explanation Note for Histogram
**Purpose**: Clarify why different residues appear at same position

**Explanation**:
> Different residues may appear at the same position because this shows **all structures combined**. 
> Each docking pose can have different protein residues contacting DNA. To see cluster-specific 
> patterns, check the **Cluster Summary** tab for consensus residues within each binding mode.

## UI Updates

### New Tab: "📊 Cluster Summary"
Contains two sections:
1. **Cluster Summary Table**
   - Shows all clusters with consensus residues
   - Consensus = residues present in >50% of cluster members
   
2. **Per-Cluster Hotspot Analysis**
   - Dropdown to select cluster
   - Genome-browser style plot for that cluster only
   - Shows binding pattern specific to that mode

### Updated Tab: "🔥 Binding Hotspots"
- Added explanatory note (yellow info box)
- Clarifies that this shows ALL structures combined
- Points users to Cluster Summary for mode-specific patterns

## Debug Improvements

### Scatter Plot Debugging
Added console logging to diagnose empty plot issue:
- Print stats DataFrame shape
- Print distance matrix shape
- Print MDS computation progress
- Print coordinate ranges
- Print number of points and clusters
- Print HTML length

### MDS Improvements
- Added `n_init=1` to avoid FutureWarning
- Added `max_iter=300` for stability
- Added debug prints for coordinate ranges
- Better error handling with traceback

## Code Changes Summary

### core/analysis.py
```python
@staticmethod
def get_cluster_summary(contact_residue_pairs, labels, threshold=50):
    """Get summary of all clusters with consensus residues"""
    # For each cluster:
    #   - Count structures
    #   - Find consensus residues (>threshold% frequency)
    #   - Get top hotspots
    # Returns DataFrame
```

### core/clusterer.py
```python
def get_cluster_summary(self, threshold=50):
    """Get summary table of all clusters"""
    return InterfaceAnalyzer.get_cluster_summary(
        self.contact_residue_pairs, self.labels, threshold
    )
```

### visualization/interactive.py
```python
# Added debug prints to create_interactive_scatter()
print(f"Computing MDS projection for {n} structures...")
print(f"✓ MDS complete. Coords range: X=[...], Y=[...]")
print(f"Creating scatter plot with {n} points, clusters: {set(...)}")
print(f"✓ Scatter plot created successfully")
```

### app_interactive.py
1. Added `cluster_summary` reactive value
2. Added cluster selector: `ui.input_select("cluster_select", ...)`
3. Added new tab: "📊 Cluster Summary"
4. Added output functions:
   - `cluster_summary_table()`: Show consensus table
   - `cluster_hotspot_plot()`: Per-cluster pattern
   - `hotspot_note()`: Explanation box
5. Enhanced error handling with traceback

## Usage

### View Cluster Consensus
1. Run clustering
2. Click "📊 Cluster Summary" tab
3. See table with consensus residues for each cluster
4. Consensus residues = present in >50% of structures

### View Cluster-Specific Pattern
1. In "📊 Cluster Summary" tab
2. Select cluster from dropdown
3. View smooth plot showing binding pattern for ONLY that cluster
4. Compare different clusters to see different binding modes

### Understanding Hotspot Histogram
- **Binding Hotspots tab** = ALL structures combined
- **Cluster Summary tab** = Per-cluster patterns
- Different residues at same position = different docking poses
- Use consensus residues to identify common patterns within a mode

## Example Interpretation

If you get 11 clusters with Jaccard:

**Cluster Summary Table** might show:
```
Cluster | Structures | Consensus Residues           | Binding Mode
--------|------------|------------------------------|-------------
0       | 25         | ARG123,LYS45,SER67,THR89    | Mode 0
1       | 25         | ARG123,GLU56,ASP78,LYS90    | Mode 1
2       | 5          | LYS45,ARG88,GLU99           | Mode 2
...
-1      | 5          | (various)                    | Noise
```

**Interpretation**:
- Mode 0 and Mode 1 both use ARG123 (common anchor)
- But Mode 0 uses LYS45+SER67, Mode 1 uses GLU56+ASP78
- This suggests two different binding orientations sharing some contacts
- Mode 2 is distinct (different residues)
- Noise cluster (-1) has no consistent pattern

## Testing

Run the app and test:
```bash
conda activate structural
shiny run app_interactive.py
```

1. **Check scatter plot**: Should now show points (debug prints in console)
2. **Check Cluster Summary**: Should show consensus residues
3. **Check per-cluster plots**: Select different clusters, compare patterns
4. **Check explanation note**: Yellow info box should explain why residues vary

## Next Steps

Potential enhancements:
- Adjustable consensus threshold (slider for 30-70%)
- Export cluster summary to CSV
- Compare two clusters side-by-side
- Highlight consensus residues in 3D viewer
- Color code residues by frequency in structure viewer
