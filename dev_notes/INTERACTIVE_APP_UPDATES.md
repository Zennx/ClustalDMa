# Interactive App Updates

## Date: 2025-11-13

## Issues Fixed

### 1. Empty Interactive Map 
**Problem**: Scatter plot was empty  
**Cause**: Using PCA on distance matrix instead of MDS  
**Solution**: Switched to MDS (Multidimensional Scaling) with `dissimilarity='precomputed'`
- MDS is designed for distance matrices
- PCA is for feature matrices
- Now uses: `MDS(n_components=2, dissimilarity='precomputed')`

### 2. Hotspot Histogram Improvements
**Problem**: Bar chart not smooth, not sorted by sequence  
**Solution**: Created genome-browser style smooth plot
- Extracts residue numbers using regex: `(\d+)` from "ARG123" → 123
- Sorts by residue position along sequence
- Displays as smooth area plot (`fill='tozeroy'`)
- Color: Red gradient (`rgb(220, 53, 69)` with 30% opacity fill)
- Hover shows: residue name, position, frequency %

**New parameters**:
```python
create_hotspot_histogram(hotspots_df, top_n=None, smooth=True)
# top_n=None shows ALL residues
# smooth=True creates area plot instead of bars
```

### 3. 3D Viewer Settings
**Problem**: No customization controls  
**Solution**: Added left sidebar with controls
- Protein style: cartoon, stick, sphere, line
- Protein color: cyan, gray, green, blue, purple
- DNA style: stick, sphere, line, cartoon
- Show surface: checkbox
- Surface opacity: slider (0.1-1.0)
- Zoom level: slider (0.5-2.0x)

**New function signature**:
```python
create_mol_viewer_html(
    pdb_file, 
    width=800, 
    height=600,
    protein_style='cartoon',
    protein_color='cyan',
    dna_style='stick',
    show_surface=False,
    surface_opacity=0.5,
    zoom_level=1.0
)
```

## Changes Made

### visualization/interactive.py
1. Added `import numpy as np` 
2. Replaced PCA with MDS in `create_interactive_scatter()`
3. Completely rewrote `create_hotspot_histogram()`:
   - Extract residue numbers for sorting
   - Sort by sequence position
   - Create smooth area plot (genome browser style)
   - Red color scheme maintained
4. Enhanced `create_mol_viewer_html()` with customization parameters:
   - Protein style options
   - Color options
   - Surface rendering
   - Adjustable zoom

### app_interactive.py
1. Updated hotspot plot to use: `create_hotspot_histogram(hotspots_df, top_n=None, smooth=True)`
2. Restructured 3D Viewer tab with sidebar:
   - Left sidebar: All display controls
   - Right panel: Structure info + viewer
3. Added `structure_info()` output function
4. Updated `structure_viewer()` to use all settings from inputs
5. Pass all style parameters to `create_mol_viewer_html()`

## User Feedback Addressed

✅ **"Interactive map is empty"**  
   - Fixed by using MDS instead of PCA

✅ **"Make histogram smooth (like genome browsers)"**  
   - Now shows smooth area plot sorted by residue position
   
✅ **"Arranged according to sequence number"**  
   - Extracts residue numbers and sorts sequentially
   
✅ **"I like the colours tho :D"**  
   - Kept red color scheme, enhanced with smooth gradient
   
✅ **"Can we have some basic settings/parameters for 3D model"**  
   - Added comprehensive controls in left sidebar:
     * Protein/DNA styles
     * Colors
     * Surface options
     * Zoom control
     
✅ **"We can have it take up the left toolbar space too"**  
   - Implemented sidebar layout for 3D viewer tab

## Testing

Run the app:
```bash
conda activate structural
shiny run app_interactive.py
```

Expected behavior:
1. Interactive Map: Should show MDS projection with points colored by cluster
2. Hotspot Histogram: Smooth red area plot sorted left-to-right by residue position
3. 3D Viewer: Sidebar with style controls, viewer updates on change

## 10 Clusters from Jaccard

User reported getting 10 clusters from Jaccard clustering. This suggests:
- Jaccard is working correctly (detecting binding modes)
- May be more granular than expected (eps=0.3 might be too strict)
- Could adjust DBSCAN parameters: increase `eps` for fewer clusters
- Or increase `min_samples` to merge small clusters

Consider testing with different parameters in future:
```python
clust.cluster_dbscan(eps=0.4, min_samples=3)  # Fewer, larger clusters
clust.cluster_dbscan(eps=0.5, min_samples=2)  # Even fewer clusters
```
