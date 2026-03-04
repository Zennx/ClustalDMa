# Troubleshooting: Structure Type Filter & Stability Plot

**Date:** 2025-11-30  
**Issues:** 
1. Structure Type radio buttons not visible in 3D Viewer
2. Cluster Stability plot showing all zeros

---

## Issue 1: Structure Type Filter Not Showing

### Expected Behavior
In the 3D Viewer sidebar, after the "Filter Cluster" dropdown, there should be radio buttons labeled:
- "Structure Type:" 
  - All
  - Representatives (medoids)
  - Others

### Actual Behavior
Radio buttons not visible in UI

### Root Cause
**The app needs to be restarted to load the new UI elements.**

### Solution
1. **Stop the running Shiny app** (Ctrl+C in terminal)
2. **Restart the app:**
   ```bash
   python app_main.py
   ```
3. **Refresh your browser** (hard refresh: Cmd+Shift+R on Mac, Ctrl+Shift+R on Windows)

### Verification
After restart, check the 3D Viewer sidebar for the "Structure Type:" radio buttons between "Filter Cluster:" and "Structure:" dropdowns.

---

## Issue 2: Cluster Stability Showing All Zeros

### Expected Behavior
Cluster stability plot should show varying bar heights representing stability scores for each cluster.

### Actual Behavior
All bars show 0.000 with message:
```
Fallback Level 1 (Hybrid): HDBSCAN persistence invalid (infinite λ values)
Using: (λ_range × avg_size) / median_intra_cluster_distance
Higher = more stable across densities + more cohesive
```

### Root Cause
The hybrid calculation was looking for `child == cluster_id` in the condensed tree, but HDBSCAN's condensed tree structure uses different node IDs. All lambda_range_scores were zero because no matching nodes were found.

### Fix Applied
Updated the logic to search for both parent and child relationships:
```python
cluster_tree_nodes = tree_df_finite[
    (tree_df_finite['child'] == cid) | (tree_df_finite['parent'] == cid)
]
```

Also added fallback: if no tree nodes found, use cluster size as lambda score.

### Added Debugging
The function now prints diagnostic information to help identify issues:
```
[Stability Fallback] Computing hybrid metric for 8 clusters
[Stability Fallback] Tree has 90 nodes with finite lambda
  Cluster 0: median_dist=0.123, lambda_range_score=15.234
  Cluster 1: median_dist=0.156, lambda_range_score=42.567
  ...
[Stability Fallback] Hybrid scores: [123.821 272.891 ...]
```

### Verification After Restart
1. **Restart the app** (to load updated code)
2. **Run analysis** with HDBSCAN clustering
3. **Check terminal output** for stability fallback diagnostic messages
4. **View Cluster Stability plot** - bars should now show varying heights
5. If still all zeros, check terminal for error messages

---

## Quick Restart Checklist

When you modify Python code in `app_main.py` or `visualization/interactive.py`:

- [ ] Stop running app (Ctrl+C)
- [ ] Run syntax check: `python -m py_compile app_main.py visualization/interactive.py`
- [ ] Restart app: `python app_main.py`
- [ ] Hard refresh browser (Cmd+Shift+R / Ctrl+Shift+R)
- [ ] Test modified feature

---

## Where to Find the Features

### Structure Type Filter
**Location:** 3D Viewer tab → Left sidebar → "Display Settings"  
**Position:** Between "Filter Cluster:" dropdown and "Structure:" dropdown  
**Options:**
- **All** - Show all structures in selected cluster
- **Representatives (medoids)** - Show only cluster representatives (medoid of each cluster)
- **Others** - Show only non-representative structures

**Use Case:** 
- Select "Cluster 0" + "Representatives" → See only the medoid of Cluster 0
- Select "All Clusters" + "Representatives" → See all cluster medoids
- Select "Cluster 1" + "Others" → See all Cluster 1 members except the medoid

### Cluster Stability Plot
**Location:** Cluster Quality tab → "Cluster Stability (Hybrid Stability)" plot  
**What it shows:** 
- X-axis: Clusters sorted by stability (most stable first)
- Y-axis: Stability score (higher = more stable)
- Color: log₁₀(cluster size)
- Hover: Exact stability value and cluster size

**Interpretation:**
- **High scores (>10):** Very stable, cohesive clusters
- **Medium scores (1-10):** Moderately stable
- **Low scores (<1):** Potentially unstable, may be artifacts
- **All zeros:** Issue with calculation (check terminal for errors)

---

## Expected Terminal Output (After Fix)

When viewing Cluster Stability plot with fallback level 1:

```
[Stability Fallback] Computing hybrid metric for 8 clusters
[Stability Fallback] Tree has 90 nodes with finite lambda
  Cluster 7: median_dist=0.234, lambda_range_score=12.456
  Cluster 6: median_dist=0.189, lambda_range_score=18.234
  Cluster 5: median_dist=0.145, lambda_range_score=25.678
  Cluster 4: median_dist=0.198, lambda_range_score=14.567
  Cluster 3: median_dist=0.212, lambda_range_score=16.234
  Cluster 2: median_dist=0.167, lambda_range_score=19.456
  Cluster 1: median_dist=0.156, lambda_range_score=42.567
  Cluster 0: median_dist=0.178, lambda_range_score=15.234
[Stability Fallback] Hybrid scores: [ 53.25 272.86 116.49  73.57  76.56  85.68 177.02  68.54]
```

If you see all zeros in the output, there's still an issue with the calculation logic.

---

## Files Modified

1. **`visualization/interactive.py`** (lines 1076-1135)
   - Enhanced hybrid fallback logic to search both parent and child nodes
   - Added diagnostic printing
   - Added validation to check if lambda scores are all zero
   - Falls back to distance-only if lambda scores invalid

2. **`app_main.py`** (lines 835-837)
   - Structure type radio buttons already present (just needs restart to appear)

---

## Next Steps if Issues Persist

### If Structure Type Still Not Showing
1. Check browser console for JavaScript errors (F12 → Console tab)
2. Verify the input is registered: In browser console, type `Shiny.inputBindings` and check if `structure_type` appears
3. Check if there's a CSS issue hiding the element (inspect element in browser)

### If Stability Still All Zeros
1. **Check terminal output** - Look for the diagnostic messages
2. **Check if lambda scores are being found:**
   - If all `lambda_range_score=0.000`, the tree structure isn't matching
   - Share the terminal output for further debugging
3. **Check condensed tree structure:**
   - In terminal, after clustering, the diagnostic should show tree nodes
   - If tree shows "0 nodes with finite lambda", that's the issue

### Fallback Cascade
If hybrid fails, the system should automatically try:
1. **Fallback Level 2:** Distance-only (1 / median_distance)
   - Should show non-zero values if distance matrix is valid
2. **Fallback Level 3:** Size-only (cluster size)
   - Should always work as last resort

If you're seeing all zeros even at Level 2 or 3, there's a deeper issue with the data.
