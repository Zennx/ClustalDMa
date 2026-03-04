# Cluster Stability Metric - Three-Tier Fallback System

**Date:** 2025-11-30  
**Status:** ✅ Implemented

## Overview

The cluster stability visualization now uses a robust three-tier fallback system to handle cases where HDBSCAN's native persistence score is invalid (due to infinite lambda values from duplicate/near-duplicate structures).

## Fallback Hierarchy

### Level 0: Primary Metric (HDBSCAN Persistence)
- **Metric:** `cluster_persistence_` from HDBSCAN
- **Formula:** λ-integral (excess of mass stability)
- **When used:** When persistence scores are valid (not all 1.0, no infinite values)
- **Interpretation:** Higher = more persistent across density thresholds
- **UI Message:** `✓ Primary Metric: HDBSCAN cluster_persistence_ (λ-integral)`

### Level 1: Hybrid Approach
- **Metric:** `(λ_range × avg_size) / median_intra_cluster_distance`
- **When used:** When persistence is invalid BUT condensed tree has finite lambda values AND distance matrix is available
- **Components:**
  - `λ_range`: max(λ) - min(λ) for cluster (from condensed tree, filtering infinite values)
  - `avg_size`: Average cluster size during lifetime
  - `median_intra_cluster_distance`: Median pairwise distance within cluster
- **Interpretation:** 
  - Higher λ_range × size = more stable across densities
  - Lower median distance = more cohesive
  - Combined: balances density-based and distance-based stability
- **UI Message:** 
  ```
  ℹ️ Fallback Level 1 (Hybrid): HDBSCAN persistence invalid (infinite λ values).
  Using: (λ_range × avg_size) / median_intra_cluster_distance
  Higher = more stable across densities + more cohesive
  ```

### Level 2: Distance-Based Only
- **Metric:** `1 / median_intra_cluster_distance`
- **When used:** When lambda values are unavailable/invalid BUT distance matrix is available
- **Interpretation:** Higher = more cohesive (lower internal distances)
- **UI Message:**
  ```
  ⚠️ Fallback Level 2 (Distance-Only): HDBSCAN persistence and λ-values unavailable.
  Using: 1 / median_intra_cluster_distance
  Higher = more cohesive cluster (lower internal distances)
  ```

### Level 3: Size-Based Proxy (Last Resort)
- **Metric:** Cluster size
- **When used:** When no distance matrix is available
- **Interpretation:** Larger clusters assumed more stable (weak assumption)
- **UI Message:**
  ```
  ⚠️ Fallback Level 3 (Size-Only): No distance matrix available.
  Using cluster size as stability proxy
  ```

## Implementation Details

### Function Signature
```python
def create_cluster_stability_plot(clusterer, labels, distance_matrix=None):
    """
    Parameters:
    -----------
    clusterer : hdbscan.HDBSCAN
        Fitted HDBSCAN clusterer
    labels : array-like
        Cluster labels
    distance_matrix : np.ndarray, optional
        Distance matrix for fallback calculations (required for levels 1-2)
    """
```

### Median Intra-Cluster Distance Calculation
For each cluster:
```python
cluster_indices = np.where(labels == cluster_id)[0]
cluster_dists = []
for i in range(len(cluster_indices)):
    for j in range(i + 1, len(cluster_indices)):
        cluster_dists.append(distance_matrix[cluster_indices[i], cluster_indices[j]])
median_dist = np.median(cluster_dists)
```

### Lambda Range Calculation
From HDBSCAN condensed tree:
```python
tree_df = clusterer.condensed_tree_.to_pandas()
tree_df_finite = tree_df[np.isfinite(tree_df['lambda_val'])]  # Filter infinite values
cluster_nodes = tree_df_finite[tree_df_finite['child'] == cluster_id]
lambda_range = cluster_nodes['lambda_val'].max() - cluster_nodes['lambda_val'].min()
avg_size = cluster_nodes['child_size'].mean()
```

## Why This Approach?

### Problem
HDBSCAN's `cluster_persistence_` becomes invalid when:
- Distance matrix contains zeros (duplicate structures)
- Lambda (λ = 1/distance) becomes infinite
- All persistence scores default to 1.0

### Solution Benefits
1. **Graceful degradation:** Tries best metric first, falls back intelligently
2. **Combines approaches:** Hybrid metric balances density-based (HDBSCAN philosophy) with distance-based (interpretable)
3. **User transparency:** Clear warnings explain which metric is being used and why
4. **Handles edge cases:** Works even without distance matrix (though less informative)

## User Interpretation Guide

### Reading the Plot
- **Y-axis:** Stability metric value (higher = more stable)
- **X-axis:** Clusters sorted by stability (descending)
- **Color:** log₁₀(cluster size) - darker = larger clusters
- **Hover:** Shows exact metric value and cluster size

### What Makes a "Good" Cluster?
- **High stability score** (regardless of metric)
- **Clear separation** from low-stability clusters (gap in bar chart)
- **Reasonable size** (not too small, but size alone doesn't guarantee quality)

### When to Trust Each Metric
- **Primary (λ-integral):** Most reliable - directly from HDBSCAN theory
- **Hybrid:** Very good - combines stability across densities with cohesiveness
- **Distance-only:** Good for cohesiveness, but ignores density hierarchy
- **Size-only:** Weakest - large clusters aren't always stable

## Testing Scenarios

### Test Case 1: Clean Dataset (No Duplicates)
- **Expected:** Level 0 (Primary metric)
- **Result:** Persistence scores valid, all clusters ranked

### Test Case 2: Dataset with Duplicates
- **Expected:** Level 1 (Hybrid) if enough finite lambda values
- **Result:** Filters infinite lambdas, computes hybrid stability

### Test Case 3: All Infinite Lambdas (Extreme Duplicates)
- **Expected:** Level 2 (Distance-only)
- **Result:** Falls back to median intra-cluster distances

### Test Case 4: No Distance Matrix Available
- **Expected:** Level 3 (Size-only)
- **Result:** Uses cluster sizes with warning

## Files Modified

1. **`visualization/interactive.py`**
   - Enhanced `create_cluster_stability_plot()` with three-tier fallback
   - Added `distance_matrix` parameter
   - Implemented median intra-cluster distance calculation
   - Added detailed UI warning messages

2. **`app_main.py`**
   - Updated call to pass `clust.distance_matrix`
   - Added safety check for distance matrix availability

## Future Enhancements

- [ ] Add toggle to view all metrics side-by-side
- [ ] Show scatter plot: stability vs size (identify outliers)
- [ ] Add cluster quality score combining stability + silhouette + size
- [ ] Export stability metrics to CSV
- [ ] Add stability threshold filter for cluster selection

## References

- HDBSCAN Paper: Campello et al. (2013) - "Density-Based Clustering Based on Hierarchical Density Estimates"
- Excess of Mass: λ-integral measures stability across density thresholds
- Cohesiveness: Lower intra-cluster distances indicate tighter grouping
