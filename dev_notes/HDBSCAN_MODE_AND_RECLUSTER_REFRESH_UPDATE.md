# HDBSCAN Mode Expansion + CLI Alignment + Recluster Refresh

**Date:** 2026-04-16

## Summary
This update aligned clustering controls across backend, GUI, and CLI while fixing incomplete UI refresh behavior after re-clustering.

## Scope

### 1) Backend (`core/clusterer.py`)
- Expanded `cluster_hdbscan()` to support:
  - `cluster_selection_epsilon`
  - `cluster_selection_method` in `{eom, leaf, mds, tsne, umap, pca}`
- Added `_project_distance_matrix_2d()` to generate 2D embeddings from precomputed distance matrices.
- Added mode routing logic:
  - `eom` / `leaf`: HDBSCAN on precomputed distances (`metric='precomputed'`)
  - `mds` / `tsne` / `umap` / `pca`: project to 2D first, then HDBSCAN on embedding (`metric='euclidean'`) with EOM extraction
- Added metadata for downstream inspection:
  - `self.clustering_input_mode`
  - `self.clustering_projection_2d`

### 2) GUI Controls (`app_main.py`)
- Added new HDBSCAN controls in sidebar:
  - `hdbscan_epsilon`
  - `hdbscan_selection_method` with options:
    - EOM (precomputed matrix)
    - Leaf (precomputed matrix)
    - MDS/t-SNE/UMAP/PCA (2D embedding + EOM)
- Passed these parameters into both:
  - full run path
  - re-cluster path

### 3) Re-cluster Refresh Behavior (`app_main.py`)
- Problem: after re-cluster, only some label-driven views updated.
- Root cause: many render functions depended on global mutable state (`clust`) without guaranteed reactive invalidation.
- Fix:
  - Added `cluster_refresh_tick` reactive value.
  - Added `bump_cluster_refresh()` helper.
  - Called refresh bump at end of full run and re-cluster.
  - Added refresh dependency (`_ = cluster_refresh_tick.get()`) to cluster-dependent renderers so cards/heatmaps/trees/hotspot views all redraw.
- Re-cluster logic was upgraded to recompute and set:
  - `rmsd_stats`
  - `interface_stats`
  - `hotspots`
  - `cluster_summary` (from canonical `get_cluster_summary(threshold=50)` path)
  - motif/TM/stability enrichment where available
- Re-cluster selector updates now also refresh intra-cluster dropdowns:
  - `jaccard_intra_cluster`
  - `rmsd_intra_cluster`

### 4) CLI Alignment (`clustal_cli.py`)
- Added GUI-aligned chain input:
  - `--protein-chains` (comma/space-separated)
- Added HDBSCAN controls:
  - `--hdbscan-epsilon`
  - `--hdbscan-selection-method` (`eom`, `leaf`, `mds`, `tsne`, `umap`, `pca`)
- Implemented chain-selection helpers in CLI similar to app behavior:
  - parse labels
  - resolve labels to MDTraj chain indices
  - fallback auto-detection for protein chains
- Maintained backward compatibility:
  - deprecated `--ligand-id` fallback path still accepted
- Updated CLI run logs and summary-report parameter cards to show:
  - protein chains
  - HDBSCAN epsilon
  - HDBSCAN method

## Documentation Updated
- `CLI_QUICKSTART.md`
  - replaced primary `--ligand-id` examples with `--protein-chains`
  - documented `--hdbscan-epsilon` and `--hdbscan-selection-method`
  - added mode behavior notes and examples
- `README.md`
  - updated CLI examples and feature summary for new HDBSCAN modes
- `QUICKSTART.md`
  - added GUI guidance for epsilon + selection method
  - noted full refresh behavior after re-cluster
- `README_APP.md`
  - updated app workflow/feature descriptions to match new controls
  - replaced outdated DBSCAN-style options block with CLI quickstart reference

## Validation
- Editor diagnostics checked for:
  - `core/clusterer.py`
  - `app_main.py`
  - `clustal_cli.py`
- No errors reported after changes.

## Notes
- Existing unrelated workspace changes (e.g., DS_Store, prior visualization edits) were intentionally left untouched.
- UMAP embedding mode requires `umap-learn`.
