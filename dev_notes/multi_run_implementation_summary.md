# Multi-Run Comparison Implementation Summary

## Overview
Implemented comprehensive multi-directory support for comparing docking results from different methods, ligands, and modes in a single analysis session.

## Changes Made

### 1. Core IO Utilities (`core/io_utils.py`)
- **New Function**: `parse_folder_metadata(folder_path)`
  - Parses folder names following pattern: `DockingMethod_ProteinModel_LigandType_DockingMode`
  - Example: `HDOCK_7vxr_dsDNA_FreeTF` → extracts all 4 metadata fields
  - Returns dict with keys: docking_method, protein_model, ligand_type, docking_mode, run_id
  - Handles non-conforming folders gracefully (sets fields to "Unknown")

- **Updated Function**: `find_pdb_files(directories, filter_pattern, extract_metadata)`
  - Added `extract_metadata` parameter (default: False for backward compatibility)
  - When True, returns tuple: (pdb_files, metadata_df)
  - metadata_df contains columns: structure, pdb_path, docking_method, protein_model, ligand_type, docking_mode, run_id
  - Structure names include run_id prefix to handle duplicate filenames: `run_id/filename.pdb`

### 2. PDB Clusterer (`core/clusterer.py`)
- **Updated Constructor**: Added `metadata` parameter
  - Stores metadata DataFrame for later use
  - Optional parameter maintains backward compatibility

- **New Method**: `align_all_structures(reference_idx, align_selection)`
  - Aligns all loaded structures to a reference by protein backbone
  - **Critical for multi-run comparisons** to ensure comparable coordinates
  - Uses MDAnalysis.analysis.align.alignto()
  - Default alignment selection: 'protein and name CA'

- **Updated Method**: `get_interface_stats()`
  - Now uses metadata-aware structure names if metadata is available
  - Merges metadata columns into results DataFrame
  - Passes structure names (not file paths) to InterfaceAnalyzer
  - Adds columns: docking_method, protein_model, ligand_type, docking_mode, run_id

### 3. Interface Analyzer (`core/analysis.py`)
- **Updated Method**: `get_interface_stats(structure_names, ...)`
  - Changed parameter from `pdb_files` to `structure_names`
  - Supports both basenames and metadata-derived names (run_id/filename)
  - More flexible for multi-run analysis

### 4. Visualization (`visualization/interactive.py`)
- **Updated Function**: `create_scatter_multimethod(df, distance_matrix, method, group_by)`
  - Added `group_by` parameter (default: 'cluster')
  - Supports grouping by: cluster, ligand_type, docking_mode, docking_method, protein_model
  - Enhanced hover text to show metadata fields
  - Automatic fallback to 'cluster' if specified column doesn't exist
  - Updates plot title to reflect grouping choice
  - Color palette automatically assigned to groups

### 5. GUI Application (`app_main.py`)

#### UI Changes:
- **Sidebar Input**: Changed from single `input_text("pdb_dir")` to `input_text_area("pdb_dirs")`
  - Accepts multiple directories, one per line
  - Added helper text explaining folder naming convention
  - Browse button now appends directories instead of replacing

- **New Control**: `input_checkbox("align_proteins", value=True)`
  - Enables/disables protein alignment
  - Includes helper text emphasizing importance for multi-run analysis

- **Interactive Maps Tab**: Added grouping selector
  - New `input_select("scatter_group_by")` with options:
    - Cluster (default)
    - Ligand Type
    - Docking Mode
    - Docking Method
    - Protein Model
  - Helper text noting metadata grouping only available for multi-run mode

#### Server Logic Changes:
- **Directory Parsing**: 
  - Splits input text by newlines
  - Validates each directory
  - Logs valid/invalid directories

- **Multi-Run Detection**:
  - Automatically detects if multiple directories provided
  - Switches to metadata extraction mode
  - Logs metadata statistics (unique methods, ligands, modes)

- **Alignment Integration**:
  - Calls `clust.align_all_structures()` if `align_proteins` is checked
  - Performed after loading structures, before distance calculations

- **Visualization Update**:
  - Passes `group_by` parameter to `create_scatter_multimethod()`
  - Uses scatter_group_by input value

- **Export Update**:
  - Includes `group_by` parameter when exporting cluster visualization

- **Browse Handler Update**:
  - Appends directories to existing list instead of replacing
  - Better for multi-run workflow

## Key Features

### 1. Metadata Extraction
- Automatic parsing from folder names
- Structured storage in DataFrame
- Available in all result tables and visualizations

### 2. Duplicate Filename Handling
- Structure names prefixed with run_id
- Example: `HDOCK_7vxr_dsDNA_FreeTF/model_1.pdb` vs `AutoDock_7vxr_dsDNA_FreeTF/model_1.pdb`
- Ensures unique identification across runs

### 3. Protein Alignment
- Essential for valid multi-run comparisons
- Aligns all structures to common reference frame
- Makes coordinates comparable across different docking runs

### 4. Flexible Grouping
- Color structures by any metadata field
- Useful for comparing:
  - Different docking methods
  - Different ligand types (dsDNA vs ssDNA)
  - Different docking modes (FreeTF vs Bound)
  - Different protein models

### 5. Backward Compatibility
- Single-directory mode still fully functional
- All new parameters have sensible defaults
- Existing scripts/workflows unaffected

## Usage Example

```python
# In the GUI:
# 1. Enter directories (one per line):
#    /data/HDOCK_7vxr_dsDNA_FreeTF
#    /data/AutoDock_7vxr_dsDNA_FreeTF
#    /data/HDOCK_7vxr_ssDNA_FreeTF

# 2. Enable "Align proteins before analysis"

# 3. Run Analysis

# 4. In Interactive Maps tab:
#    - Select "Color by: Ligand Type" to compare dsDNA vs ssDNA
#    - Or "Color by: Docking Method" to compare HDOCK vs AutoDock
```

## Testing Recommendations

1. **Single Directory**: Verify backward compatibility
2. **Two Directories**: Test basic multi-run with same ligand, different method
3. **Three+ Directories**: Test with varied metadata (different ligands, modes)
4. **Non-Standard Names**: Verify graceful handling of folders not following convention
5. **Alignment**: Compare results with/without alignment to verify impact

## Future Enhancements

Potential improvements for future versions:
- Support for custom metadata via CSV file upload
- More sophisticated alignment options (e.g., align by specific domains)
- Statistical comparison of clusters across runs
- Export metadata-grouped plots automatically
- Support for loading pre-aligned structures (skip alignment step)

## Documentation

Created `MULTI_RUN_USAGE.md` with:
- Folder naming convention
- Step-by-step usage guide
- Example workflows
- Troubleshooting tips
- Important notes on alignment
