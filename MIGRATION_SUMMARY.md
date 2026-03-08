# ClustalDM - AlphaFold2 Migration Summary

## Overview
This document summarizes the key changes made to convert ClustalDM from an HDOCK-focused docking analysis tool to an AlphaFold2 model clustering application.

## Major Changes

### 1. File Input System (COMPLETED)

**Changed:**
- Removed HDOCK .out file parsing
- Removed receptor/ligand PDB upload requirements
- Removed HDOCKParser import and dependencies

**Added:**
- Directory-based input for AlphaFold2 models
- Support for both PDB and mmCIF file formats
- Optional reference model upload for RMSD superposition

**Files Modified:**
- `core/io_utils.py`: Updated `find_pdb_files()` to search for .pdb and .cif files
- `app_main.py`: Replaced file upload UI with directory text input and reference model upload

### 2. RMSD Calculation (COMPLETED)

**Changed:**
- Default RMSD selection changed from `nucleic and name P` to `protein and name CA`
- RMSD now uses alignment by default (suitable for AlphaFold2 models)
- Reference index can be specified for superposition

**Added:**
- Optional reference model parameter
- Reference model uploaded can be used for consistent superposition
- Support for standard protein structure alignment

**Files Modified:**
- `app_main.py`: Updated RMSD calculation parameters in analysis workflow
- UI: Added reference model upload field

### 3. SASA Removal (COMPLETED)

**Removed:**
- All ΔSASA (buried surface area) calculations
- freesasa dependency
- SASA worker function for parallel computation
- SASA-related UI elements (sasa_mode radio buttons, max_poses slider)
- SASA columns from data exports
- SASA metrics from interface quality cards
- SASA statistics from cluster summary tables

**Replaced:**
- Interface quality cards now show contact-based metrics instead of SASA
- Focus shifted to protein/nucleic acid contact counts

**Files Modified:**
- `app_main.py`: Removed SASA computation, worker function, UI elements, and display logic
- `requirements.txt`: Removed freesasa>=2.1.0 dependency

### 4. UI Updates (COMPLETED)

**Changed:**
- App title: "ClustalDM - AlphaFold2 Model Clustering"
- Upload section header: "Upload AlphaFold2 Models"
- Removed pose extraction controls (max_poses slider)
- Removed SASA mode selection
- Removed duplicate threshold slider
- Updated help text throughout to reflect AlphaFold2 focus

**Added:**
- Directory text input for model location
- Reference model file upload
- Clearer separation of required vs optional inputs

**Files Modified:**
- `app_main.py`: Updated all UI elements and help text

### 5. Documentation (COMPLETED)

**Updated:**
- `README_APP.md`: Rewrote to focus on AlphaFold2 workflows
- `requirements.txt`: Removed HDOCK-specific dependencies
- `QUICKSTART.md`: Created new quick start guide for AlphaFold2 usage

**Key Documentation Changes:**
- Removed all HDOCK references
- Added AlphaFold2-specific usage examples
- Updated feature list to reflect current capabilities
- Emphasized Jaccard contact similarity and pLDDT support

## Analysis Workflow Changes

### Old Workflow (HDOCK)
1. Upload .out files + receptor + ligand PDBs
2. Extract poses using createpl
3. Compute Jaccard contact matrix
4. Compute RMSD (nucleic P atoms, no alignment)
5. Compute ΔSASA per structure
6. Cluster by Jaccard
7. Display SASA-colored visualizations

### New Workflow (AlphaFold2)
1. Specify directory containing AF2 models
2. (Optional) Upload reference model
3. Load PDB/mmCIF files from directory
4. Compute Jaccard contact matrix
5. Compute RMSD (protein CA atoms, with alignment)
6. Cluster by Jaccard
7. Display contact-based visualizations
8. Support pLDDT visualization in Molstar viewer

## Removed Dependencies
- `freesasa>=2.1.0` (SASA calculation)
- HDOCK createpl binary dependency
- Docker integration for HDOCK tools

## Future Enhancements (Not Implemented)
- pLDDT score parsing from AlphaFold2 models
- pLDDT-based coloring in Molstar viewer
- PAE (Predicted Aligned Error) matrix integration
- AlphaFold-Multimer specific features

## Testing Recommendations

1. **Input Validation**
   - Test with directory containing only PDB files
   - Test with directory containing only mmCIF files
   - Test with mixed PDB/mmCIF directory
   - Test with invalid directory path
   - Test with empty directory

2. **Reference Model**
   - Test with valid reference PDB
   - Test with valid reference CIF
   - Test without reference model (should use first structure)
   - Test reference model superposition accuracy

3. **Clustering**
   - Test with small dataset (5-10 models)
   - Test with medium dataset (50-100 models)
   - Test with large dataset (200+ models)
   - Verify HDBSCAN parameters work correctly

4. **Visualization**
   - Verify all plots render without SASA data
   - Check interface quality cards show contact metrics
   - Ensure cluster summary table excludes SASA columns
   - Test Molstar viewer with AlphaFold2 models

## Migration Notes

- The app no longer requires HDOCK installation
- Processing is faster without SASA calculations
- Memory usage is reduced (no SASA precomputation)
- The focus is on interface contacts rather than buried surface area
- AlphaFold2's pLDDT scores provide uncertainty estimation (replaces SASA as quality metric)

## Breaking Changes

⚠️ **Users of the previous version should note:**
1. HDOCK .out files are no longer supported
2. SASA metrics are not available
3. pose extraction/createpl functionality removed
4. Different input format (directory instead of file uploads)
5. Default RMSD selection changed (protein CA instead of nucleic P)

## Files Modified Summary

```
core/
  io_utils.py          - Updated to support PDB/mmCIF search
app_main.py            - Major refactor (HDOCK→AF2, removed SASA)
requirements.txt       - Removed freesasa, updated description
README_APP.md          - Comprehensive rewrite for AF2
QUICKSTART.md          - New quick start guide
```

## Completion Status

✅ All planned changes completed
✅ HDOCK dependencies removed
✅ AlphaFold2 input system implemented
✅ SASA calculations and display removed
✅ Documentation updated
✅ Reference model support added
