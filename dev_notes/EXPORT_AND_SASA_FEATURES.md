# Export and ΔSASA Features

## Export Functionality

### UI Components (app_main.py)
- **Output Directory**: Text input with browse button (uses tkinter file dialog)
- **Data Export Checkboxes**:
  - Cluster Info (cluster_info.txt)
  - Cluster Statistics (cluster_statistics.csv)
  - Cluster Summary (cluster_summary.csv)
  - RMSD Summary (cluster_rmsd_summary.csv)
- **Plot Export Checkboxes**:
  - Cluster Visualization
  - Size Distribution
  - k-NN Distance Plot
  - Heatmaps (Jaccard + RMSD)

### Export Handler
- Creates output directory if it doesn't exist
- Copies selected CSV files to export directory
- Exports selected plots as PNG (1000-1200px width)
- Shows detailed status with timestamp and file list
- Error handling with user-friendly messages

### Requirements
- Added `kaleido` to requirements.txt for plotly static image export

## ΔSASA Calculation

### Purpose
Calculate buried surface area upon protein-DNA binding as additional QC metric beyond RMSD.

### Formula
```
ΔSASA = SASA(complex) - [SASA(protein) + SASA(DNA)]
```

### Two Modes

#### 1. Whole Protein ΔSASA
- Calculates total buried surface area for entire protein-DNA interface
- Returns absolute value (area is negative, we want magnitude)
- Components: sasa_complex, sasa_protein, sasa_nucleic

#### 2. Motif-Specific ΔSASA
- Calculates buried surface area for specific residue ranges
- Input format: `"Chain:Start-End, Chain:Start-End"` 
- Example: `"A:10-20, B:30-35"`
- Useful for focusing on known binding motifs or DNA-binding domains

### Implementation (core/analysis.py)

#### Functions Added:
1. **`calculate_delta_sasa(pdb_file, protein_chains, nucleic_chains, motif_residues=None)`**
   - Uses FreeSASA library
   - Returns dict with total_delta_sasa, motif_delta_sasa, and component SASAs
   - Handles chain selections for protein and nucleic acid
   
2. **`parse_motif_residues(motif_string)`**
   - Parses user input: "A:10-20, B:30-35"
   - Returns dict: {'A': [10,11,...,20], 'B': [30,31,...,35]}
   - Handles ranges and individual residues

### Usage in App

```python
from core.analysis import InterfaceAnalyzer

# Parse motif input
motif_dict = InterfaceAnalyzer.parse_motif_residues("A:10-20")

# Calculate ΔSASA
sasa_results = InterfaceAnalyzer.calculate_delta_sasa(
    pdb_file="model.pdb",
    protein_chains=['A', 'B'],
    nucleic_chains=['C', 'D'],
    motif_residues=motif_dict  # Optional
)

# Access results
total_delta = sasa_results['total_delta_sasa']  # Full interface
motif_delta = sasa_results['motif_delta_sasa']  # Motif only (or None)
```

### Next Steps

To fully integrate ΔSASA into the app:

1. **Add to analysis pipeline** - Calculate ΔSASA during initial PDB processing
2. **Store results** - Add to quality metrics dataframe
3. **Display in UI** - Show in Quality Metrics cards and cluster summary table
4. **Add to exports** - Include ΔSASA in cluster statistics CSV
5. **Add tooltips** - Explain what ΔSASA means and how to interpret values

### Interpretation

- **Higher ΔSASA** = More surface area buried = Stronger binding interface
- **Typical values**: 500-2000 Ų for protein-DNA complexes
- **Motif ΔSASA**: Compare to total to see motif's contribution
- **Quality check**: ΔSASA should be consistent within a cluster if structures are similar

### Dependencies
- FreeSASA 2.1.0 (already in requirements.txt)
- Works with standard PDB format
- Chain IDs must match PDB file

## Testing

### Export Feature
1. Run analysis to generate results
2. Go to Export tab
3. Select output directory (or browse)
4. Check desired files
5. Click Export
6. Verify files in output directory

### ΔSASA Calculation
```python
# Test with a sample PDB
from core.analysis import InterfaceAnalyzer

# Whole protein
result = InterfaceAnalyzer.calculate_delta_sasa(
    "test.pdb",
    protein_chains=['A'],
    nucleic_chains=['B']
)
print(f"Total ΔSASA: {result['total_delta_sasa']:.2f} Ų")

# With motif
motif = InterfaceAnalyzer.parse_motif_residues("A:50-70")
result = InterfaceAnalyzer.calculate_delta_sasa(
    "test.pdb",
    protein_chains=['A'],
    nucleic_chains=['B'],
    motif_residues=motif
)
print(f"Motif ΔSASA: {result['motif_delta_sasa']:.2f} Ų")
```

## Files Modified

1. **app_main.py**
   - Added Export tab UI (checkboxes, directory input, browse button)
   - Added motif_residues input field
   - Added export_handler() with file copying and plot export
   - Added browse_export_dir() handler
   - Added export_status output

2. **core/analysis.py**
   - Added calculate_delta_sasa() method
   - Added parse_motif_residues() helper

3. **requirements.txt**
   - Added kaleido for plotly image export
