# HDOCK Output Parser

## Overview

This module parses HDOCK docking output files (`.out` format) and generates PDB structures for each docking pose. It enables direct analysis of HDOCK results without manual conversion.

## HDOCK Output Format

HDOCK `.out` files contain:
- **Header (5 lines)**:
  1. Grid spacing for translational search
  2. Angle step for rotational search  
  3. Initial rotation (optional)
  4. Receptor file and center of geometry (x, y, z)
  5. Ligand file and center of geometry (x, y, z)

- **Poses (line 6 onwards)**: Each line represents one docking solution
  - 3 translations (x, y, z)
  - 3 rotations (Euler angles in radians, ZYZ convention)
  - Binding score
  - RMSD from initial ligand orientation
  - Translational ID

## Usage

### Quick Start

```python
from core.hdock_parser import parse_hdock_output

# Parse HDOCK output and generate PDB files
parser, pdb_files = parse_hdock_output(
    out_file='hdock.out',
    receptor_pdb='protein.pdb',
    ligand_pdb='dna.pdb',
    output_dir='./poses/',
    max_poses=100  # Optional: limit to top 100 poses
)

print(f"Generated {len(pdb_files)} PDB files")
```

### Command Line

```bash
python demo_hdock_parser.py hdock.out receptor.pdb ligand.pdb output_dir/ [max_poses]
```

Example:
```bash
python demo_hdock_parser.py \\
    hdock_results.out \\
    protein.pdb \\
    dsDNA.pdb \\
    ./docking_poses/ \\
    500
```

### Advanced Usage

```python
from core.hdock_parser import HDOCKParser

# Parse HDOCK output
parser = HDOCKParser('hdock.out')

# Get summary statistics
summary = parser.get_summary()
print(f"Total poses: {summary['total_poses']}")
print(f"Score range: {summary['score_range']}")

# Generate specific pose
pdb_file = parser.generate_pose_pdb(
    receptor_pdb='protein.pdb',
    ligand_pdb='dna.pdb',
    pose_idx=0,  # Top-scoring pose
    output_dir='./poses/'
)

# Get pose information
info = parser.get_pose_info(0)
print(f"Pose 1 score: {info['score']}")
print(f"Pose 1 RMSD: {info['rmsd']}")

# Generate all poses
pdb_files = parser.generate_all_poses(
    receptor_pdb='protein.pdb',
    ligand_pdb='dna.pdb',
    output_dir='./all_poses/',
    max_poses=1000
)
```

## Integration with ClustalDM

Once PDB files are generated, analyze them with ClustalDM:

```bash
# 1. Generate poses from HDOCK output
python demo_hdock_parser.py hdock.out protein.pdb dna.pdb ./poses/ 500

# 2. Launch ClustalDM app
./launch_app.sh

# 3. In the app:
#    - Set PDB Directory: ./poses/
#    - Optional: Define motif residues for screening
#    - Click "Run Analysis"
```

## Workflow Example

Complete workflow for HDOCK → ClustalDM analysis:

```python
from core.hdock_parser import parse_hdock_output
from core.clusterer import PDBClusterer

# Step 1: Parse HDOCK output and generate PDBs
parser, pdb_files = parse_hdock_output(
    out_file='hdock.out',
    receptor_pdb='7vxr_protein.pdb',
    ligand_pdb='dsDNA.pdb',
    output_dir='./hdock_poses/',
    max_poses=500  # Top 500 poses
)

# Step 2: Cluster poses using ClustalDM
clust = PDBClusterer(pdb_files)

# Optional: Filter by motif before clustering
motif = {'A': [10, 11, 12], 'B': [30, 31, 32]}
clust.compute_jaccard_contact_matrix(n_jobs=-1, motif_residues=motif)

# Cluster with DBSCAN
clust.cluster_dbscan(eps=0.3, min_samples=3)

# Export results
clust.export_clusters('./clustered_poses/')
```

## Performance Tips

### For Large HDOCK Outputs (>1000 poses)

1. **Limit initial generation**: Start with top 500-1000 poses
   ```python
   parser, pdb_files = parse_hdock_output(..., max_poses=500)
   ```

2. **Use motif screening**: Filter off-target poses early
   ```python
   motif = {'A': [10, 11, 12]}  # Key binding residues
   clust.compute_jaccard_contact_matrix(motif_residues=motif)
   ```

3. **Parallel processing**: Use all CPU cores
   ```python
   clust.compute_jaccard_contact_matrix(n_jobs=-1)
   ```

### Memory Optimization

For 4000+ poses, the parser generates PDBs incrementally to avoid loading all structures in memory at once.

## Troubleshooting

### Common Issues

**Issue**: "Rotation produces incorrect structures"
- **Solution**: Check HDOCK Euler angle convention. Current implementation uses ZYZ convention. Adjust in `apply_transformation()` if needed.

**Issue**: "Generated PDBs are misaligned"
- **Solution**: Verify receptor and ligand PDB files match those used in HDOCK docking. Centers of geometry must be consistent.

**Issue**: "Out of memory with 4000 poses"
- **Solution**: Generate poses in batches:
  ```python
  for i in range(0, 4000, 500):
      parser.generate_all_poses(..., max_poses=i+500)
  ```

## Technical Details

### Coordinate Transformation

HDOCK transformation applied to ligand:

1. **Center at origin**: `coords - ligand_center`
2. **Apply rotation**: Euler angles (ZYZ convention, radians)
3. **Translate**: `rotated + translation`
4. **Move to receptor**: `translated + receptor_center`

### File Structure

Generated PDB files:
- `model_1.pdb`: Top-scoring pose
- `model_2.pdb`: Second-best pose
- ...
- `model_N.pdb`: Nth pose

Each PDB contains:
- Receptor atoms (original coordinates)
- Ligand atoms (transformed coordinates)
- END record

## Requirements

- MDAnalysis
- NumPy
- SciPy (for Rotation)

## See Also

- HDOCK server: http://hdock.phys.hust.edu.cn/
- HDOCK paper: Yan et al., Bioinformatics, 2020
- ClustalDM documentation: See README.md
