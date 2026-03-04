"""
HDOCK output parser and PDB generator

Parses HDOCK .out files and generates PDB structures for each docking pose
"""

import numpy as np
import os
import tempfile
from typing import Dict, List, Tuple
import MDAnalysis as mda
from scipy.spatial.transform import Rotation


class HDOCKParser:
    """Parse HDOCK output files and generate PDB structures for poses"""
    
    def __init__(self, hdock_out_file: str):
        """
        Initialize parser with HDOCK output file
        
        Parameters:
        -----------
        hdock_out_file : str
            Path to HDOCK .out file
        """
        self.out_file = hdock_out_file
        self.grid_spacing = None
        self.angle_step = None
        self.initial_rotation = None
        self.receptor_file = None
        self.receptor_center = None
        self.ligand_file = None
        self.ligand_center = None
        self.poses = []
        
        self._parse_header()
        self._parse_poses()
    
    def _parse_header(self):
        """Parse header lines (first 5 lines) of HDOCK output"""
        with open(self.out_file, 'r') as f:
            lines = f.readlines()
        
        # Line 1: Grid spacing
        parts = lines[0].strip().split()
        self.grid_spacing = float(parts[-1])
        
        # Line 2: Angle step
        parts = lines[1].strip().split()
        self.angle_step = float(parts[-1])
        
        # Line 3: Initial rotation (optional, skip for now)
        parts = lines[2].strip().split()
        self.initial_rotation = [float(x) for x in parts[-3:]]
        
        # Line 4: Receptor file and center
        parts = lines[3].strip().split()
        self.receptor_file = parts[0]
        self.receptor_center = np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
        
        # Line 5: Ligand file and center
        parts = lines[4].strip().split()
        self.ligand_file = parts[0]
        self.ligand_center = np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
    
    def _parse_poses(self):
        """Parse all pose lines (line 6 onwards)"""
        with open(self.out_file, 'r') as f:
            lines = f.readlines()[5:]  # Skip header
        
        for line_num, line in enumerate(lines, start=1):
            parts = line.strip().split()
            if len(parts) < 9:
                continue  # Skip invalid lines
            
            try:
                pose = {
                    'pose_id': line_num,
                    'translation': np.array([float(parts[0]), float(parts[1]), float(parts[2])]),
                    'rotation': np.array([float(parts[3]), float(parts[4]), float(parts[5])]),
                    'score': float(parts[6]),
                    'rmsd': float(parts[7]),
                    'trans_id': float(parts[8])
                }
                self.poses.append(pose)
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line {line_num + 5}: {e}")
                continue
    
    def apply_transformation(self, coords: np.ndarray, pose: dict) -> np.ndarray:
        """
        Apply HDOCK transformation to coordinates
        
        HDOCK transformation (reverse engineered from createpl output):
        - Translation values: x, y, z (ALREADY in Angstroms, NOT grid units!)
        - Rotation values: phi, theta, psi (Euler angles in radians, YZY sequence, NEGATED)
        
        Transformation steps:
        1. Center ligand at origin
        2. Apply rotation (YZY Euler angles, NEGATED: -phi, -theta, -psi)
        3. Apply direct translation (already in Angstroms)
        4. Add receptor center
        
        Parameters:
        -----------
        coords : np.ndarray
            Original ligand coordinates (N x 3)
        pose : dict
            Pose dictionary with 'translation' and 'rotation'
        
        Returns:
        --------
        np.ndarray : Transformed coordinates in receptor frame
        """
        # Step 1: Center ligand at its center of geometry
        centered = coords - self.ligand_center
        
        # Step 2: Apply rotation (YZY Euler angles, NEGATED!)
        rotation = Rotation.from_euler('YZY', -pose['rotation'], degrees=False)
        rotated = rotation.apply(centered)
        
        # Step 3: Apply translation + receptor center (no scaling!)
        translated = rotated + pose['translation'] + self.receptor_center
        
        return translated
    
    def generate_pose_pdb(self, receptor_pdb: str, ligand_pdb: str, pose_idx: int, 
                          output_dir: str = None) -> str:
        """
        Generate a PDB file for a specific pose
        
        Parameters:
        -----------
        receptor_pdb : str
            Path to receptor PDB file
        ligand_pdb : str
            Path to ligand PDB file
        pose_idx : int
            Index of pose to generate (0-based)
        output_dir : str, optional
            Output directory for PDB. If None, uses temp directory
        
        Returns:
        --------
        str : Path to generated PDB file
        """
        if pose_idx >= len(self.poses):
            raise ValueError(f"Pose index {pose_idx} out of range (max: {len(self.poses)-1})")
        
        pose = self.poses[pose_idx]
        
        # Load receptor and ligand
        receptor = mda.Universe(receptor_pdb)
        ligand = mda.Universe(ligand_pdb)
        
        # Transform ligand coordinates
        ligand_coords_transformed = self.apply_transformation(
            ligand.atoms.positions, pose
        )
        
        # Create output directory
        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix="hdock_poses_")
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate output filename
        output_file = os.path.join(output_dir, f"model_{pose_idx + 1}.pdb")
        
        # Standard chain renaming: Ligand chains → L, M, N, ...
        # This provides clear separation between receptor and ligand
        ligand_chains = sorted(set(ligand.atoms.chainIDs))
        ligand_chain_letters = [chr(ord('L') + i) for i in range(len(ligand_chains))]
        ligand_chain_mapping = dict(zip(ligand_chains, ligand_chain_letters))
        
        # Create merged structure with standard conventions
        with open(output_file, 'w') as out:
            # Write header with metadata
            out.write(f"REMARK   1 HDOCK POSE {pose_idx + 1}\n")
            out.write(f"REMARK   2 SCORE: {pose['score']:.2f}\n")
            out.write(f"REMARK   3 RECEPTOR CHAINS: {', '.join(sorted(set(receptor.atoms.chainIDs)))}\n")
            out.write(f"REMARK   4 LIGAND CHAINS: {', '.join(ligand_chain_letters)}\n")
            
            # Write receptor (keep original chains)
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
                receptor.atoms.write(tmp.name)
                with open(tmp.name, 'r') as rec_f:
                    for line in rec_f:
                        if line.startswith(('ATOM', 'HETATM')):
                            out.write(line)
                os.unlink(tmp.name)
            
            # Write TER to separate receptor from ligand
            out.write("TER\n")
            
            # Write transformed ligand (renumber chains to L, M, N...)
            ligand.atoms.positions = ligand_coords_transformed
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
                ligand.atoms.write(tmp.name)
                with open(tmp.name, 'r') as lig_f:
                    for line in lig_f:
                        if line.startswith(('ATOM', 'HETATM')):
                            # Update chain ID (PDB column 22, 0-indexed = 21)
                            old_chain = line[21] if len(line) > 21 else ' '
                            new_chain = ligand_chain_mapping.get(old_chain, 'L')
                            line = line[:21] + new_chain + line[22:]
                            out.write(line)
                os.unlink(tmp.name)
            
            out.write("END\n")
        
        return output_file
    
    def generate_all_poses(self, receptor_pdb: str, ligand_pdb: str, 
                          output_dir: str, max_poses: int = None) -> List[str]:
        """
        Generate PDB files for all (or top N) poses
        
        Parameters:
        -----------
        receptor_pdb : str
            Path to receptor PDB file
        ligand_pdb : str
            Path to ligand PDB file
        output_dir : str
            Output directory for PDB files
        max_poses : int, optional
            Maximum number of poses to generate (default: all)
        
        Returns:
        --------
        List[str] : List of generated PDB file paths
        """
        os.makedirs(output_dir, exist_ok=True)
        
        n_poses = len(self.poses) if max_poses is None else min(max_poses, len(self.poses))
        
        print(f"Generating {n_poses} pose PDB files...")
        pdb_files = []
        
        for i in range(n_poses):
            try:
                pdb_file = self.generate_pose_pdb(receptor_pdb, ligand_pdb, i, output_dir)
                pdb_files.append(pdb_file)
                
                if (i + 1) % 100 == 0:
                    print(f"  Generated {i + 1}/{n_poses} poses...")
            except Exception as e:
                print(f"  Warning: Failed to generate pose {i + 1}: {e}")
                continue
        
        print(f"✓ Generated {len(pdb_files)} PDB files in {output_dir}")
        return pdb_files
    
    def get_pose_info(self, pose_idx: int) -> Dict:
        """Get information about a specific pose"""
        if pose_idx >= len(self.poses):
            raise ValueError(f"Pose index {pose_idx} out of range")
        
        pose = self.poses[pose_idx]
        return {
            'pose_id': pose['pose_id'],
            'score': pose['score'],
            'rmsd': pose['rmsd'],
            'translation': pose['translation'].tolist(),
            'rotation': pose['rotation'].tolist()
        }
    
    def get_summary(self) -> Dict:
        """Get summary statistics of all poses"""
        if not self.poses:
            return {}
        
        scores = [p['score'] for p in self.poses]
        rmsds = [p['rmsd'] for p in self.poses]
        
        return {
            'total_poses': len(self.poses),
            'receptor': self.receptor_file,
            'ligand': self.ligand_file,
            'score_range': (min(scores), max(scores)),
            'score_mean': np.mean(scores),
            'rmsd_range': (min(rmsds), max(rmsds)),
            'rmsd_mean': np.mean(rmsds)
        }


# Convenience function for quick parsing
def parse_hdock_output(out_file: str, receptor_pdb: str, ligand_pdb: str, 
                       output_dir: str, max_poses: int = None) -> Tuple[HDOCKParser, List[str]]:
    """
    Parse HDOCK output and generate PDB files in one call
    
    Parameters:
    -----------
    out_file : str
        Path to HDOCK .out file
    receptor_pdb : str
        Path to receptor PDB file
    ligand_pdb : str
        Path to ligand PDB file
    output_dir : str
        Output directory for generated PDBs
    max_poses : int, optional
        Maximum number of poses to generate
    
    Returns:
    --------
    Tuple[HDOCKParser, List[str]] : Parser object and list of generated PDB files
    """
    parser = HDOCKParser(out_file)
    summary = parser.get_summary()
    
    print(f"\n{'='*60}")
    print(f"HDOCK Output Summary")
    print(f"{'='*60}")
    print(f"Receptor: {summary['receptor']}")
    print(f"Ligand: {summary['ligand']}")
    print(f"Total poses: {summary['total_poses']}")
    print(f"Score range: {summary['score_range'][0]:.2f} to {summary['score_range'][1]:.2f}")
    print(f"RMSD range: {summary['rmsd_range'][0]:.2f} to {summary['rmsd_range'][1]:.2f}")
    print(f"{'='*60}\n")
    
    pdb_files = parser.generate_all_poses(receptor_pdb, ligand_pdb, output_dir, max_poses)
    
    return parser, pdb_files
