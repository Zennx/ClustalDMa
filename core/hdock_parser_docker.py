"""
HDOCK Parser using official createpl binary via Docker

Since HDOCK's transformation format is proprietary and cannot be reliably 
reverse-engineered, this parser uses the official createpl binary via Docker.
"""

import os
import subprocess
import shutil
from typing import List
import MDAnalysis as mda


class HDOCKParser:
    """Parse HDOCK output using the official createpl tool via Docker"""
    
    def __init__(self, hdock_out_file: str, createpl_path: str = None):
        """
        Initialize parser
        
        Parameters:
        -----------
        hdock_out_file : str
            Path to HDOCK .out file
        createpl_path : str, optional
            Path to createpl binary (default: auto-detect)
        """
        self.out_file = os.path.abspath(hdock_out_file)
        
        # Find createpl binary
        if createpl_path:
            self.createpl = os.path.abspath(createpl_path)
        else:
            # Try multiple locations
            search_paths = [
                '/app/hdock_tools/createpl_linux',  # Docker container
                os.path.join(os.path.dirname(self.out_file), 'createpl_linux'),
                os.path.join(os.path.dirname(os.path.dirname(self.out_file)), 'createpl_linux'),
                './createpl_linux',
            ]
            
            self.createpl = None
            for path in search_paths:
                if os.path.exists(path):
                    self.createpl = os.path.abspath(path)
                    break
        
        if not self.createpl or not os.path.exists(self.createpl):
            raise FileNotFoundError(
                f"createpl_linux binary not found. Searched: {search_paths}\n"
                f"Please provide createpl_path or place binary in hdock_tools/"
            )
        
        # Check if we're inside Docker or need Docker
        self.use_docker = not os.path.exists('/app/hdock_tools')  # Outside container
        
        # Parse header to get file info
        self._parse_header()
    
    def _parse_header(self):
        """Parse HDOCK output header"""
        with open(self.out_file, 'r') as f:
            lines = f.readlines()
        
        # Get receptor and ligand filenames from header
        self.receptor_file = lines[3].split()[0]
        self.ligand_file = lines[4].split()[0]
        
        # Count poses
        self.num_poses = len(lines) - 5  # Header is 5 lines
        
        # Parse first pose to get score
        if self.num_poses > 0:
            first_pose = lines[5].split()
            self.best_score = float(first_pose[6])
    
    def generate_all_poses(self, output_dir: str, max_poses: int = None,
                          receptor_pdb: str = None, ligand_pdb: str = None) -> List[str]:
        """
        Generate PDB files for all poses using createpl via Docker
        
        Parameters:
        -----------
        output_dir : str
            Output directory for PDB files
        max_poses : int, optional
            Maximum number of poses to generate (default: all)
        receptor_pdb : str, optional
            Path to receptor PDB (if different from .out header)
        ligand_pdb : str, optional
            Path to ligand PDB (if different from .out header)
        
        Returns:
        --------
        List[str] : List of generated PDB file paths
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Determine full paths to receptor and ligand
        out_dir = os.path.dirname(self.out_file)
        
        if receptor_pdb is None:
            receptor_pdb = os.path.join(out_dir, self.receptor_file)
        if ligand_pdb is None:
            ligand_pdb = os.path.join(out_dir, self.ligand_file)
        
        # Run createpl (directly if in container, via Docker if outside)
        nmax = max_poses if max_poses else self.num_poses
        
        # Prepare working directory (use temp dir for createpl execution)
        import tempfile
        with tempfile.TemporaryDirectory() as work_dir:
            out_basename = os.path.basename(self.out_file)
            
            # Copy input files to work dir with correct names
            # createpl expects filenames matching the .out header
            shutil.copy(self.out_file, work_dir)
            shutil.copy(receptor_pdb, os.path.join(work_dir, self.receptor_file))
            shutil.copy(ligand_pdb, os.path.join(work_dir, self.ligand_file))
            
            if self.use_docker:
                # Outside container - use Docker
                self._run_via_docker(work_dir, out_basename, nmax)
            else:
                # Inside container - run directly
                self._run_native(work_dir, out_basename, nmax)
            
            # Move generated models to output directory
            generated_files = []
            for i in range(1, nmax + 1):
                src = os.path.join(work_dir, f'model_{i}.pdb')
                if os.path.exists(src):
                    # Rename chains to standard: receptor stays, ligand → L
                    dst = os.path.join(output_dir, f'model_{i}.pdb')
                    self._standardize_chains(src, dst)
                    generated_files.append(dst)
            
            return generated_files
    
    def _run_native(self, work_dir: str, out_basename: str, nmax: int):
        """Run createpl directly (inside container)"""
        cmd = [
            self.createpl,
            out_basename,
            'model.pdb',
            '-complex',
            '-models',
            '-nmax', str(nmax)
        ]
        
        print(f"Generating {nmax} poses using createpl...")
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"createpl failed: {result.stderr}")
    
    def _run_via_docker(self, work_dir: str, out_basename: str, nmax: int):
        """Run createpl via Docker (outside container)"""
        # Check if Docker is available
        try:
            subprocess.run(['docker', 'version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(
                "Docker is required. Please install Docker or run inside container."
            )
        
        cmd = [
            'docker', 'run', '--rm',
            '-v', f'{work_dir}:/data',
            '-w', '/data',
            '--platform', 'linux/amd64',
            'centos:7',
            '/data/createpl_linux',
            out_basename,
            'model.pdb',
            '-complex',
            '-models',
            '-nmax', str(nmax)
        ]
        
        print(f"Generating {nmax} poses using createpl via Docker...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"createpl failed: {result.stderr}")
    
    def _standardize_chains(self, src_pdb: str, dst_pdb: str):
        """
        Standardize chain IDs: keep receptor chains, rename ligand to L
        """
        # Load structure
        u = mda.Universe(src_pdb)
        
        # Separate protein/DNA (receptor) from ligand by residue names
        # HDOCK output has all chains as 'A', need to distinguish by position
        # Ligand is typically at the end
        
        # For now, just copy the file and add chain L designation in remarks
        with open(src_pdb, 'r') as f:
            lines = f.readlines()
        
        with open(dst_pdb, 'w') as out:
            # Get receptor and ligand chain info from REMARK
            for line in lines:
                if line.startswith('REMARK'):
                    out.write(line)
            
            # Write atoms - ligand atoms come after TER usually
            in_ligand = False
            for line in lines:
                if line.startswith('MODEL'):
                    continue  # Skip MODEL lines
                elif line.startswith('ENDMDL'):
                    continue  # Skip ENDMDL
                elif line.startswith('ATOM') or line.startswith('HETATM'):
                    # Change chain to L for nucleic acids (ligand)
                    resname = line[17:20].strip()
                    if resname.startswith('D'):  # DNA residues
                        line = line[:21] + 'L' + line[22:]
                    out.write(line)
                elif line.startswith('TER'):
                    out.write(line)
                elif line.startswith('END'):
                    out.write(line)
        
        print(f"  ✓ {os.path.basename(dst_pdb)}")


def parse_hdock_output(hdock_file: str) -> HDOCKParser:
    """
    Convenience function to create parser
    
    Parameters:
    -----------
    hdock_file : str
        Path to HDOCK .out file
    
    Returns:
    --------
    HDOCKParser instance
    """
    return HDOCKParser(hdock_file)
