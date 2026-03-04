"""
File I/O utilities for ClustalDM
"""

import os
import glob


def find_pdb_files(directories, filter_pattern=None):
    """
    Find all PDB files in given directories
    
    Parameters:
    -----------
    directories : list or str
        List of directory paths to search for PDB files, or single directory path
    filter_pattern : str, optional
        If provided, only include PDB files with this string in the filename
        (e.g., "model" to only get model_1.pdb, model_2.pdb, etc.)
    
    Returns:
    --------
    list : List of paths to PDB files
    """
    # Handle single directory string
    if isinstance(directories, str):
        directories = [directories]
    
    pdb_files = []
    
    for directory in directories:
        if os.path.isdir(directory):
            # Find all .pdb files in directory
            pattern = os.path.join(directory, "*.pdb")
            found_files = glob.glob(pattern)
            
            # Apply filter if specified
            if filter_pattern:
                found_files = [f for f in found_files if filter_pattern.lower() in os.path.basename(f).lower()]
                print(f"Found {len(found_files)} PDB files matching '{filter_pattern}' in {directory}")
            else:
                print(f"Found {len(found_files)} PDB files in {directory}")
            
            pdb_files.extend(found_files)
        else:
            print(f"Warning: {directory} is not a valid directory")
    
    return sorted(pdb_files)
