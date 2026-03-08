"""
File I/O utilities for ClustalDM
"""

import os
import glob
import re
from pathlib import Path


def parse_alphafold_filename(filename):
    """
    Extract receptor start/end residue numbers from AlphaFold2 multimer output filename.
    
    Filename pattern:
    RECEPTOR_NAME_START_ENDLIGAND_NAME_START_END_a878f_unrelaxed_rank_XXX_alphafold2_multimer_v3_model_X_seed_XXX.pdb
    
    Example:
    ELOA1_HUMAN_1_437SNAI1_HUMAN_1_265_a878f_unrelaxed_rank_003_alphafold2_multimer_v3_model_3_seed_000.pdb
    
    Parameters:
    -----------
    filename : str
        AlphaFold2 output filename (basename only, not full path)
    
    Returns:
    --------
    dict or None : Dictionary with 'receptor_start', 'receptor_end', 'receptor_name' if parsing succeeds, else None
    """
    # Pattern: RECEPTOR_NAME_START_ENDLIGAND_NAME...
    # The END number is merged with LIGAND_NAME (no underscore between them)
    # We match until we hit a capital letter indicating the start of LIGAND_NAME
    pattern = r'^([A-Z0-9_]+)_(\d+)_(\d+)([A-Z][A-Z0-9_]+)_'
    
    match = re.match(pattern, filename)
    if match:
        receptor_name = match.group(1)
        receptor_start = int(match.group(2))
        receptor_end = int(match.group(3))
        
        return {
            'receptor_name': receptor_name,
            'receptor_start': receptor_start,
            'receptor_end': receptor_end,
            'offset': receptor_start - 1  # PDB residues start at 1, so offset = start - 1
        }
    
    # Fallback: try simpler pattern for non-standard names
    # Just look for two numbers separated by underscores before "alphafold" keyword
    simple_pattern = r'_(\d+)_(\d+)[A-Z].*alphafold'
    match = re.search(simple_pattern, filename)
    if match:
        receptor_start = int(match.group(1))
        receptor_end = int(match.group(2))
        return {
            'receptor_name': 'UNKNOWN',
            'receptor_start': receptor_start,
            'receptor_end': receptor_end,
            'offset': receptor_start - 1
        }
    
    # No match - return None (will skip offset correction)
    return None


def read_fasta_sequence(fasta_path):
    """
    Read sequence from FASTA file.
    
    Parameters:
    -----------
    fasta_path : str
        Path to FASTA file
    
    Returns:
    --------
    str : Amino acid sequence (single-letter codes, uppercase)
    """
    sequence = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                # Skip header lines
                continue
            sequence.append(line.upper())
    
    return ''.join(sequence)


def extract_sequence_from_pdb(pdb_path, chain_id='A'):
    """
    Extract amino acid sequence from PDB file for a specific chain.
    
    Parameters:
    -----------
    pdb_path : str
        Path to PDB file
    chain_id : str
        Chain ID to extract (default: 'A' for receptor)
    
    Returns:
    --------
    str : Amino acid sequence (single-letter codes)
    """
    import MDAnalysis as mda
    
    # Three-letter to one-letter code mapping
    aa_codes = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'MSE': 'M'  # Selenomethionine
    }
    
    try:
        u = mda.Universe(pdb_path)
        # Select protein residues from specified chain
        selection = u.select_atoms(f"protein and chainID {chain_id}")
        
        # Get unique residues in order
        residues = selection.residues
        sequence = []
        
        for res in residues:
            resname = res.resname.strip()
            if resname in aa_codes:
                sequence.append(aa_codes[resname])
            else:
                sequence.append('X')  # Unknown residue
        
        return ''.join(sequence)
    except Exception as e:
        print(f"Warning: Could not extract sequence from {pdb_path}: {e}")
        return ""


def align_sequences(seq1, seq2, mode='local'):
    """
    Perform pairwise sequence alignment.
    Uses Smith-Waterman (local) for chopped sequences or Needleman-Wunsch (global).
    
    Parameters:
    -----------
    seq1 : str
        First sequence (e.g., reference full-length)
    seq2 : str
        Second sequence (e.g., chopped AlphaFold fragment)
    mode : str
        'local' for Smith-Waterman (better for fragments), 'global' for Needleman-Wunsch
    
    Returns:
    --------
    dict : {'score': float, 'identity': float, 'start_pos': int, 'alignment': str}
           start_pos is the 1-based position in seq1 where seq2 aligns best
    """
    try:
        from Bio import Align
        
        # Use Biopython's aligner
        aligner = Align.PairwiseAligner()
        aligner.mode = mode  # 'local' or 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        
        alignments = list(aligner.align(seq1, seq2))
        if len(alignments) == 0:
            raise ValueError("No alignment found")
        
        best_alignment = alignments[0]
        
        # Extract aligned sequences
        alignment_str = str(best_alignment)
        lines = alignment_str.split('\n')
        
        # Calculate identity (matches / alignment length)
        if len(lines) >= 3:
            aligned_seq1 = lines[0]
            aligned_seq2 = lines[2]
            matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
            alignment_length = len([c for c in aligned_seq2 if c != '-'])
            identity = (matches / alignment_length * 100) if alignment_length > 0 else 0.0
        else:
            identity = 0.0
        
        # Find start position in seq1 (where alignment begins, 1-based)
        # For local alignment, this is where the local match starts
        aligned_seq1_full = str(best_alignment).split('\n')[0]
        start_pos = 1
        for i, char in enumerate(aligned_seq1_full):
            if char != '-':
                start_pos = i + 1
                break
        
        return {
            'score': best_alignment.score,
            'identity': identity,
            'start_pos': start_pos,
            'alignment': alignment_str,
            'mode': mode
        }
    
    except ImportError:
        # Fallback: simple exact substring search
        print("Warning: Biopython not available, using exact substring matching")
        pos = seq1.find(seq2)
        if pos >= 0:
            return {
                'score': len(seq2),
                'identity': 100.0,
                'start_pos': pos + 1,  # 1-based
                'alignment': 'Exact match',
                'mode': 'exact'
            }
        else:
            # No exact match - try to find best approximate position
            max_matches = 0
            best_pos = 0
            for i in range(len(seq1) - len(seq2) + 1):
                matches = sum(1 for a, b in zip(seq1[i:i+len(seq2)], seq2) if a == b)
                if matches > max_matches:
                    max_matches = matches
                    best_pos = i
            
            identity = max_matches / len(seq2) * 100
            return {
                'score': max_matches,
                'identity': identity,
                'start_pos': best_pos + 1,
                'alignment': f'Approximate match at position {best_pos + 1}',
                'mode': 'approximate'
            }


def validate_and_correct_residue_offset(pdb_path, reference_sequence=None, cache=None):
    """
    Validate residue numbering using filename parsing and optional sequence alignment.
    
    Strategy:
    1. Parse receptor start/end from AlphaFold filename (fast)
    2. If reference sequence provided, extract PDB sequence and align to validate
    3. Cache results to avoid re-computing for duplicate structures
    
    Parameters:
    -----------
    pdb_path : str
        Path to PDB file
    reference_sequence : str, optional
        Reference full-length amino acid sequence for validation
    cache : dict, optional
        Dictionary to cache results (key: pdb_basename, value: offset_info)
    
    Returns:
    --------
    dict : {
        'offset': int,           # Residue number offset (0 if no correction needed)
        'method': str,          # 'filename' or 'alignment'
        'validated': bool,      # True if alignment confirmed filename parsing
        'identity': float,      # Sequence identity % (if alignment performed)
        'receptor_start': int,  # Expected start position
        'receptor_end': int     # Expected end position
    }
    """
    filename = os.path.basename(pdb_path)
    
    # Check cache first
    if cache is not None and filename in cache:
        return cache[filename]
    
    # Parse filename
    parsed = parse_alphafold_filename(filename)
    
    if parsed is None:
        # Could not parse filename - no offset correction
        result = {
            'offset': 0,
            'method': 'none',
            'validated': False,
            'identity': 0.0,
            'receptor_start': 1,
            'receptor_end': 0
        }
        if cache is not None:
            cache[filename] = result
        return result
    
    # If no reference sequence, just use filename offset
    if reference_sequence is None:
        result = {
            'offset': parsed['offset'],
            'method': 'filename',
            'validated': False,
            'identity': 0.0,
            'receptor_start': parsed['receptor_start'],
            'receptor_end': parsed['receptor_end']
        }
        if cache is not None:
            cache[filename] = result
        return result
    
    # Extract sequence from PDB and align to reference
    pdb_sequence = extract_sequence_from_pdb(pdb_path, chain_id='A')
    
    if not pdb_sequence:
        # Could not extract sequence - fall back to filename
        result = {
            'offset': parsed['offset'],
            'method': 'filename',
            'validated': False,
            'identity': 0.0,
            'receptor_start': parsed['receptor_start'],
            'receptor_end': parsed['receptor_end']
        }
        if cache is not None:
            cache[filename] = result
        return result
    
    # Align to reference (use local alignment for chopped fragments)
    alignment_result = align_sequences(reference_sequence, pdb_sequence, mode='local')
    
    # Check if alignment confirms filename parsing
    expected_start = parsed['receptor_start']
    aligned_start = alignment_result['start_pos']
    
    # Allow small discrepancy (±2 residues) due to alignment uncertainty
    validated = abs(expected_start - aligned_start) <= 2
    
    # Use alignment-based offset if identity is high, otherwise use filename
    if alignment_result['identity'] >= 90.0:
        final_offset = aligned_start - 1
        method = 'alignment'
    else:
        final_offset = parsed['offset']
        method = 'filename'
        print(f"Warning: Low sequence identity ({alignment_result['identity']:.1f}%) for {filename}, using filename offset")
    
    result = {
        'offset': final_offset,
        'method': method,
        'validated': validated,
        'identity': alignment_result['identity'],
        'receptor_start': parsed['receptor_start'],
        'receptor_end': parsed['receptor_end'],
        'aligned_start': aligned_start
    }
    
    if cache is not None:
        cache[filename] = result
    
    return result


def find_pdb_files(directories, filter_pattern=None):
    """
    Find all PDB and mmCIF files in given directories (AlphaFold2 models)
    Recursively searches subdirectories.
    
    Parameters:
    -----------
    directories : list or str
        List of directory paths to search for structure files, or single directory path
    filter_pattern : str, optional
        If provided, only include files with this string in the filename
        (e.g., "model" to only get model_1.pdb, model_2.pdb, etc.)
    
    Returns:
    --------
    list : List of paths to structure files (PDB and mmCIF)
    """
    # Handle single directory string
    if isinstance(directories, str):
        directories = [directories]
    
    structure_files = []
    
    for directory in directories:
        if os.path.isdir(directory):
            # Recursively find all .pdb and .cif files in directory and subdirectories
            pdb_pattern = os.path.join(directory, "**", "*.pdb")
            cif_pattern = os.path.join(directory, "**", "*.cif")
            found_files = glob.glob(pdb_pattern, recursive=True) + glob.glob(cif_pattern, recursive=True)
            
            # Apply filter if specified
            if filter_pattern:
                found_files = [f for f in found_files if filter_pattern.lower() in os.path.basename(f).lower()]
                print(f"Found {len(found_files)} structure files matching '{filter_pattern}' in {directory} (including subdirectories)")
            else:
                print(f"Found {len(found_files)} structure files (PDB/mmCIF) in {directory} (including subdirectories)")
            
            structure_files.extend(found_files)
        else:
            print(f"Warning: {directory} is not a valid directory")
    
    return sorted(structure_files)
