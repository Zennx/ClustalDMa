"""
Sequence alignment visualization with secondary structure and pLDDT confidence
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
from pathlib import Path


def is_alphafold_file(structure_path):
    """
    Heuristic check whether a structure file is from AlphaFold.
    Used to decide whether B-factor column should be treated/displayed as pLDDT.
    """
    if not structure_path:
        return False

    filename = os.path.basename(str(structure_path)).lower()
    if (
        'alphafold' in filename
        or 'af2' in filename
        or 'af3' in filename
        or filename.startswith('af-')
        or 'ranked_' in filename
        or 'unrelaxed' in filename
        or 'relaxed' in filename
    ):
        return True

    # Fallback to project parser if available
    try:
        from core.io_utils import parse_alphafold_filename
        return parse_alphafold_filename(filename) is not None
    except Exception:
        return False


def get_protein_chains(pdb_path, ligand_chain=None):
    """
    Get list of protein chain IDs from a structure file.
    
    Parameters:
    -----------
    pdb_path : str
        Path to PDB/CIF file
    ligand_chain : str, optional
        Chain to exclude (e.g., ligand chain)
    
    Returns:
    --------
    list : Chain IDs of protein chains
    """
    from Bio.PDB import PDBParser, MMCIFParser
    
    try:
        # Use appropriate parser
        if str(pdb_path).lower().endswith('.cif'):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        
        structure = parser.get_structure('protein', pdb_path)
        model = structure[0]
        
        # Standard amino acids
        aa_codes = {
            'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',
            'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER',
            'THR', 'VAL', 'TRP', 'TYR', 'MSE'
        }
        
        protein_chains = []
        for chain in model.get_chains():
            if ligand_chain and chain.id == ligand_chain:
                continue
            # Check if chain has protein residues
            has_protein = any(res.resname.strip() in aa_codes and res.id[0] == ' ' 
                            for res in chain)
            if has_protein:
                protein_chains.append(chain.id)
        
        return protein_chains
    except Exception as e:
        print(f"Warning: Could not extract chain info from {pdb_path}: {e}")
        return []


def extract_secondary_structure_dssp(pdb_path, chain_id='A', ligand_chain=None):
    """
    Extract secondary structure using MDAnalysis DSSP.
    
    Parameters:
    -----------
    pdb_path : str
        Path to PDB file
    chain_id : str or list
        Chain ID(s) to analyze (default: 'A')
        Can be single chain 'A' or list ['A', 'B', 'C'] for multi-chain receptors
    ligand_chain : str, optional
        If provided, analyze all protein chains EXCEPT this one (for multi-chain receptors)
    
    Returns:
    --------
    dict : {
        'sequence': str,
        'ss': list of str (secondary structure codes),
        'residue_numbers': list of int,
        'plddt': list of float (from B-factor)
    }
    """
    try:
        # Check if DSSP executable is available
        import shutil
        dssp_exe = shutil.which('mkdssp') or shutil.which('dssp')
        
        if dssp_exe:
            # Try BioPython DSSP first (more accurate - preserves H/G/I distinction)
            print(f"  Using BioPython DSSP (found: {dssp_exe})")
            from Bio.PDB import PDBParser, MMCIFParser, DSSP as BioDSSP
            
            # Use appropriate parser
            if str(pdb_path).lower().endswith('.cif'):
                parser = MMCIFParser(QUIET=True)
            else:
                parser = PDBParser(QUIET=True)
            
            structure = parser.get_structure('protein', pdb_path)
            model = structure[0]
            
            # Run BioPython DSSP
            dssp = BioDSSP(model, pdb_path)
            
            # Extract data
            sequence = []
            ss_codes = []
            residue_numbers = []
            plddt_values = []
            
            aa_codes = {
                'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
                'MSE': 'M'
            }
            
            # Determine which chains to process
            if ligand_chain:
                chains_to_process = [ch for ch in model.get_chains() if ch.id != ligand_chain]
            elif isinstance(chain_id, list):
                chains_to_process = [model[ch] for ch in chain_id if ch in model]
            else:
                # Single chain - try segid first, then chainID
                if chain_id in model:
                    chains_to_process = [model[chain_id]]
                else:
                    chains_to_process = list(model.get_chains())[:1]
            
            chain_ids = []  # Track which chain each residue belongs to
            
            for chain in chains_to_process:
                for residue in chain:
                    if residue.id[0] == ' ':  # Standard residue
                        dssp_key = (chain.id, residue.id)
                        if dssp_key in dssp:
                            ss_code = dssp[dssp_key][2]  # Secondary structure code (H/G/I/E/B/T/S/C/-)
                            resname = residue.resname.strip()
                            if resname in aa_codes:
                                sequence.append(aa_codes[resname])
                                ss_codes.append(ss_code)
                                residue_numbers.append(residue.id[1])
                                chain_ids.append(chain.id)  # Track chain ID
                                # Get B-factor from DSSP or atoms
                                bfactors = [atom.get_bfactor() for atom in residue.get_atoms()]
                                plddt_values.append(np.mean(bfactors) if bfactors else 0.0)
            
            print(f"  ✓ BioPython DSSP: {len(sequence)} residues, first 30 SS: {''.join(ss_codes[:30])}...")
            return {
                'sequence': ''.join(sequence),
                'ss': ss_codes,
                'residue_numbers': residue_numbers,
                'plddt': plddt_values,
                'chain_ids': chain_ids
            }
        else:
            # DSSP executable not found, skip to MDAnalysis
            print(f"  DSSP executable not found, using MDAnalysis (note: lumps H/G/I together)")
            raise FileNotFoundError("mkdssp not available")
    
    except Exception as bio_e:
        # BioPython DSSP failed or not available, try MDAnalysis DSSP
        if not isinstance(bio_e, FileNotFoundError):
            print(f"  BioPython DSSP failed: {bio_e}")
        print(f"  Trying MDAnalysis DSSP...")
        try:
            from MDAnalysis import Universe
            from MDAnalysis.analysis.dssp import DSSP as MDA_DSSP
            u = Universe(str(pdb_path))
            
            # Select protein atoms, excluding ligand chain if specified
            if ligand_chain:
                protein = u.select_atoms(f'protein and not (chainID {ligand_chain})')
            else:
                protein = u.select_atoms('protein')
            
            dssp = MDA_DSSP(protein).run()
            
            # Extract sequence and secondary structure
            sequence = []
            ss_codes = []
            residue_numbers = []
            plddt_values = []
            chain_ids = []  # Track chain IDs
            
            for residue in protein.residues:
                sequence.append(residue.resname)
                residue_numbers.append(residue.resnum)
                chain_id = getattr(residue, 'segid', '')
                if not chain_id:
                    chain_id = getattr(residue, 'chainID', '')
                if not chain_id and hasattr(residue, 'atoms') and len(residue.atoms) > 0:
                    atom_chain = getattr(residue.atoms[0], 'chainID', '')
                    atom_segid = getattr(residue.atoms[0], 'segid', '')
                    chain_id = atom_chain or atom_segid
                chain_ids.append(chain_id if chain_id else '?')  # Track chain robustly
                # Get average B-factor (pLDDT) for this residue
                plddt_values.append(np.mean(residue.atoms.tempfactors))
            
            # Convert 3-letter to 1-letter amino acid codes
            sequence_1letter = []
            for resname in sequence:
                aa_map = {
                    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
                    'MSE': 'M'
                }
                sequence_1letter.append(aa_map.get(resname, 'X'))
            
            # Get DSSP results
            ss_array = dssp.results.dssp
            if len(ss_array) > 0:
                # DSSP array has one entry per residue
                for ss_code in ss_array:
                    ss_codes.append(ss_code if ss_code else 'C')
            else:
                ss_codes = ['C'] * len(sequence_1letter)
            
            print(f"  ✓ MDAnalysis DSSP: {len(sequence_1letter)} residues, first 30 SS: {''.join(ss_codes[:30])}...")
            return {
                'sequence': ''.join(sequence_1letter),
                'ss': ss_codes,
                'residue_numbers': residue_numbers,
                'plddt': plddt_values,
                'chain_ids': chain_ids
            }
        
        except Exception as mda_e:
            print(f"Warning: MDAnalysis DSSP failed for {pdb_path}: {mda_e}")
            print(f"  Using sequence-only fallback...")
            return extract_sequence_and_plddt_fallback(pdb_path, chain_id, ligand_chain)


def extract_sequence_and_plddt_fallback(pdb_path, chain_id='A', ligand_chain=None):
    """
    Fallback: Extract sequence and pLDDT without DSSP.
    Supports multi-chain receptors via ligand_chain parameter.
    """
    from Bio.PDB import PDBParser, MMCIFParser
    
    # Use appropriate parser based on file extension
    if str(pdb_path).lower().endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    
    structure = parser.get_structure('protein', pdb_path)
    model = structure[0]
    
    # Three-letter to one-letter code
    aa_codes = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'MSE': 'M'
    }
    
    sequence = []
    ss_codes = []
    residue_numbers = []
    plddt_values = []
    chain_ids = []  # Track chain IDs
    
    # Determine which chains to process
    if ligand_chain:
        # Multi-chain receptor: process all PROTEIN chains except ligand
        # Filter out non-protein chains (ions, waters, small molecules)
        all_chains = list(model.get_chains())
        chains_to_process = []
        for ch in all_chains:
            if ch.id == ligand_chain:
                continue
            # Check if chain has protein residues (standard amino acids)
            has_protein = any(res.resname.strip() in aa_codes and res.id[0] == ' ' for res in ch)
            if has_protein:
                chains_to_process.append(ch)
    elif isinstance(chain_id, list):
        # Multiple specific chains
        chains_to_process = [model[ch] for ch in chain_id if ch in model]
    else:
        # Single chain
        chains_to_process = [model[chain_id]] if chain_id in model else []
    
    for chain in chains_to_process:
        for residue in chain:
            if residue.id[0] == ' ':  # Standard residue
                resname = residue.resname.strip()
                if resname in aa_codes:
                    sequence.append(aa_codes[resname])
                    ss_codes.append('-')  # No SS info
                    residue_numbers.append(residue.id[1])
                    chain_ids.append(chain.id)  # Track chain ID
                    
                    # Get B-factor
                    bfactors = [atom.get_bfactor() for atom in residue.get_atoms()]
                    plddt_values.append(np.mean(bfactors) if bfactors else 0.0)
    
    return {
        'sequence': ''.join(sequence),
        'ss': ss_codes,
        'residue_numbers': residue_numbers,
        'plddt': plddt_values,
        'chain_ids': chain_ids
    }


def create_alignment_visualization_medoids(
    pdb_files,
    labels,
    reference_pdb=None,
    reference_sequence=None,
    cluster_summary=None,
    ligand_chain='B',
    reference_is_alphafold=False,
    protein_chains=None,
):
    """
    Create sequence alignment visualization comparing cluster medoids.
    
    Mode 1: Compare all medoids to reference (or first medoid if no reference).
    
    Parameters:
    -----------
    pdb_files : list
        List of PDB file paths
    labels : array-like
        Cluster labels
    reference_pdb : str, optional
        Reference PDB structure (if provided, use as comparison baseline)
    reference_sequence : str, optional
        Reference full-length sequence
    cluster_summary : pd.DataFrame, optional
        Cluster summary with consensus residues
    protein_chains : list of str, optional
        Receptor protein chain labels to visualize.
        If provided, these chains are used directly.
        If None, fallback behavior uses ligand_chain exclusion.
    reference_is_alphafold : bool, optional
        If True, treat reference B-factor as pLDDT in the reference lane.
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    from core.io_utils import parse_alphafold_filename, align_sequences
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import multiprocessing
    
    # Find medoid for each cluster
    unique_clusters = sorted(set(labels))
    unique_clusters = [c for c in unique_clusters if c != -1]  # Exclude noise
    
    if len(unique_clusters) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No clusters found",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    # Create intermediate data table
    all_data = []
    
    # Add reference sequence first if available
    print(f"\n=== Reference Sequence Check ===")
    print(f"  reference_pdb provided: {reference_pdb is not None} (value={reference_pdb})")
    print(f"  reference_sequence type: {type(reference_sequence)}")
    print(f"  reference_sequence is None: {reference_sequence is None}")
    print(f"  reference_sequence bool: {bool(reference_sequence) if reference_sequence is not None else 'N/A'}")
    
    if reference_sequence and len(reference_sequence) > 0:
        print(f"  Adding reference sequence: length={len(reference_sequence)}")
        print(f"  First 50 chars: {reference_sequence[:50]}")
        
        # If reference PDB is provided, extract its DSSP data
        ref_ss = None
        ref_plddt = None
        ref_dssp_data = None
        if reference_pdb and os.path.exists(reference_pdb):
            print(f"  Extracting DSSP from reference PDB: {reference_pdb}")
            try:
                if protein_chains:
                    ref_dssp_data = extract_secondary_structure_dssp(reference_pdb, chain_id=protein_chains)
                else:
                    ref_dssp_data = extract_secondary_structure_dssp(reference_pdb, ligand_chain=ligand_chain)
                if ref_dssp_data:
                    ref_ss = ref_dssp_data['ss']
                    if reference_is_alphafold:
                        ref_plddt = ref_dssp_data['plddt']
                    else:
                        ref_plddt = None
                        print("  Reference marked non-AlphaFold: hiding pLDDT lane")
                    print(f"  Reference DSSP extracted: ss_len={len(ref_ss)}, plddt_len={len(ref_plddt) if ref_plddt else 0}")
                    # Check if this is a crystal structure (B-factors likely not pLDDT)
                    if ref_plddt and len(ref_plddt) > 0:
                        avg_bfactor = sum(ref_plddt) / len(ref_plddt)
                        if avg_bfactor < 50:  # Crystal structures typically have B-factors 10-50
                            print(f"  Note: Reference appears to be crystal structure (avg B-factor={avg_bfactor:.1f})")
                            print(f"       B-factors shown are thermal factors, not pLDDT confidence")
            except Exception as e:
                print(f"  Could not extract reference DSSP: {e}")
                print(f"  Reference will show sequence only (no SS or B-factor data)")
        
        all_data.append({
            'name': 'Reference',
            'sequence': reference_sequence,
            'ss': ref_ss,
            'plddt': ref_plddt,
            'chain_ids': ref_dssp_data.get('chain_ids') if ref_dssp_data else None,
            'offset': 0
        })
    else:
        print(f"  No reference sequence to add")
    
    # Extract data for each medoid (in parallel for speed)
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    def extract_medoid_data(cluster_id, pdb_path, reference_seq, lig_ch):
        """Extract DSSP data for a single medoid"""
        # Get combined data for all chains
        if protein_chains:
            dssp_data = extract_secondary_structure_dssp(pdb_path, chain_id=protein_chains)
        else:
            dssp_data = extract_secondary_structure_dssp(pdb_path, ligand_chain=lig_ch)
        filename_info = parse_alphafold_filename(os.path.basename(pdb_path))
        
        # Calculate offset
        aligned_offset = 0
        if reference_seq:
            try:
                alignment = align_sequences(reference_seq, dssp_data['sequence'])
                if alignment and 'target_start' in alignment:
                    aligned_offset = alignment['target_start']
                else:
                    aligned_offset = filename_info['offset'] if filename_info else 0
            except Exception as e:
                aligned_offset = filename_info['offset'] if filename_info else 0
        else:
            aligned_offset = filename_info['offset'] if filename_info else 0
        
        return {
            'cluster_id': cluster_id,
            'name': f'Cluster {cluster_id} Medoid',
            'sequence': dssp_data['sequence'],
            'ss': dssp_data['ss'],
            'plddt': dssp_data['plddt'],
            'chain_ids': dssp_data.get('chain_ids'),
            'offset': aligned_offset
        }
    
    # Collect medoid data - use parallel processing with ThreadPoolExecutor
    medoid_jobs = []
    for cluster_id in unique_clusters:
        cluster_indices = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
        if len(cluster_indices) == 0:
            continue
        
        # Use first structure as medoid
        medoid_idx = cluster_indices[0]
        pdb_path = pdb_files[medoid_idx]
        medoid_jobs.append((cluster_id, pdb_path))
    
    # Process medoids in parallel using ThreadPoolExecutor
    print(f"\n=== Extracting DSSP for {len(medoid_jobs)} medoids (parallel) ===")
    medoid_results = []
    
    # Use ThreadPoolExecutor for parallel DSSP extraction (avoids pickling issues)
    max_workers = min(multiprocessing.cpu_count(), len(medoid_jobs))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_cluster = {
            executor.submit(extract_medoid_data, cluster_id, pdb_path, reference_sequence, ligand_chain): cluster_id
            for cluster_id, pdb_path in medoid_jobs
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_cluster):
            cluster_id = future_to_cluster[future]
            try:
                result = future.result()  # Now returns a single dict again
                medoid_results.append(result)
                print(f"Cluster {cluster_id}: seq_len={len(result['sequence'])}, "
                      f"ss_len={len(result['ss']) if result['ss'] else 0}, "
                      f"offset={result['offset']}")
            except Exception as e:
                print(f"Error extracting data for cluster {cluster_id}: {e}")
                import traceback
                traceback.print_exc()
    
    # Sort by cluster ID and add to all_data
    medoid_results.sort(key=lambda x: x['cluster_id'])
    for result in medoid_results:
        all_data.append(result)
    
    # Debug: Show intermediate data table head
    print(f"\n=== Intermediate Data Table ===")
    print(f"  Total entries: {len(all_data)}")
    if len(all_data) > 0:
        for idx, entry in enumerate(all_data[:3]):  # Show first 3 entries
            print(f"  Entry {idx}: {entry['name']}")
            print(f"    seq length: {len(entry['sequence'])}") 
            print(f"    ss: {type(entry['ss'])}, length={len(entry['ss']) if entry['ss'] else 'None'}")
            if entry['ss']:
                print(f"    First 10 SS codes: {entry['ss'][:10]}")
            print(f"    plddt: {type(entry['plddt'])}, length={len(entry['plddt']) if entry['plddt'] else 'None'}")
            print(f"    offset: {entry['offset']}")
    
    # Create figure with subplots
    n_rows = len(all_data)
    
    # Determine the maximum span for alignment (from reference or all sequences)
    max_span = 0
    reference_length = 0
    if len(all_data) > 0 and all_data[0]['name'] == 'Reference':
        reference_length = len(all_data[0]['sequence'])
        max_span = reference_length
    
    # Also check all medoid spans to ensure we capture everything
    for data in all_data:
        seq_end = data['offset'] + len(data['sequence'])
        max_span = max(max_span, seq_end)
    
    print(f"\n=== Alignment Parameters ===")
    print(f"  Reference length: {reference_length}")
    print(f"  Max span: {max_span}")
    
    # Extract consensus residues per cluster from summary
    consensus_by_cluster = {}
    print(f"\n=== CONSENSUS RESIDUE EXTRACTION ===")
    print(f"  cluster_summary is None: {cluster_summary is None}")
    if cluster_summary is not None:
        print(f"  cluster_summary columns: {list(cluster_summary.columns)}")
        print(f"  cluster_summary shape: {cluster_summary.shape}")
        print(f"  First few rows:\n{cluster_summary.head()}")
        
    if cluster_summary is not None and 'consensus_residues' in cluster_summary.columns:
        for _, row in cluster_summary.iterrows():
            cluster_id = row['cluster']
            if cluster_id != -1 and 'consensus_residues' in row:
                consensus_str = str(row['consensus_residues'])
                print(f"  Cluster {cluster_id} raw consensus string: '{consensus_str}'")
                if consensus_str and consensus_str != 'nan':
                    # Parse comma-separated residue identifiers (format: A:ASN561)
                    # Extract just the residue number from each entry
                    consensus_residues = []
                    for r in consensus_str.split(','):
                        r = r.strip()
                        if ':' in r:
                            # Format: A:ASN561 -> extract 561
                            res_part = r.split(':')[1]  # Get 'ASN561'
                            # Extract digits from 'ASN561' -> '561'
                            res_num_str = ''.join(filter(str.isdigit, res_part))
                            if res_num_str:
                                consensus_residues.append(int(res_num_str))
                    
                    if consensus_residues:
                        consensus_by_cluster[cluster_id] = consensus_residues
                        print(f"  Cluster {cluster_id} consensus ({len(consensus_residues)} residues): {consensus_residues[:10]}{'...' if len(consensus_residues) > 10 else ''}")
    
    print(f"  Total clusters with consensus: {len(consensus_by_cluster)}")
    print(f"===================================\n")
    
    # Plotly constraint: vertical_spacing <= 1 / (rows - 1)
    # Clamp spacing dynamically for large numbers of rows.
    if n_rows > 1:
        max_allowed_spacing = (1.0 / (n_rows - 1)) - 1e-6
        # Keep spacing small for many rows; large spacing makes lanes unusable
        base_spacing = 0.012 if n_rows >= 6 else 0.03
        vertical_spacing = min(base_spacing, max_allowed_spacing)
    else:
        vertical_spacing = 0.0

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        subplot_titles=[d['name'] for d in all_data],
        vertical_spacing=vertical_spacing,
        row_heights=[1.0] * n_rows
    )
    
    # Add tracks for each entry in data table
    print(f"\n=== Generating track data in parallel for {len(all_data)} rows ===")
    
    # Prepare arguments for parallel execution
    track_jobs = []
    for i, data in enumerate(all_data):
        consensus_residues = None
        if 'cluster_id' in data and data['cluster_id'] in consensus_by_cluster:
            consensus_residues = consensus_by_cluster[data['cluster_id']]
        
        track_jobs.append({
            'sequence': data['sequence'],
            'ss_codes': data['ss'],
            'plddt_values': data['plddt'],
            'chain_ids': data.get('chain_ids'),  # May be None for reference
            'row': i + 1,
            'offset': data['offset'],
            'reference_length': reference_length,
            'max_span': max_span,
            'consensus_residues': consensus_residues
        })
    
    # Generate all track data in parallel
    track_results = []
    max_workers = min(multiprocessing.cpu_count(), len(track_jobs))
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                generate_sequence_track_data,
                job['sequence'], job['ss_codes'], job['plddt_values'],
                job['chain_ids'], job['row'], job['offset'], job['reference_length'],
                job['max_span'], job['consensus_residues']
            ): job['row'] for job in track_jobs
        }
        
        for future in as_completed(futures):
            row_num = futures[future]
            try:
                result = future.result()
                track_results.append(result)
                print(f"  ✓ Row {row_num} complete: {len(result['shapes'])} shapes")
            except Exception as e:
                print(f"  ✗ Row {row_num} failed: {e}")
                import traceback
                traceback.print_exc()
    
    # Sort results by row number
    track_results.sort(key=lambda x: x['row'])
    
    # Combine all shapes, annotations, and traces
    print(f"\n=== Combining results and updating figure ===")
    all_shapes = []
    all_annotations = []
    
    for result in track_results:
        all_shapes.extend(result['shapes'])
        all_annotations.extend(result['annotations'])
        if result['trace']:
            fig.add_trace(result['trace'], row=result['row'], col=1)
    
    # Update figure with all shapes and annotations at once
    print(f"  Adding {len(all_shapes)} shapes and {len(all_annotations)} annotations...")
    if all_shapes:
        fig.update_layout(shapes=tuple(all_shapes))
    if all_annotations:
        fig.update_layout(annotations=tuple(all_annotations))
    
    # Update all axes
    for i in range(1, len(all_data) + 1):
        is_last_row = i == len(all_data)
        fig.update_xaxes(
            title_text="Residue Position" if is_last_row else None,
            row=i,
            col=1,
            showticklabels=is_last_row,
            tickmode='auto'
        )
        fig.update_yaxes(
            showticklabels=False, range=[-0.1, 1.1], row=i, col=1
        )
    
    fig.update_layout(
        title='Cluster Medoid Sequence Alignment',
        height=max(500, 95 * n_rows + 120),
        showlegend=False  # No legend needed
    )
    
    # Set x-axis range and y-axis for all subplots
    for i in range(1, n_rows + 1):
        xaxis_key = 'xaxis' if i == 1 else f'xaxis{i}'
        yaxis_key = 'yaxis' if i == 1 else f'yaxis{i}'
        fig.update_layout({
            xaxis_key: dict(
                range=[0.5, reference_length + 0.5],
                constrain='domain'
            ),
            yaxis_key: dict(
                range=[-0.1, 1.1],
                fixedrange=False,
                showticklabels=False
            )
        })
    
    return fig


def create_alignment_visualization_within_cluster(
    pdb_files,
    labels,
    cluster_id,
    reference_pdb=None,
    protein_chains=None,
    ligand_chain='B'
):
    """
    Create sequence alignment visualization for all structures within a cluster.
    
    Mode 2: Compare all structures in a cluster to the cluster medoid.
    
    Parameters:
    -----------
    pdb_files : list
        List of PDB file paths
    labels : array-like
        Cluster labels
    cluster_id : int
        Cluster to visualize
    reference_pdb : str, optional
        Reference PDB structure
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    from core.io_utils import parse_alphafold_filename
    
    # Get structures in this cluster
    cluster_indices = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
    
    if len(cluster_indices) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text=f"No structures in cluster {cluster_id}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    # Extract data for all structures
    structure_data = []
    
    for idx in cluster_indices:
        pdb_path = pdb_files[idx]
        if protein_chains:
            dssp_data = extract_secondary_structure_dssp(pdb_path, chain_id=protein_chains)
        else:
            dssp_data = extract_secondary_structure_dssp(pdb_path, ligand_chain=ligand_chain)
        filename_info = parse_alphafold_filename(os.path.basename(pdb_path))
        
        structure_data.append({
            'index': idx,
            'pdb_path': pdb_path,
            'filename': os.path.basename(pdb_path),
            'sequence': dssp_data['sequence'],
            'ss': dssp_data['ss'],
            'residue_numbers': dssp_data['residue_numbers'],
            'plddt': dssp_data['plddt'],
            'chain_ids': dssp_data.get('chain_ids'),
            'filename_info': filename_info
        })
    
    # Create figure
    n_structures = len(structure_data)
    # Plotly constraint: vertical_spacing <= 1 / (rows - 1)
    if n_structures > 1:
        max_allowed_spacing = (1.0 / (n_structures - 1)) - 1e-6
        base_spacing = 0.012 if n_structures >= 6 else 0.03
        vertical_spacing = min(base_spacing, max_allowed_spacing)
    else:
        vertical_spacing = 0.0

    fig = make_subplots(
        rows=n_structures,
        cols=1,
        subplot_titles=[f'{s["filename"][:50]}...' if len(s["filename"]) > 50 else s["filename"] 
                       for s in structure_data],
        vertical_spacing=vertical_spacing,
        row_heights=[1.0] * n_structures
    )
    
    # Add tracks
    for i, struct in enumerate(structure_data):
        add_sequence_track(
            fig,
            struct['sequence'],
            struct['ss'],
            struct['plddt'],
            row=i + 1,
            offset=struct['filename_info']['offset'] if struct['filename_info'] else 0,
            reference_length=0,  # No alignment for within-cluster view
            max_span=0,
            consensus_residues=None,  # No consensus markers in within-cluster view
            chain_ids=struct.get('chain_ids')
        )
    
    fig.update_layout(
        title=f'Cluster {cluster_id} - All Structures',
        height=max(420, 90 * n_structures + 100),
        showlegend=True
    )
    
    return fig


def generate_sequence_track_data(sequence, ss_codes, plddt_values, chain_ids, row, offset=0, reference_length=0, max_span=0, consensus_residues=None):
    """
    Generate shape/annotation data for a sequence track (parallelizable).
    
    Returns a dictionary with shapes, annotations, and trace data.
    This function is pure and can be run in parallel for multiple rows.
    
    Parameters:
    -----------
    sequence : str
        Amino acid sequence
    ss_codes : list of str or None
        Secondary structure codes (DSSP format)
    plddt_values : list of float or None
        pLDDT confidence scores (0-100)
    chain_ids : list of str or None
        Chain IDs for each residue
    row : int
        Subplot row number
    offset : int
        Residue numbering offset
    reference_length : int
        Length of reference sequence (for gap filling)
    max_span : int
        Maximum span across all sequences (for x-axis range)
    consensus_residues : list of int, optional
        List of consensus contact residue numbers to highlight
    
    Returns:
    --------
    dict : {
        'shapes': list of shape dicts,
        'annotations': list of annotation dicts,
        'trace': plotly trace object,
        'row': int
    }
    """
    print(f"  generate_track_data: row={row}, seq_len={len(sequence) if sequence else 0}, "
          f"offset={offset}, ref_len={reference_length}, consensus={len(consensus_residues) if consensus_residues else 0}")
    
    if not sequence or len(sequence) == 0:
        print(f"  WARNING: Empty sequence for row {row}")
        return {'shapes': [], 'annotations': [], 'trace': None, 'row': row}
    
    seq_len = len(sequence)
    seq_start = offset + 1
    seq_end = offset + seq_len
    # Total length should accommodate both reference AND the current sequence
    total_length = max(reference_length, seq_end) if reference_length > 0 else seq_end
    
    # Detect whether this row has meaningful pLDDT values.
    # If not, collapse the pLDDT lane and show an inline note.
    has_plddt_data = False
    if plddt_values is not None and len(plddt_values) > 0:
        has_plddt_data = any(v is not None for v in plddt_values)

    # Dynamic lane geometry
    chain_y0, chain_y1 = 0.9, 1.0
    if has_plddt_data:
        plddt_y0, plddt_y1 = 0.6, 0.9
        ss_y0, ss_y1 = 0.33, 0.6
        aa_y = 0.465
    else:
        # Keep SS lane size consistent with rows that have pLDDT
        # (pLDDT lane is omitted and replaced by note)
        plddt_y0, plddt_y1 = None, None
        ss_y0, ss_y1 = 0.33, 0.6
        aa_y = 0.465

    # === BUILD DATA ARRAYS ===
    plddt_array = np.full(total_length, -10.0)
    ss_array = np.full(total_length, -1.0)
    
    # SS code to numeric mapping (DSSP codes) - using distinct integers for distinct colors
    # H = alpha helix, G = 3-10 helix, I = pi helix
    # E = extended strand (beta sheet), B = beta bridge
    # T = turn, S = bend, C = coil, - = no structure
    ss_map = {
        'H': 1,   # Alpha helix → Red
        'G': 2,   # 3-10 helix → Orange
        'I': 3,   # Pi helix → Pink
        'E': 4,   # Extended strand (beta sheet) → Blue
        'B': 5,   # Beta bridge → Light blue
        'T': 6,   # Turn → Green
        'S': 7,   # Bend → Yellow
        'C': 0,   # Coil → Gray
        '-': 0    # No structure → Gray
    }
    
    # Fill in actual sequence data
    for i in range(seq_len):
        pos_idx = seq_start + i - 1
        
        # Bounds check: skip if position is out of range
        if pos_idx < 0 or pos_idx >= total_length:
            continue
        
        if has_plddt_data and i < len(plddt_values):
            plddt_array[pos_idx] = plddt_values[i]
        
        if ss_codes is not None and i < len(ss_codes):
            ss_code = str(ss_codes[i]).strip().upper()
            ss_array[pos_idx] = ss_map.get(ss_code, 0)
        else:
            ss_array[pos_idx] = 0
    
    x_positions = np.arange(1, total_length + 1)

    # Plotly axis refs: first subplot uses "x"/"y" (not "x1"/"y1")
    xref = 'x' if row == 1 else f'x{row}'
    yref = 'y' if row == 1 else f'y{row}'
    
    # === GENERATE SHAPES (rectangles) ===
    shapes = []
    
    # Chain ID lane (top: 0.9-1.0) - only show where we have actual sequence data
    chain_color_map = {
        'A': 'rgba(100, 150, 255, 0.7)', 
        'B': 'rgba(255, 150, 100, 0.7)',
        'C': 'rgba(100, 255, 150, 0.7)',
        'D': 'rgba(255, 100, 255, 0.7)'
    }
    if chain_ids:
        for i in range(seq_len):
            pos_idx = seq_start + i - 1
            if pos_idx >= 0 and pos_idx < total_length and i < len(chain_ids):
                chain_id = chain_ids[i]
                if chain_id:
                    color = chain_color_map.get(chain_id, 'rgba(150, 150, 150, 0.5)')
                    x_pos = pos_idx + 1  # Convert to 1-based
                    shapes.append(dict(
                        type="rect", xref=xref, yref=yref,
                        x0=x_pos - 0.5, y0=chain_y0, x1=x_pos + 0.5, y1=chain_y1,
                        fillcolor=color, line=dict(width=0)
                    ))
    
    # pLDDT rectangles (middle-top lane) when available
    if has_plddt_data:
        for x_pos, val in zip(x_positions, plddt_array):
            if val <= -5:
                color = 'rgba(0, 0, 0, 0.0)'
            elif val < 50:
                color = 'rgba(255, 125, 69, 0.8)'
            elif val < 70:
                color = 'rgba(255, 219, 19, 0.8)'
            elif val < 90:
                color = 'rgba(101, 203, 243, 0.8)'
            else:
                color = 'rgba(0, 83, 214, 0.8)'

            shapes.append(dict(
                type="rect", xref=xref, yref=yref,
                x0=x_pos - 0.5, y0=plddt_y0, x1=x_pos + 0.5, y1=plddt_y1,
                fillcolor=color, line=dict(width=0)
            ))
    
    # SS rectangles (middle: 0.33-0.6)
    ss_color_lookup = {
        -1: 'rgba(180, 180, 180, 0.5)',  # Gap
        0: 'rgba(200, 200, 200, 0.5)',   # Coil/no structure (C, -)
        1: 'rgba(255, 0, 0, 0.8)',       # Alpha helix (H) - Red
        2: 'rgba(255, 140, 0, 0.8)',     # 3-10 helix (G) - Orange
        3: 'rgba(255, 150, 150, 0.8)',   # Pi helix (I) - Pink
        4: 'rgba(0, 0, 255, 0.8)',       # Beta sheet (E) - Blue
        5: 'rgba(100, 100, 255, 0.8)',   # Beta bridge (B) - Light blue
        6: 'rgba(0, 200, 0, 0.8)',       # Turn (T) - Green
        7: 'rgba(200, 200, 0, 0.8)'      # Bend (S) - Yellow
    }
    
    for x_pos, val in zip(x_positions, ss_array):
        color = ss_color_lookup.get(int(val), 'rgba(200, 200, 200, 0.5)')
        shapes.append(dict(
            type="rect", xref=xref, yref=yref,
            x0=x_pos - 0.5, y0=ss_y0, x1=x_pos + 0.5, y1=ss_y1,
            fillcolor=color, line=dict(width=0)
        ))
    
    # === ANNOTATIONS (amino acids every 10 residues) ===
    annotations = []
    positions = np.arange(seq_start, seq_end + 1)
    label_stride = 20 if seq_len > 120 else 10
    for i, (pos, aa) in enumerate(zip(positions, sequence)):
        if i % label_stride == 0:
            annotations.append(dict(
                x=pos, y=aa_y, xref=xref, yref=yref,
                text=aa, showarrow=False,
                font=dict(size=7, color='white', family='monospace')
            ))

    # Inline note for structures without pLDDT lane (e.g., non-AlphaFold reference)
    if not has_plddt_data:
        annotations.append(dict(
            x=seq_start,
            y=1.04,
            xref=xref,
            yref=yref,
            xanchor='left',
            text='Not an AlphaFold predicted model - No pLDDT data to display',
            showarrow=False,
            font=dict(size=10, color='rgba(70,70,70,0.9)')
        ))
    
    # === CONSENSUS MARKERS (bottom third: 0.0-0.33) ===
    if consensus_residues:
        print(f"    🎯 Adding {len(consensus_residues)} consensus markers to row {row}")
        for res_num in consensus_residues:
            # Single-style consensus marker: small magenta tick at very bottom to minimize gap
            shapes.append(dict(
                type="line", xref=xref, yref=yref,
                x0=res_num, y0=0.1, x1=res_num, y1=0.35,
                line=dict(color='rgba(214, 51, 132, 0.85)', width=1.8)
            ))
    
    # === HOVER DATA (scatter trace) ===
    hover_text = []
    for i, aa in enumerate(sequence):
        text_parts = [f"AA: {aa}", f"Pos: {int(positions[i])}"]
        if chain_ids and i < len(chain_ids):
            text_parts.append(f"Chain: {chain_ids[i]}")
        if plddt_values is not None and i < len(plddt_values):
            text_parts.append(f"pLDDT: {plddt_values[i]:.1f}")
        if ss_codes is not None and i < len(ss_codes):
            ss = str(ss_codes[i]).strip().upper()
            text_parts.append(f"SS: {ss}")
        hover_text.append("<br>".join(text_parts))
    
    trace = go.Scattergl(
        x=positions, y=[0.5] * len(positions),
        mode='markers', marker=dict(size=10, color='rgba(0,0,0,0)'),
        showlegend=False, hovertemplate='%{text}<extra></extra>',
        text=hover_text
    )
    
    print(f"  Track data generated: {len(shapes)} shapes + {len(annotations)} annotations")
    
    return {
        'shapes': shapes,
        'annotations': annotations,
        'trace': trace,
        'row': row
    }


def add_sequence_track(fig, sequence, ss_codes, plddt_values, row, offset=0, reference_length=0, max_span=0, consensus_residues=None, chain_ids=None):
    """
    Add a sequence track to the figure (wrapper for backward compatibility).
    
    This now just calls generate_sequence_track_data and adds results to figure.
    """
    track_data = generate_sequence_track_data(
        sequence, ss_codes, plddt_values, chain_ids, row, offset,
        reference_length, max_span, consensus_residues
    )
    
    # Add shapes and annotations
    if track_data['shapes']:
        fig.update_layout(shapes=fig.layout.shapes + tuple(track_data['shapes']))
    if track_data['annotations']:
        fig.update_layout(annotations=fig.layout.annotations + tuple(track_data['annotations']))
    
    # Add trace
    if track_data['trace']:
        fig.add_trace(track_data['trace'], row=row, col=1)
    
    # Update axes
    fig.update_xaxes(
        title_text="Residue Position", row=row, col=1,
        showticklabels=True, tickmode='auto'
    )
    fig.update_yaxes(
        showticklabels=False, range=[-0.1, 1.1], row=row, col=1
    )

