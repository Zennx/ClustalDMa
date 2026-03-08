"""
Sequence alignment visualization with secondary structure and pLDDT confidence
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
from pathlib import Path


def extract_secondary_structure_dssp(pdb_path, chain_id='A'):
    """
    Extract secondary structure using MDAnalysis DSSP.
    
    Parameters:
    -----------
    pdb_path : str
        Path to PDB file
    chain_id : str
        Chain ID to analyze
    
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
        import MDAnalysis as mda
        from MDAnalysis.analysis.dssp import DSSP as MDA_DSSP
        
        # Load structure
        u = mda.Universe(pdb_path)
        
        # Select protein chain
        protein = u.select_atoms(f'protein and segid {chain_id}')
        if len(protein) == 0:
            # Try without segid (some PDBs use chainID)
            protein = u.select_atoms(f'protein and chainID {chain_id}')
        
        if len(protein) == 0:
            print(f"Warning: Chain {chain_id} not found, using all protein atoms")
            protein = u.select_atoms('protein')
        
        # Run DSSP
        dssp = MDA_DSSP(protein).run()
        
        # Extract data per residue
        sequence = []
        ss_codes = []
        residue_numbers = []
        plddt_values = []
        
        # DSSP results can be 2D array with shape (1, n_residues) or 1D with shape (n_residues,)
        # Flatten to 1D first
        dssp_array = dssp.results.dssp.flatten() if hasattr(dssp.results.dssp, 'flatten') else dssp.results.dssp
        
        # Convert DSSP array to string (MDAnalysis returns array of characters)
        # Join all characters into a single string for easy indexing
        dssp_string = "".join(dssp_array)
        print(f"  Extracted DSSP: {len(dssp_string)} residues, first 30: '{dssp_string[:30]}...'")
        print(f"  DSSP string length: {len(dssp_string)}")
        
        for res_idx, residue in enumerate(protein.residues):
            # Get amino acid (1-letter code)
            aa = residue.resname
            aa_map = {
                'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
                'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
            }
            aa_letter = aa_map.get(aa, 'X')
            
            # Get secondary structure from DSSP string
            if res_idx < len(dssp_string):
                ss = dssp_string[res_idx]
                if not ss or ss == '':
                    ss = '-'
            else:
                ss = '-'
            
            # Get B-factor (pLDDT) - average over all atoms in residue
            bfactors = residue.atoms.tempfactors
            avg_bfactor = np.mean(bfactors) if len(bfactors) > 0 else 0.0
            
            sequence.append(aa_letter)
            ss_codes.append(ss)
            residue_numbers.append(residue.resnum)
            plddt_values.append(avg_bfactor)
        
        return {
            'sequence': ''.join(sequence),
            'ss': ss_codes,
            'residue_numbers': residue_numbers,
            'plddt': plddt_values
        }
    
    except Exception as e:
        print(f"Warning: MDAnalysis DSSP failed for {pdb_path}: {e}")
        import traceback
        traceback.print_exc()
        # Fallback: extract sequence and B-factors without DSSP
        return extract_sequence_and_plddt_fallback(pdb_path, chain_id)


def extract_sequence_and_plddt_fallback(pdb_path, chain_id='A'):
    """
    Fallback: Extract sequence and pLDDT without DSSP.
    """
    from Bio.PDB import PDBParser
    
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
    
    if chain_id in model:
        for residue in model[chain_id]:
            if residue.id[0] == ' ':  # Standard residue
                resname = residue.resname.strip()
                if resname in aa_codes:
                    sequence.append(aa_codes[resname])
                    ss_codes.append('-')  # No SS info
                    residue_numbers.append(residue.id[1])
                    
                    # Get B-factor
                    bfactors = [atom.get_bfactor() for atom in residue.get_atoms()]
                    plddt_values.append(np.mean(bfactors) if bfactors else 0.0)
    
    return {
        'sequence': ''.join(sequence),
        'ss': ss_codes,
        'residue_numbers': residue_numbers,
        'plddt': plddt_values
    }


def create_alignment_visualization_medoids(pdb_files, labels, reference_pdb=None, reference_sequence=None, cluster_summary=None):
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
        if reference_pdb and os.path.exists(reference_pdb):
            print(f"  Extracting DSSP from reference PDB: {reference_pdb}")
            try:
                ref_dssp_data = extract_secondary_structure_dssp(reference_pdb, chain_id='A')
                if ref_dssp_data:
                    ref_ss = ref_dssp_data['ss']
                    ref_plddt = ref_dssp_data['plddt']
                    print(f"  Reference DSSP extracted: ss_len={len(ref_ss)}, plddt_len={len(ref_plddt)}")
            except Exception as e:
                print(f"  Could not extract reference DSSP: {e}")
        
        all_data.append({
            'name': 'Reference',
            'sequence': reference_sequence,
            'ss': ref_ss,
            'plddt': ref_plddt,
            'offset': 0
        })
    else:
        print(f"  No reference sequence to add")
    
    # Extract data for each medoid (in parallel for speed)
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    def extract_medoid_data(cluster_id, pdb_path, reference_seq):
        """Extract DSSP data for a single medoid"""
        dssp_data = extract_secondary_structure_dssp(pdb_path, chain_id='A')
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
            executor.submit(extract_medoid_data, cluster_id, pdb_path, reference_sequence): cluster_id
            for cluster_id, pdb_path in medoid_jobs
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_cluster):
            cluster_id = future_to_cluster[future]
            try:
                result = future.result()
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
    
    fig = make_subplots(
        rows=n_rows,
        cols=1,
        subplot_titles=[d['name'] for d in all_data],
        vertical_spacing=0.12,  # Increased spacing to prevent overlap
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
                job['row'], job['offset'], job['reference_length'],
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
        fig.update_xaxes(
            title_text="Residue Position", row=i, col=1,
            showticklabels=True, tickmode='auto'
        )
        fig.update_yaxes(
            showticklabels=False, range=[-0.1, 1.1], row=i, col=1
        )
    
    fig.update_layout(
        title='Cluster Medoid Sequence Alignment',
        height=200 * n_rows,  # Increased height per row for better spacing
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
def create_alignment_visualization_within_cluster(pdb_files, labels, cluster_id, reference_pdb=None):
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
        dssp_data = extract_secondary_structure_dssp(pdb_path, chain_id='A')
        filename_info = parse_alphafold_filename(os.path.basename(pdb_path))
        
        structure_data.append({
            'index': idx,
            'pdb_path': pdb_path,
            'filename': os.path.basename(pdb_path),
            'sequence': dssp_data['sequence'],
            'ss': dssp_data['ss'],
            'residue_numbers': dssp_data['residue_numbers'],
            'plddt': dssp_data['plddt'],
            'filename_info': filename_info
        })
    
    # Create figure
    n_structures = len(structure_data)
    fig = make_subplots(
        rows=n_structures,
        cols=1,
        subplot_titles=[f'{s["filename"][:50]}...' if len(s["filename"]) > 50 else s["filename"] 
                       for s in structure_data],
        vertical_spacing=0.01,
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
            consensus_residues=None  # No consensus markers in within-cluster view
        )
    
    fig.update_layout(
        title=f'Cluster {cluster_id} - All Structures',
        height=100 * n_structures,
        showlegend=True
    )
    
    return fig


def generate_sequence_track_data(sequence, ss_codes, plddt_values, row, offset=0, reference_length=0, max_span=0, consensus_residues=None):
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
    total_length = reference_length if reference_length > 0 else seq_end
    
    # === BUILD DATA ARRAYS ===
    plddt_array = np.full(total_length, -10.0)
    ss_array = np.full(total_length, -1.0)
    
    # SS code to numeric mapping
    ss_map = {
        'H': 1, 'G': 2, 'I': 3, 'E': 4, 'B': 5, 'T': 6, 'S': 7, 'C': 0, '-': 0
    }
    
    # Fill in actual sequence data
    for i in range(seq_len):
        pos_idx = seq_start + i - 1
        
        if plddt_values is not None and i < len(plddt_values):
            plddt_array[pos_idx] = plddt_values[i]
        else:
            plddt_array[pos_idx] = 50
        
        if ss_codes is not None and i < len(ss_codes):
            ss_code = str(ss_codes[i]).strip().upper()
            ss_array[pos_idx] = ss_map.get(ss_code, 0)
        else:
            ss_array[pos_idx] = 0
    
    x_positions = np.arange(1, total_length + 1)
    
    # === GENERATE SHAPES (rectangles) ===
    shapes = []
    
    # pLDDT rectangles (top half: 0.5-1.0)
    for x_pos, val in zip(x_positions, plddt_array):
        if val <= -5:
            color = 'rgba(180, 180, 180, 0.5)'
        elif val < 50:
            color = 'rgba(255, 125, 69, 0.8)'
        elif val < 70:
            color = 'rgba(255, 219, 19, 0.8)'
        elif val < 90:
            color = 'rgba(101, 203, 243, 0.8)'
        else:
            color = 'rgba(0, 83, 214, 0.8)'
        
        shapes.append(dict(
            type="rect", xref=f"x{row}", yref=f"y{row}",
            x0=x_pos - 0.5, y0=0.5, x1=x_pos + 0.5, y1=1.0,
            fillcolor=color, line=dict(width=0)
        ))
    
    # SS rectangles (bottom half: 0.0-0.5)
    ss_color_lookup = {
        -1: 'rgba(180, 180, 180, 0.5)', 0: 'rgba(200, 200, 200, 0.5)',
        1: 'rgba(255, 0, 0, 0.8)', 2: 'rgba(255, 100, 0, 0.8)',
        3: 'rgba(255, 150, 150, 0.8)', 4: 'rgba(0, 0, 255, 0.8)',
        5: 'rgba(100, 100, 255, 0.8)', 6: 'rgba(0, 200, 0, 0.8)',
        7: 'rgba(200, 200, 0, 0.8)'
    }
    
    for x_pos, val in zip(x_positions, ss_array):
        color = ss_color_lookup.get(int(val), 'rgba(200, 200, 200, 0.5)')
        shapes.append(dict(
            type="rect", xref=f"x{row}", yref=f"y{row}",
            x0=x_pos - 0.5, y0=0.0, x1=x_pos + 0.5, y1=0.5,
            fillcolor=color, line=dict(width=0)
        ))
    
    # === ANNOTATIONS (amino acids every 10 residues) ===
    annotations = []
    positions = np.arange(seq_start, seq_end + 1)
    for i, (pos, aa) in enumerate(zip(positions, sequence)):
        if i % 10 == 0:
            annotations.append(dict(
                x=pos, y=0.75, xref=f"x{row}", yref=f"y{row}",
                text=aa, showarrow=False,
                font=dict(size=7, color='white', family='monospace')
            ))
    
    # === CONSENSUS MARKERS ===
    if consensus_residues:
        print(f"    🎯 Adding {len(consensus_residues)} consensus markers to row {row}")
        for res_num in consensus_residues:
            shapes.append(dict(
                type="line", xref=f"x{row}", yref=f"y{row}",
                x0=res_num, y0=0.0, x1=res_num, y1=1.0,
                line=dict(color='rgba(255, 0, 255, 0.4)', width=2)
            ))
            annotations.append(dict(
                x=res_num, y=-0.05, xref=f"x{row}", yref=f"y{row}",
                text="▼", showarrow=False,
                font=dict(size=8, color='rgba(255, 0, 255, 0.6)')
            ))
    
    # === HOVER DATA (scatter trace) ===
    hover_text = []
    for i, aa in enumerate(sequence):
        text_parts = [f"AA: {aa}", f"Pos: {int(positions[i])}"]
        if plddt_values is not None and i < len(plddt_values):
            text_parts.append(f"pLDDT: {plddt_values[i]:.1f}")
        if ss_codes is not None and i < len(ss_codes):
            ss = str(ss_codes[i]).strip().upper()
            text_parts.append(f"SS: {ss}")
        hover_text.append("<br>".join(text_parts))
    
    trace = go.Scattergl(
        x=positions, y=[0.5] * len(positions),
        mode='markers', marker=dict(size=0.1, color='rgba(0,0,0,0)'),
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


def add_sequence_track(fig, sequence, ss_codes, plddt_values, row, offset=0, reference_length=0, max_span=0, consensus_residues=None):
    """
    Add a sequence track to the figure (wrapper for backward compatibility).
    
    This now just calls generate_sequence_track_data and adds results to figure.
    """
    track_data = generate_sequence_track_data(
        sequence, ss_codes, plddt_values, row, offset,
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

