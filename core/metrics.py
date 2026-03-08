"""
Distance metric computation for structural clustering
"""

import os
import numpy as np
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from joblib import Parallel, delayed
from scipy.spatial.distance import pdist, squareform


class DistanceMetrics:
    """Collection of distance metric computation methods"""
    
    @staticmethod
    def _process_structure_contacts(pdb_file, distance_cutoff, protein_selection, nucleic_selection, motif_residues=None, residue_offset=0):
        """
        Helper function to process a single structure for contact detection.
        Used for parallel processing.
        
        MOTIF SCREENING:
        - If motif_residues provided: Quick motif contact screening (local search, fast)
        - If no motif: Process all protein-nucleic contacts
        
        RESIDUE OFFSET CORRECTION:
        - For chopped AlphaFold models, residue numbers in PDB don't match full-length positions
        - residue_offset is added to all protein residue numbers to correct this
        
        Parameters:
        -----------
        motif_residues : dict, optional
            Dict mapping chain_id -> list of residue numbers
            Example: {'A': [10, 11, 12], 'B': [30, 31, 32]}
            If provided, structures with no motif contacts are filtered out
        residue_offset : int
            Offset to add to protein residue numbers (for chopped sequences)
            Example: if PDB has residues 1-100 but they're actually 501-600, offset=500
        
        Returns:
        --------
        tuple : (contact_set, contact_list)
            Returns (set(), []) if filtering fails (off-target)
        """
        try:
            # Load structure with MDTraj
            traj = md.load(pdb_file)
            topology = traj.topology
            
            # Select protein and nucleic atoms
            protein_atoms = topology.select(protein_selection)
            nucleic_atoms = topology.select(nucleic_selection)
            
            if len(protein_atoms) == 0 or len(nucleic_atoms) == 0:
                return set(), []
            
            # MOTIF SCREENING (if provided)
            if motif_residues is not None and len(motif_residues) > 0:
                # Build selection for motif residues only
                motif_atom_indices = []
                for chain_id, residue_nums in motif_residues.items():
                    for res in topology.residues:
                        if res.chain.chain_id == chain_id and res.resSeq in residue_nums:
                            motif_atom_indices.extend([atom.index for atom in res.atoms])
                
                if len(motif_atom_indices) > 0:
                    # Quick screening: check if ANY motif atom contacts nucleic acid
                    motif_nucleic_pairs = np.array([(m, n) for m in motif_atom_indices for n in nucleic_atoms])
                    
                    if len(motif_nucleic_pairs) > 0:
                        motif_distances = md.compute_distances(traj, motif_nucleic_pairs)[0]
                        motif_has_contact = np.any(motif_distances <= (distance_cutoff / 10.0))
                        
                        if not motif_has_contact:
                            # OFF-TARGET: Motif makes no contact, skip this pose entirely
                            return set(), []
            
            # FULL CONTACT SEARCH: Process all protein-nucleic contacts
            # Create all possible atom pairs
            pairs = np.array([(p, n) for p in protein_atoms for n in nucleic_atoms])
            
            # Compute distances (vectorized)
            distances = md.compute_distances(traj, pairs)[0]
            
            # Find contacts below cutoff
            contact_mask = distances <= (distance_cutoff / 10.0)  # Å to nm
            contact_atom_pairs = pairs[contact_mask]
            
            # Convert to protein residue contacts
            protein_residues_in_contact = set()
            contact_list = []
            
            # Track minimum distance for each protein residue
            residue_min_distances = {}
            
            for idx, (p_idx, n_idx) in enumerate(contact_atom_pairs):
                p_res = topology.atom(p_idx).residue
                n_res = topology.atom(n_idx).residue
                
                # Apply residue offset for chopped AlphaFold models (only to protein chain A)
                # This corrects residue numbering to match full-length sequence positions
                if p_res.chain.chain_id == 'A':
                    corrected_resnum = p_res.resSeq + residue_offset
                else:
                    corrected_resnum = p_res.resSeq
                
                # CRITICAL: Include chain ID to distinguish different chains and L-R symmetry
                protein_res_id = f"{p_res.chain.chain_id}:{p_res.name}{corrected_resnum}"
                nucleic_res_id = f"{n_res.chain.chain_id}:{n_res.name}{n_res.resSeq}"
                
                # Get distance for this contact (convert nm to Å)
                contact_distance = distances[contact_mask][idx] * 10.0
                
                # Track minimum distance for each protein residue
                if protein_res_id not in residue_min_distances:
                    residue_min_distances[protein_res_id] = contact_distance
                else:
                    residue_min_distances[protein_res_id] = min(residue_min_distances[protein_res_id], contact_distance)
                
                # Add to set of protein residues in contact
                protein_residues_in_contact.add(protein_res_id)
                
                # Store full pair for export with distance
                contact_list.append({
                    'protein_residue': protein_res_id,
                    'nucleic_residue': nucleic_res_id,
                    'protein_res_idx': p_res.index,
                    'nucleic_res_idx': n_res.index,
                    'protein_chain': p_res.chain.chain_id,
                    'nucleic_chain': n_res.chain.chain_id,
                    'distance': contact_distance,  # Store actual distance in Å
                    'min_residue_distance': residue_min_distances[protein_res_id]  # Minimum for this protein residue
                })
            
            return protein_residues_in_contact, contact_list
            
        except Exception as e:
            print(f"  Warning: Error processing {pdb_file}: {e}")
            return set(), []
    
    @staticmethod
    def _compute_rmsd_pair(args):
        """
        Helper function to compute RMSD between two structures.
        Used for parallel processing.
        
        Parameters:
        -----------
        args : tuple
            (i, j, pdb_file_i, pdb_file_j, selection, align_structures)
        
        Returns:
        --------
        tuple : (i, j, rmsd_value)
        """
        i, j, pdb_i, pdb_j, selection, do_align = args
        try:
            # Load universes fresh in each worker (avoids pickling overhead)
            u1 = mda.Universe(pdb_i)
            u2 = mda.Universe(pdb_j)
            
            # Align if requested (but usually we don't for nucleic P atoms)
            if do_align:
                align.AlignTraj(u1, u2, select=selection, in_memory=True).run()
            
            # Compute RMSD
            rmsd_value = rms.rmsd(
                u1.select_atoms(selection).positions,
                u2.select_atoms(selection).positions,
                superposition=False
            )
            return (i, j, rmsd_value)
        except Exception as e:
            return (i, j, np.nan)
    
    @staticmethod
    def compute_rmsd_matrix(universes, selection='protein and name CA', 
                           align_structures=True, reference_idx=0, normalize=False, n_jobs=1):
        """
        Compute pairwise RMSD distance matrix with VECTORIZED computation
        
        Parameters:
        -----------
        universes : list
            List of MDAnalysis Universe objects (or will use pdb_files if available)
        selection : str
            Selection string for RMSD calculation
        align_structures : bool
            Whether to align structures before RMSD
        reference_idx : int
            Index of reference structure for alignment
        normalize : bool
            Whether to normalize to [0, 1]
        n_jobs : int
            IGNORED - vectorized computation is single-threaded but much faster
        
        Returns:
        --------
        np.ndarray : Distance matrix
        """
        n = len(universes)
        
        if align_structures:
            print("Aligning structures to reference...")
            reference = universes[reference_idx]
            for u in universes:
                try:
                    align.AlignTraj(u, reference, select=selection, in_memory=True).run()
                except Exception as e:
                    print(f"Warning: Could not align structure: {e}")
        else:
            print("Skipping alignment - computing RMSD on original coordinates")
        
        print(f"Computing pairwise RMSD matrix (VECTORIZED - {n} structures)...")
        
        # Extract all coordinates into a single 3D array: (n_structures, n_atoms, 3)
        print(f"  Extracting coordinates for {n} structures...")
        coords_list = []
        for i, u in enumerate(universes):
            try:
                coords = u.select_atoms(selection).positions
                coords_list.append(coords)
            except Exception as e:
                print(f"Warning: Could not extract coords from structure {i}: {e}")
                coords_list.append(None)
        
        # Check all have same number of atoms
        n_atoms_list = [c.shape[0] for c in coords_list if c is not None]
        if len(set(n_atoms_list)) > 1:
            print(f"Warning: Structures have different numbers of atoms: {set(n_atoms_list)}")
            print("Falling back to pairwise computation...")
            # Fallback to pairwise
            distance_matrix = np.zeros((n, n))
            for i in range(n):
                for j in range(i+1, n):
                    if coords_list[i] is not None and coords_list[j] is not None:
                        try:
                            # Handle different sizes by using minimum
                            min_atoms = min(coords_list[i].shape[0], coords_list[j].shape[0])
                            diff = coords_list[i][:min_atoms] - coords_list[j][:min_atoms]
                            rmsd_value = np.sqrt((diff ** 2).sum() / min_atoms)
                            distance_matrix[i, j] = rmsd_value
                            distance_matrix[j, i] = rmsd_value
                        except:
                            distance_matrix[i, j] = np.nan
                            distance_matrix[j, i] = np.nan
                    else:
                        distance_matrix[i, j] = np.nan
                        distance_matrix[j, i] = np.nan
            return distance_matrix
        
        # All structures have same number of atoms - proceed with vectorization
        n_atoms = n_atoms_list[0]
        print(f"  Structures have {n_atoms} atoms each")
        
        # Stack into 3D array: shape (n_structures, n_atoms, 3)
        coords_array = np.stack(coords_list, axis=0)
        
        print(f"  Computing RMSD matrix using vectorized operations...")
        print(f"  Array shape: {coords_array.shape}")
        
        # MEMORY-EFFICIENT VECTORIZED RMSD COMPUTATION
        # For large n, computing full (n, n, n_atoms, 3) array can use too much memory
        # Instead, compute in chunks
        chunk_size = 500  # Process 500 structures at a time
        rmsd_matrix = np.zeros((n, n))
        
        print(f"  Computing in chunks of {chunk_size} to manage memory...")
        for i_start in range(0, n, chunk_size):
            i_end = min(i_start + chunk_size, n)
            
            # Compute RMSD for chunk [i_start:i_end] against all structures
            coords_i = coords_array[i_start:i_end, np.newaxis, :, :]  # (chunk, 1, n_atoms, 3)
            coords_j = coords_array[np.newaxis, :, :, :]  # (1, n, n_atoms, 3)
            
            # Compute differences for this chunk
            diff = coords_i - coords_j  # (chunk, n, n_atoms, 3)
            squared_diff = diff ** 2
            sum_squared = squared_diff.sum(axis=(2, 3))  # (chunk, n)
            rmsd_chunk = np.sqrt(sum_squared / n_atoms)
            
            # Store in matrix
            rmsd_matrix[i_start:i_end, :] = rmsd_chunk
            
            if (i_start // chunk_size + 1) % 5 == 0:
                print(f"    Progress: {i_end}/{n} structures ({100*i_end/n:.1f}%)")
        
        # Ensure diagonal is exactly 0
        np.fill_diagonal(rmsd_matrix, 0.0)
        
        print(f"  ✓ RMSD matrix computed!")
        
        if normalize:
            max_dist = np.nanmax(rmsd_matrix)
            if max_dist > 0:
                rmsd_matrix = rmsd_matrix / max_dist
        
        return rmsd_matrix
    
    @staticmethod
    def compute_jaccard_matrix(pdb_files, distance_cutoff=4.5, 
                              protein_selection='protein', 
                              nucleic_selection='nucleic',
                              n_jobs=1, motif_residues=None, reference_sequence=None):
        """
        Compute Jaccard distance matrix based on protein-DNA interface contacts
        
        Key insight: Compares protein residues in contact, NOT residue pairs.
        This makes binding modes comparable even with different DNA sequences.
        
        RESIDUE OFFSET CORRECTION (for chopped AlphaFold models):
        - Parses AlphaFold filename to extract receptor start position
        - Optionally validates against reference sequence via alignment
        - Applies offset to all protein residue numbers before Jaccard comparison
        - This ensures "residue 5" means the same position across all structures
        
        PERFORMANCE OPTIMIZED:
        - Step 0: Motif screening (optional) - filters off-target poses before full contact search
        - Step 1: Parallelized contact detection (I/O intensive, benefits from parallelization)
        - Step 2: Vectorized distance calculation (CPU intensive, scipy pdist is fastest)
        
        Parameters:
        -----------
        pdb_files : list
            List of PDB file paths
        distance_cutoff : float
            Contact distance cutoff in Angstroms
        protein_selection : str
            MDTraj selection for protein
        nucleic_selection : str
            MDTraj selection for nucleic acid
        n_jobs : int
            Number of parallel jobs for contact detection (-1 for all cores, 1 for serial)
        motif_residues : dict, optional
            Dict mapping chain_id -> list of residue numbers for motif screening
            Example: {'A': [10, 11, 12], 'B': [30, 31, 32]}
            If provided, poses with no motif contacts are filtered out (marked as off-target)
        reference_sequence : str, optional
            Reference full-length amino acid sequence (single-letter codes)
            Used to validate/correct residue numbering for chopped AlphaFold models
        
        Returns:
        --------
        tuple : (distance_matrix, contact_sets, contact_residue_pairs, pdb_files_filtered, filtered_indices)
        """
        n = len(pdb_files)
        print(f"Computing Jaccard contact distance matrix (cutoff={distance_cutoff} Å)...")
        
        # STEP 0: Compute residue offsets for chopped AlphaFold models
        from core.io_utils import validate_and_correct_residue_offset
        residue_offsets = []
        offset_cache = {}
        
        if reference_sequence:
            print(f"  📏 RESIDUE CORRECTION: Validating residue numbering against reference sequence")
            print(f"  → Reference length: {len(reference_sequence)} residues")
        
        for pdb_file in pdb_files:
            offset_info = validate_and_correct_residue_offset(
                pdb_file, 
                reference_sequence=reference_sequence,
                cache=offset_cache
            )
            residue_offsets.append(offset_info['offset'])
            
            # Log validation results for debugging
            if reference_sequence and offset_info.get('validated'):
                method = offset_info.get('method', 'unknown')
                identity = offset_info.get('identity', 0)
                if method == 'alignment' and identity < 95:
                    print(f"  ⚠️  Low identity ({identity:.1f}%) for {os.path.basename(pdb_file)}")
        
        # Summary of offset corrections
        unique_offsets = set(residue_offsets)
        if len(unique_offsets) > 1:
            print(f"  ✓ Applied {len(unique_offsets)} different offsets across structures")
        elif list(unique_offsets)[0] != 0:
            print(f"  ✓ Applied uniform offset of +{list(unique_offsets)[0]} residues")
        else:
            print(f"  → No offset correction needed (all structures start at residue 1)")
        
        if motif_residues:
            print(f"  🎯 MOTIF MODE: Motif screening enabled")
            print(f"  → Off-target poses (no motif contact) will be filtered")
        else:
            print(f"  🌍 GLOBAL MODE: Processing all protein-nucleic contacts")
        print(f"Step 1/2: Computing interface contacts for all structures...")
        if n_jobs != 1:
            print(f"  Using {n_jobs} parallel jobs for contact detection")
        
        # STEP 1: Parallel contact detection (I/O bound - benefits from parallelization)
        # With optional motif screening to filter off-target poses
        # Pass residue offsets for chopped AlphaFold model correction
        results = Parallel(n_jobs=n_jobs, verbose=5 if n_jobs != 1 else 0)(
            delayed(DistanceMetrics._process_structure_contacts)(
                pdb_file, distance_cutoff, protein_selection, nucleic_selection, motif_residues, offset
            )
            for pdb_file, offset in zip(pdb_files, residue_offsets)
        )
        
        # Unpack results
        contact_sets = [r[0] for r in results]
        contact_residue_pairs = [r[1] for r in results]
        
        # Compute motif match scores (% of motif residues making contact)
        motif_match_scores = []
        if motif_residues is not None and len(motif_residues) > 0:
            # Count total motif residues
            total_motif_residues = sum(len(res_nums) for res_nums in motif_residues.values())
            
            for contact_set in contact_sets:
                # Count how many motif residues are in the contact set
                matched_motif_count = 0
                for chain_id, res_nums in motif_residues.items():
                    for res_num in res_nums:
                        # Check if any contact in contact_set matches this motif residue
                        # Contact format: "CHAIN:RESNAME123"
                        for contact_res in contact_set:
                            if contact_res.startswith(f"{chain_id}:"):
                                # Extract residue number from contact string
                                try:
                                    contact_res_num = int(''.join(filter(str.isdigit, contact_res.split(':')[1])))
                                    if contact_res_num == res_num:
                                        matched_motif_count += 1
                                        break
                                except (ValueError, IndexError):
                                    continue
                
                # Compute match score as percentage
                match_score = (matched_motif_count / total_motif_residues * 100) if total_motif_residues > 0 else 0.0
                motif_match_scores.append(match_score)
        else:
            # No motif provided - all scores are None
            motif_match_scores = [None] * len(contact_sets)
        
        # CRITICAL FIX: Remove structures with no contacts (off-target for motif search)
        # These inflate metrics and create artificial clusters at distance=0
        valid_indices = [i for i, cs in enumerate(contact_sets) if len(cs) > 0]
        filtered_indices = [i for i, cs in enumerate(contact_sets) if len(cs) == 0]
        
        filtered_count = len(filtered_indices)
        retained_count = len(valid_indices)
        
        if filtered_count > 0:
            print(f"  ✓ Filtering: {retained_count}/{n} poses retained ({filtered_count} filtered out)")
            print(f"  → Filtered poses (no contacts): {[os.path.basename(pdb_files[i]) for i in filtered_indices[:5]]}" + 
                  (f" ... and {filtered_count - 5} more" if filtered_count > 5 else ""))
            
            # Filter out structures with no contacts
            contact_sets = [contact_sets[i] for i in valid_indices]
            contact_residue_pairs = [contact_residue_pairs[i] for i in valid_indices]
            motif_match_scores = [motif_match_scores[i] for i in valid_indices]
            pdb_files_filtered = [pdb_files[i] for i in valid_indices]
            n = len(valid_indices)
            
            if n == 0:
                raise ValueError("All poses were filtered out (no contacts found). Check motif residues or SASA threshold.")
        else:
            print(f"  ✓ All {n} poses have contacts (none filtered)")
            # No filtering needed - use original lists
            pdb_files_filtered = pdb_files
        
        # Print motif match score summary if available
        if motif_residues is not None and len(motif_match_scores) > 0 and any(s is not None for s in motif_match_scores):
            valid_scores = [s for s in motif_match_scores if s is not None]
            if len(valid_scores) > 0:
                print(f"  ✓ Motif match scores: mean={np.mean(valid_scores):.1f}%, min={np.min(valid_scores):.1f}%, max={np.max(valid_scores):.1f}%")
        
        # STEP 2: Vectorized Jaccard distance computation (CPU bound - scipy is optimal)
        print(f"Step 2/2: Computing Jaccard distances (vectorized with scipy.pdist)...")
        
        contact_counts = [len(cs) for cs in contact_sets]
        print(f"  Contact statistics: min={min(contact_counts) if contact_counts else 0}, "
              f"max={max(contact_counts) if contact_counts else 0}, "
              f"mean={np.mean(contact_counts) if contact_counts else 0:.1f}")
        
        # VECTORIZED JACCARD COMPUTATION using scipy.spatial.distance.pdist
        # Convert contact sets to binary matrix representation
        # Get all unique residues across all structures
        all_residues = sorted(set().union(*contact_sets))
        residue_to_idx = {res: idx for idx, res in enumerate(all_residues)}
        
        # Create binary contact matrix: rows = structures, cols = residues
        # binary_matrix[i, j] = 1 if residue j is in contact in structure i
        binary_matrix = np.zeros((n, len(all_residues)), dtype=bool)
        for i, contact_set in enumerate(contact_sets):
            for residue in contact_set:
                binary_matrix[i, residue_to_idx[residue]] = True
        
        # Compute Jaccard distance using scipy's pdist with 'jaccard' metric
        # This is MUCH faster than nested loops (vectorized C implementation)
        # pdist returns condensed distance matrix (upper triangle only)
        jaccard_distances = pdist(binary_matrix, metric='jaccard')
        
        # Convert to square form (full symmetric matrix)
        jaccard_matrix = squareform(jaccard_distances)
        
        # No need to handle NaN values anymore - all structures have contacts (filtered above)
        # Diagonal should be 0 (identical to self)
        np.fill_diagonal(jaccard_matrix, 0.0)
        
        # Summary stats
        non_diag = jaccard_matrix[np.triu_indices(n, k=1)]
        if len(non_diag) > 0:
            print(f"Jaccard distance range: {non_diag.min():.3f} - {non_diag.max():.3f}")
            print(f"  Mean distance: {non_diag.mean():.3f} (closer to 0 = more similar)")
        
        # Store motif scores for later retrieval if needed
        DistanceMetrics._last_motif_match_scores = motif_match_scores
        
        return jaccard_matrix, contact_sets, contact_residue_pairs, pdb_files_filtered, filtered_indices
