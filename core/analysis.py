"""
Analysis tools for interface characterization
"""

import pandas as pd
import os
from collections import Counter


class InterfaceAnalyzer:
    """Analyze protein-DNA interfaces from clustering results"""
    
    @staticmethod
    def get_interface_stats(pdb_files, contact_residue_pairs, labels=None, representative_indices=None):
        """
        Get interface statistics for all structures
        
        Parameters:
        -----------
        pdb_files : list
            List of PDB file paths
        contact_residue_pairs : list
            List of contact details for each structure
        labels : array-like, optional
            Cluster labels
        representative_indices : dict, optional
            Dict mapping cluster_id -> structure index
        
        Returns:
        --------
        pd.DataFrame with interface statistics
        """
        rows = []
        for idx, pdb_file in enumerate(pdb_files):
            structure_name = os.path.basename(pdb_file)
            cluster_id = labels[idx] if labels is not None else -1
            
            # Get contacts for this structure
            contacts = contact_residue_pairs[idx]
            
            # Extract protein and nucleic residues
            protein_res = set()
            nucleic_res = set()
            for contact in contacts:
                protein_res.add(contact['protein_residue'])
                nucleic_res.add(contact['nucleic_residue'])
            
            # Check if representative
            is_rep = False
            if representative_indices:
                # representative_indices is a list of indices
                is_rep = idx in representative_indices
            
            rows.append({
                'structure': structure_name,
                'cluster': cluster_id,
                'n_protein_residues': len(protein_res),
                'n_nucleic_residues': len(nucleic_res),
                'protein_residues': ','.join(sorted(protein_res)),
                'nucleic_residues': ','.join(sorted(nucleic_res)),
                'representative': is_rep
            })
        
        return pd.DataFrame(rows)
    
    @staticmethod
    def get_binding_hotspots(contact_residue_pairs, labels=None, cluster_id=None):
        """
        Compute binding hotspots - protein residues most frequently in contact
        
        Parameters:
        -----------
        contact_residue_pairs : list
            List of contact details for each structure
        labels : array-like, optional
            Cluster labels
        cluster_id : int, optional
            If provided, compute hotspots only for this cluster
        
        Returns:
        --------
        pd.DataFrame with columns: residue, frequency, percentage
        """
        # Filter structures by cluster if requested
        if cluster_id is not None:
            if labels is None:
                raise ValueError("No clustering labels available")
            indices = [i for i, label in enumerate(labels) if label == cluster_id]
        else:
            indices = range(len(contact_residue_pairs))
        
        # Count protein residue occurrences
        residue_counts = Counter()
        for idx in indices:
            contacts = contact_residue_pairs[idx]
            protein_res_in_struct = set(c['protein_residue'] for c in contacts)
            for res in protein_res_in_struct:
                residue_counts[res] += 1
        
        # Create DataFrame
        total_structures = len(indices)
        hotspot_data = []
        for residue, count in residue_counts.most_common():
            hotspot_data.append({
                'residue': residue,
                'frequency': count,
                'percentage': 100 * count / total_structures if total_structures > 0 else 0
            })
        
        # Return empty DataFrame with correct columns if no data
        if not hotspot_data:
            return pd.DataFrame(columns=['residue', 'frequency', 'percentage'])
        
        return pd.DataFrame(hotspot_data)
    
    @staticmethod
    def get_cluster_signature(contact_residue_pairs, labels, cluster_id, threshold=50):
        """
        Get consensus interface residues for a cluster
        
        Parameters:
        -----------
        contact_residue_pairs : list
            List of contact details
        labels : array-like
            Cluster labels
        cluster_id : int
            Cluster ID
        threshold : float
            Percentage threshold for consensus (default: 50%)
        
        Returns:
        --------
        set : Protein residues present in >threshold% of cluster members
        """
        hotspots_df = InterfaceAnalyzer.get_binding_hotspots(
            contact_residue_pairs, labels, cluster_id
        )
        
        consensus_residues = set(
            hotspots_df[hotspots_df['percentage'] > threshold]['residue'].tolist()
        )
        
        return consensus_residues
    
    @staticmethod
    def get_cluster_summary(contact_residue_pairs, labels, threshold=50):
        """
        Get summary of all clusters with their consensus residues
        
        Parameters:
        -----------
        contact_residue_pairs : list
            List of contact details
        labels : array-like
            Cluster labels
        threshold : float
            Percentage threshold for consensus (default: 50%)
        
        Returns:
        --------
        pd.DataFrame with columns: cluster, n_structures, consensus_residues, n_consensus
        """
        unique_clusters = sorted(set(labels))
        rows = []
        
        for cluster_id in unique_clusters:
            # Count structures
            n_structures = sum(1 for label in labels if label == cluster_id)
            
            # Get consensus residues
            consensus = InterfaceAnalyzer.get_cluster_signature(
                contact_residue_pairs, labels, cluster_id, threshold
            )
            
            # Get all hotspots for this cluster
            hotspots_df = InterfaceAnalyzer.get_binding_hotspots(
                contact_residue_pairs, labels, cluster_id
            )
            
            # Get top hotspots (>threshold%)
            top_hotspots = hotspots_df[hotspots_df['percentage'] > threshold].sort_values(
                'percentage', ascending=False
            )
            
            rows.append({
                'cluster': cluster_id,
                'n_structures': n_structures,
                'n_consensus': len(consensus),
                'consensus_residues': ','.join(sorted(consensus)),
                'top_residues': ','.join(top_hotspots.head(10)['residue'].tolist()),
                'binding_mode': 'Noise' if cluster_id == -1 else f'Mode {cluster_id}'
            })
        
        return pd.DataFrame(rows)

    @staticmethod
    def calculate_delta_sasa(pdb_file, protein_chains, nucleic_chains, motif_residues=None, split_by_chain=False,
                            precomputed_protein_sasa=None, precomputed_nucleic_sasa=None,
                            precomputed_structures=None):
        """
        Calculate ΔSASA (buried surface area) for protein-nucleic acid complex.
        
        PERFORMANCE: For multiple poses with the same protein/nucleic structures,
        pass precomputed_protein_sasa and precomputed_nucleic_sasa to avoid redundant calculations.
        
        Three modes:
        1. Default (no motif): Total protein-DNA interface ΔSASA
        2. Motif mode: ΔSASA for specific motif residues only
        3. Split by chain: Per-chain ΔSASA (requires motif to be meaningful)
        
        ΔSASA = (SASA_isolated + SASA_partner) - SASA_complex
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        protein_chains : list
            List of protein chain IDs
        nucleic_chains : list
            List of nucleic acid chain IDs
        motif_residues : dict, optional
            Dict mapping chain_id -> list of residue numbers for motif-specific ΔSASA
            Example: {'A': [10, 11, 12], 'B': [30, 31, 32]}
        split_by_chain : bool, optional
            If True, calculate ΔSASA for each motif chain separately.
            Note: Split mode calculates per-chain contributions to the motif interface.
        precomputed_protein_sasa : dict, optional
            Pre-computed protein SASA {'sasa': float, 'structure': freesasa.Structure, 'result': freesasa.Result}
        precomputed_nucleic_sasa : dict, optional
            Pre-computed nucleic SASA {'sasa': float, 'structure': freesasa.Structure, 'result': freesasa.Result}
        precomputed_structures : dict, optional
            Pre-computed protein/nucleic PDB paths {'protein_pdb': str, 'nucleic_pdb': str}
        
        Returns:
        --------
        dict : {'total_delta_sasa': float, 'motif_delta_sasa': float or None,
                'sasa_complex': float, 'sasa_protein': float, 'sasa_nucleic': float,
                'chain_sasa': dict (if split_by_chain=True and motif provided)}
        """
        try:
            import freesasa
            import MDAnalysis as mda
            import tempfile
            import os
            import warnings
        except ImportError as e:
            raise ImportError(f"Required library not available: {e}")
        
        # Silence FreeSASA warnings about unknown atoms
        freesasa.setVerbosity(freesasa.silent)
        
        # Load structure with MDAnalysis
        u = mda.Universe(pdb_file)
        
        # Calculate SASA for complex with nucleic acid parameters
        # Use protor classifier which includes nucleic acids
        structure_complex = freesasa.Structure(pdb_file)
        params = freesasa.Parameters()
        params.setProbeRadius(1.4)  # Standard water probe
        result_complex = freesasa.calc(structure_complex, params)
        sasa_complex = result_complex.totalArea()
        
        # Create temporary PDB files for separated structures (keep for debugging)
        basename = os.path.basename(pdb_file).replace('.pdb', '')
        debug_dir = os.path.join(os.path.dirname(pdb_file), 'sasa_debug')
        os.makedirs(debug_dir, exist_ok=True)
        
        protein_pdb = os.path.join(debug_dir, f"{basename}_protein.pdb")
        nucleic_pdb = os.path.join(debug_dir, f"{basename}_nucleic.pdb")
        
        try:
            # Use precomputed values if provided (PERFORMANCE OPTIMIZATION)
            if precomputed_protein_sasa is not None and precomputed_nucleic_sasa is not None:
                # Use cached values
                sasa_protein = precomputed_protein_sasa['sasa']
                sasa_nucleic = precomputed_nucleic_sasa['sasa']
                structure_protein = precomputed_protein_sasa['structure']
                result_protein = precomputed_protein_sasa['result']
                structure_nucleic = precomputed_nucleic_sasa['structure']
                result_nucleic = precomputed_nucleic_sasa['result']
                
                if precomputed_structures:
                    protein_pdb = precomputed_structures['protein_pdb']
                    nucleic_pdb = precomputed_structures['nucleic_pdb']
            else:
                # Compute from scratch (original behavior)
                # Select protein and nucleic acid by type (not by chain)
                protein_atoms = u.select_atoms("protein")
                nucleic_atoms = u.select_atoms("nucleic")
                
                # Write separated structures (will be saved to debug folder)
                protein_atoms.write(protein_pdb)
                nucleic_atoms.write(nucleic_pdb)
                
                # Calculate SASA for isolated protein
                structure_protein = freesasa.Structure(protein_pdb)
                result_protein = freesasa.calc(structure_protein, params)
                sasa_protein = result_protein.totalArea()
                
                # Calculate SASA for isolated nucleic acid with nucleic acid parameters
                structure_nucleic = freesasa.Structure(nucleic_pdb)
                result_nucleic = freesasa.calc(structure_nucleic, params)
                sasa_nucleic = result_nucleic.totalArea()
            
            # Calculate ΔSASA: buried surface area
            delta_sasa_total = (sasa_protein + sasa_nucleic) - sasa_complex
            
            # Calculate per-chain motif ΔSASA if requested
            # This calculates ΔSASA for motif residues on each chain separately
            chain_sasa = {}
            if split_by_chain and motif_residues:
                # print(f"\n  Split-by-chain mode: Calculating per-chain motif ΔSASA")
                
                # Only process chains that have motif residues defined
                for chain_id in motif_residues.keys():
                    residue_nums = motif_residues[chain_id]
                    
                    # Create selection for motif residues on this chain only
                    # MDAnalysis resid selection syntax: "resid 10 11 12" or "resid 10:12"
                    resid_str = " ".join(str(r) for r in residue_nums)
                    chain_motif_sel = f"segid {chain_id} and resid {resid_str}"
                    chain_motif_atoms = u.select_atoms(chain_motif_sel)
                    
                    # Fallback to chainID if segid didn't work
                    if chain_motif_atoms.n_atoms == 0:
                        chain_motif_sel = f"chainID {chain_id} and resid {resid_str}"
                        chain_motif_atoms = u.select_atoms(chain_motif_sel)
                    
                    # print(f"  Chain {chain_id} motif: residues {residue_nums}")
                    # print(f"    Selection: '{chain_motif_sel}'")
                    # print(f"    Motif atoms selected: {chain_motif_atoms.n_atoms}")
                    
                    # DEBUG: Show what residues were actually selected
                    if chain_motif_atoms.n_atoms > 0:
                        selected_resnums = set(chain_motif_atoms.resnums)
                        # print(f"    Actually selected residues: {sorted(selected_resnums)}")
                    
                    if chain_motif_atoms.n_atoms > 0:
                        # Create debug file for this chain
                        chain_pdb = os.path.join(debug_dir, f"{basename}_chain_{chain_id}_motif.pdb")
                        chain_motif_atoms.write(chain_pdb)
                        # print(f"    Saved to: {chain_pdb}")
                        
                        # Calculate motif SASA for this chain using selectArea
                        # Build FreeSASA selection command for this chain's motif residues
                        # Format: "s1, chain A and resi 11+31+32+33"
                        # IMPORTANT: Filter out nucleic acid residues since DNA/RNA may share chain labels
                        resi_str = "+".join(str(r) for r in residue_nums)
                        selection_cmd = f"s1, chain {chain_id} and resi {resi_str} and not resn DA+DT+DC+DG+DA5+DT3+A+U+C+G+DA3+DC3+DG3+DC5+DG5"
                        
                        # Use selectArea on isolated protein to get motif SASA without DNA
                        selection_result_protein = freesasa.selectArea([selection_cmd], structure_protein, result_protein)
                        sasa_chain_in_protein = selection_result_protein.get('s1', 0.0)
                        
                        # Use selectArea on complex to get motif SASA with DNA
                        selection_result_complex = freesasa.selectArea([selection_cmd], structure_complex, result_complex)
                        sasa_chain_in_complex = selection_result_complex.get('s1', 0.0)
                        
                        # print(f"    DEBUG: Selection command: {selection_cmd}")
                        # print(f"    DEBUG: Atoms in motif selection: {chain_motif_atoms.n_atoms}")
                        
                        # Chain motif ΔSASA
                        chain_sasa[chain_id] = sasa_chain_in_protein - sasa_chain_in_complex
                        
                        # print(f"    In protein: {sasa_chain_in_protein:.2f} Å²")
                        # print(f"    In complex: {sasa_chain_in_complex:.2f} Å²")
                        # print(f"    Chain {chain_id} motif ΔSASA: {chain_sasa[chain_id]:.2f} Å²")
                    else:
                        # print(f"    WARNING: No atoms found for chain {chain_id} motif, ΔSASA = 0")
                        chain_sasa[chain_id] = 0.0
            elif split_by_chain and not motif_residues:
                pass  # print(f"\n  WARNING: Split-by-chain mode requires motif residues to be defined. Skipping.")
            
            # DEBUG: Print component values (commented for performance with large datasets)
            # basename_short = os.path.basename(pdb_file)
            # print(f"\nDEBUG ΔSASA for {basename_short}:")
            # print(f"  Complex SASA: {sasa_complex:.2f} Å²")
            # print(f"  Isolated Protein SASA: {sasa_protein:.2f} Å² (should be same for all models)")
            # print(f"  Isolated Nucleic SASA: {sasa_nucleic:.2f} Å² (should be same for all models)")
            # print(f"  Sum of isolated: {sasa_protein + sasa_nucleic:.2f} Å²")
            # print(f"  ΔSASA (buried): {delta_sasa_total:.2f} Å²")
            # if split_by_chain and chain_sasa:
            #     print(f"  Per-chain motif ΔSASA:")
            #     for chain_id, sasa_val in chain_sasa.items():
            #         print(f"    Chain {chain_id}: {sasa_val:.2f} Å²")
            # print()
            
        finally:
            # Keep debug files for inspection - don't delete them
            # User can manually inspect sasa_debug folder to verify splits
            pass
        
        # Motif-specific ΔSASA if requested - ONLY if motif_residues is provided
        delta_sasa_motif = None
        sasa_motif_in_protein = None
        
        if motif_residues is not None and len(motif_residues) > 0:
            # print(f"\n  Motif residues provided: {motif_residues}")
            # print(f"  Starting motif ΔSASA calculation...")
            try:
                # Create selection for motif residues
                motif_selections = []
                for chain_id, residue_nums in motif_residues.items():
                    for resnum in residue_nums:
                        motif_selections.append(f"(segid {chain_id} and resid {resnum})")
                
                if motif_selections:
                    motif_sel = ' or '.join(motif_selections)
                    
                    # Try to select motif atoms
                    motif_atoms = u.select_atoms(motif_sel)
                    
                    # Fallback to chainID if segid didn't work
                    if motif_atoms.n_atoms == 0:
                        motif_selections = []
                        for chain_id, residue_nums in motif_residues.items():
                            for resnum in residue_nums:
                                motif_selections.append(f"(chainID {chain_id} and resid {resnum})")
                        motif_sel = ' or '.join(motif_selections)
                        motif_atoms = u.select_atoms(motif_sel)
                    
                    # print(f"\n  Motif selection: {motif_residues}")
                    # print(f"  Motif atoms selected: {motif_atoms.n_atoms}")
                    
                    if motif_atoms.n_atoms > 0:
                        # Write motif-only PDB for debugging
                        motif_pdb = os.path.join(debug_dir, f"{basename}_motif.pdb")
                        motif_atoms.write(motif_pdb)
                        # print(f"  Saved motif to: {motif_pdb}")
                        
                        # Use conventional ΔSASA calculation: motif + nucleic - complex
                        # Build FreeSASA selection commands for all chains with motif residues
                        # IMPORTANT: Filter out nucleic acid residues since DNA/RNA may share chain labels with protein
                        selection_cmds = []
                        for chain_id, residue_nums in motif_residues.items():
                            resi_str = "+".join(str(r) for r in residue_nums)
                            # Exclude common nucleic acid residue names
                            selection_cmds.append(f"s_{chain_id}, chain {chain_id} and resi {resi_str} and not resn DA+DT+DC+DG+DA5+DT3+A+U+C+G+DA3+DC3+DG3+DC5+DG5")
                        
                        # Get motif SASA from complex (this is in the context of protein+DNA together)
                        selection_result_complex = freesasa.selectArea(selection_cmds, structure_complex, result_complex)
                        sasa_motif_in_complex = sum(selection_result_complex.values())
                        
                        # Get motif SASA from isolated protein (this is in the context of protein alone)
                        selection_result_protein = freesasa.selectArea(selection_cmds, structure_protein, result_protein)
                        sasa_motif_in_protein = sum(selection_result_protein.values())
                        
                        # Motif ΔSASA = (motif_in_protein + nucleic) - motif_in_complex
                        # But we need to account for the nucleic SASA contribution
                        # Simpler: just use the difference - buried surface on motif
                        delta_sasa_motif = sasa_motif_in_protein - sasa_motif_in_complex
                        
                        # print(f"  DEBUG: Selection commands: {selection_cmds}")
                        # print(f"  DEBUG: Motif in isolated protein: {selection_result_protein} = {sasa_motif_in_protein:.2f} Ų")
                        # print(f"  DEBUG: Motif in complex: {selection_result_complex} = {sasa_motif_in_complex:.2f} Ų")
                        # print(f"  DEBUG: structure_protein nAtoms: {structure_protein.nAtoms()}")
                        # print(f"  DEBUG: structure_complex nAtoms: {structure_complex.nAtoms()}")
                        
                        # print(f"  Motif in protein: {sasa_motif_in_protein:.2f} Å²")
                        # print(f"  Motif in complex: {sasa_motif_in_complex:.2f} Å²")
                        # print(f"  Motif ΔSASA: {delta_sasa_motif:.2f} Å²")
                    else:
                        pass  # print(f"  WARNING: No atoms found for motif selection!")
                        delta_sasa_motif = 0.0
                        sasa_motif_in_protein = 0.0
                            
            except Exception as e:
                # If motif calculation fails, just skip it
                # print(f"  ERROR in motif ΔSASA calculation: {e}")
                delta_sasa_motif = None
                sasa_motif_in_protein = None
        # else:
            # print(f"\n  No motif residues provided - skipping motif ΔSASA calculation")
        
        result = {
            'total_delta_sasa': delta_sasa_total,
            'motif_delta_sasa': delta_sasa_motif,
            'motif_sasa_protein': sasa_motif_in_protein,  # Available surface area of motif
            'sasa_complex': sasa_complex,
            'sasa_protein': sasa_protein,
            'sasa_nucleic': sasa_nucleic
        }
        
        # Add per-chain ΔSASA if calculated
        if split_by_chain and chain_sasa:
            result['chain_sasa'] = chain_sasa
        
        return result
    
    @staticmethod
    def parse_motif_residues(motif_string):
        """
        Parse motif residues string into dict
        
        Format: "Chain:Start-End, Chain:Start-End, ..."
        Example: "A:10-20, B:30-35"
        
        Returns:
        --------
        dict : {chain_id: [residue_numbers]}
        """
        if not motif_string or not motif_string.strip():
            return None
        
        motif_dict = {}
        parts = motif_string.split(',')
        
        for part in parts:
            part = part.strip()
            if ':' not in part:
                continue
            
            chain_part, range_part = part.split(':', 1)
            chain_id = chain_part.strip()
            
            if '-' in range_part:
                start, end = range_part.split('-')
                residues = list(range(int(start.strip()), int(end.strip()) + 1))
            else:
                residues = [int(range_part.strip())]
            
            if chain_id in motif_dict:
                motif_dict[chain_id].extend(residues)
            else:
                motif_dict[chain_id] = residues
        
        return motif_dict if motif_dict else None

    @staticmethod
    def calculate_per_residue_sasa(pdb_file, protein_chains, nucleic_chains, output_pdb=None):
        """
        Calculate per-ATOM ΔSASA and optionally write to PDB with values in B-factor column.
        
        This calculates atomic-level SASA for proper visualization - each atom gets its own ΔSASA value.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        protein_chains : list
            List of protein chain IDs
        nucleic_chains : list
            List of nucleic acid chain IDs
        output_pdb : str, optional
            Path to output PDB file with per-atom ΔSASA in B-factor column
            
        Returns:
        --------
        dict : {chain_id: {resnum: delta_sasa}} (averaged per residue for compatibility)
        """
        import MDAnalysis as mda
        import freesasa
        import os
        
        # Load structure
        u = mda.Universe(pdb_file)
        
        # Setup FreeSASA parameters
        params = freesasa.Parameters()
        params.setProbeRadius(1.4)
        freesasa.setVerbosity(freesasa.silent)
        
        # Calculate SASA for full complex
        structure_complex = freesasa.Structure(pdb_file)
        result_complex = freesasa.calc(structure_complex, params)
        
        # Calculate SASA for isolated protein
        protein_atoms = u.select_atoms("protein")
        temp_protein_pdb = pdb_file.replace('.pdb', '_temp_protein.pdb')
        protein_atoms.write(temp_protein_pdb)
        
        structure_protein = freesasa.Structure(temp_protein_pdb)
        result_protein = freesasa.calc(structure_protein, params)
        
        # Get per-ATOM SASA using FreeSASA's atomArea() function
        # This gives us individual SASA values for each atom
        per_atom_sasa = {}  # {(chain, resnum, atom_index): delta_sasa}
        per_residue_sasa = {}  # {chain_id: {resnum: avg_delta_sasa}}
        
        # Build atom index mapping for protein structure
        # FreeSASA atoms are indexed sequentially
        atom_idx = 0
        for atom in protein_atoms:
            chain_id = atom.chainID
            resnum = atom.resnum
            
            # Get atomic SASA from isolated protein
            try:
                sasa_protein = result_protein.atomArea(atom_idx)
            except:
                sasa_protein = 0.0
            
            # Find corresponding atom in complex structure
            # We need to match by chain, residue, and atom name
            # FreeSASA indexes atoms in order they appear
            try:
                # For complex, we need to find the atom index
                # This is tricky - we'll use the full atom list
                complex_atom_idx = None
                complex_idx = 0
                for catom in u.atoms:
                    if (catom.chainID == chain_id and 
                        catom.resnum == resnum and 
                        catom.name == atom.name):
                        complex_atom_idx = complex_idx
                        break
                    complex_idx += 1
                
                if complex_atom_idx is not None:
                    sasa_complex = result_complex.atomArea(complex_atom_idx)
                else:
                    sasa_complex = 0.0
            except:
                sasa_complex = 0.0
            
            # Calculate per-atom ΔSASA
            delta_sasa = sasa_protein - sasa_complex
            per_atom_sasa[(chain_id, resnum, atom_idx)] = delta_sasa
            
            # Accumulate for per-residue average
            if chain_id not in per_residue_sasa:
                per_residue_sasa[chain_id] = {}
            if resnum not in per_residue_sasa[chain_id]:
                per_residue_sasa[chain_id][resnum] = []
            per_residue_sasa[chain_id][resnum].append(delta_sasa)
            
            atom_idx += 1
        
        # Average per-residue values for return compatibility
        for chain_id in per_residue_sasa:
            for resnum in per_residue_sasa[chain_id]:
                values = per_residue_sasa[chain_id][resnum]
                per_residue_sasa[chain_id][resnum] = sum(values) / len(values) if values else 0.0
        
        # Clean up temp file
        if os.path.exists(temp_protein_pdb):
            os.remove(temp_protein_pdb)
        
        # Write to PDB with per-ATOM ΔSASA in B-factor column if requested
        if output_pdb:
            # Build atom-level lookup for fast access
            atom_sasa_lookup = {}  # {(chain, resnum, atom_name): delta_sasa}
            protein_atom_idx = 0
            for atom in protein_atoms:
                chain_id = atom.chainID
                resnum = atom.resnum
                atom_name = atom.name
                
                if (chain_id, resnum, protein_atom_idx) in per_atom_sasa:
                    atom_sasa_lookup[(chain_id, resnum, atom_name)] = per_atom_sasa[(chain_id, resnum, protein_atom_idx)]
                
                protein_atom_idx += 1
            
            # Read original PDB and modify B-factors with per-atom values
            with open(pdb_file, 'r') as f:
                lines = f.readlines()
            
            with open(output_pdb, 'w') as f:
                for line in lines:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # Parse PDB line
                        chain = line[21:22].strip()
                        atom_name = line[12:16].strip()
                        try:
                            resnum = int(line[22:26].strip())
                        except:
                            resnum = None
                        
                        # Get per-ATOM ΔSASA
                        if (chain, resnum, atom_name) in atom_sasa_lookup:
                            delta_sasa = atom_sasa_lookup[(chain, resnum, atom_name)]
                            # Keep raw ΔSASA value for accurate coloring
                            bfactor = delta_sasa
                            
                            # Ensure line is long enough (pad if necessary)
                            if len(line) < 66:
                                line = line.rstrip('\n')
                                line = line + ' ' * (66 - len(line)) + '\n'
                            
                            # Replace B-factor: columns 61-66 (0-indexed: 60-66)
                            new_line = line[:60] + f"{bfactor:6.2f}" + line[66:]
                            f.write(new_line)
                        else:
                            f.write(line)
                    else:
                        f.write(line)
        
        return per_residue_sasa
    def create_sasa_viewer_pdb(pdb_file, motif_residues=None, contacts=None):
        """
        Create a PDB file with ΔSASA values in B-factor column for visualization
        
        This is a lightweight function that:
        1. Extracts protein and nucleic acid chains from contact data
        2. Calculates per-residue ΔSASA
        3. Writes a new PDB with ΔSASA in B-factor column
        4. Returns the path to the new PDB
        
        Note: motif_residues is ignored - we only visualize total ΔSASA
        
        Parameters:
        -----------
        pdb_file : str
            Path to input PDB file
        motif_residues : dict, optional
            Ignored (kept for compatibility)
        contacts : list, optional
            Contact residue pairs from clustering analysis
            Each contact is a dict with 'protein_residue' and 'nucleic_residue'
        
        Returns:
        --------
        str : Path to output PDB with ΔSASA in B-factors
        """
        import MDAnalysis as mda
        
        # Extract protein and nucleic chains from contact data
        protein_chains = set()
        nucleic_chains = set()
        
        if contacts:
            # Use contact data (most reliable)
            for contact in contacts:
                prot_res = contact['protein_residue']
                nuc_res = contact['nucleic_residue']
                
                if ':' in prot_res:
                    protein_chains.add(prot_res.split(':')[0])
                if ':' in nuc_res:
                    nucleic_chains.add(nuc_res.split(':')[0])
        else:
            # Fallback: Auto-detect from PDB (less reliable for mixed chains)
            print("   ⚠️ No contact data - attempting auto-detection")
            u = mda.Universe(pdb_file)
            
            # Nucleic acid residues
            nucleic_resnames = ['DA', 'DT', 'DG', 'DC', 'A', 'U', 'G', 'C', 
                               'DA5', 'DT3', 'DA3', 'DC3', 'DG3', 'DC5', 'DG5']
            
            # Try to detect chains by content
            # This won't work well if protein and nucleic share the same chain ID
            for chain in u.segments:
                chain_id = chain.segid
                chain_atoms = u.select_atoms(f"segid {chain_id}")
                
                nucleic_atoms = chain_atoms.select_atoms(f"resname {' '.join(nucleic_resnames)}")
                if len(nucleic_atoms) > 0:
                    nucleic_chains.add(chain_id)
                
                protein_atoms = chain_atoms.select_atoms(f"protein")
                if len(protein_atoms) > 0:
                    protein_chains.add(chain_id)
            
            # If no segid, try chainID
            if len(protein_chains) == 0 and len(nucleic_chains) == 0:
                for chain in set(u.atoms.chainIDs):
                    if not chain:
                        continue
                    chain_atoms = u.select_atoms(f"chainID {chain}")
                    
                    nucleic_atoms = chain_atoms.select_atoms(f"resname {' '.join(nucleic_resnames)}")
                    if len(nucleic_atoms) > 0:
                        nucleic_chains.add(chain)
                    
                    protein_atoms = chain_atoms.select_atoms(f"protein")
                    if len(protein_atoms) > 0:
                        protein_chains.add(chain)
        
        protein_chains = list(protein_chains)
        nucleic_chains = list(nucleic_chains)
        
        if not protein_chains or not nucleic_chains:
            print(f"   ✗ Could not detect chains. Protein: {protein_chains}, Nucleic: {nucleic_chains}")
            return None
        
        # Create output path in a temp directory to avoid polluting the PDB directory
        import tempfile
        temp_dir = tempfile.gettempdir()
        base_name = os.path.basename(pdb_file).replace('.pdb', '_sasa_viewer.pdb')
        output_pdb = os.path.join(temp_dir, base_name)
        
        print(f"\n🔬 Creating ΔSASA visualization PDB:")
        print(f"   Input: {os.path.basename(pdb_file)}")
        print(f"   Protein chains: {protein_chains}")
        print(f"   Nucleic chains: {nucleic_chains}")
        print(f"   Output: {output_pdb} (temporary)")
        
        # Calculate per-residue ΔSASA and write directly to PDB
        # This function does both calculation and PDB writing in one step
        per_residue_sasa = InterfaceAnalyzer.calculate_per_residue_sasa(
            pdb_file,
            protein_chains=protein_chains,
            nucleic_chains=nucleic_chains,
            output_pdb=output_pdb
        )
        
        if per_residue_sasa:
            # Show statistics about the ΔSASA values
            all_values = []
            for chain_dict in per_residue_sasa.values():
                all_values.extend(chain_dict.values())
            
            if all_values:
                min_val = min(all_values)
                max_val = max(all_values)
                mean_val = sum(all_values)/len(all_values)
                
                print(f"   ✓ ΔSASA statistics:")
                print(f"     - Residues with data: {len(all_values)}")
                print(f"     - Range: {min_val:.2f} - {max_val:.2f} Ų")
                print(f"     - Mean: {mean_val:.2f} Ų")
                
                # Write metadata to a companion file for py3Dmol to use
                meta_file = output_pdb.replace('.pdb', '_meta.txt')
                with open(meta_file, 'w') as f:
                    f.write(f"min={min_val:.2f}\n")
                    f.write(f"max={max_val:.2f}\n")
                    f.write(f"mean={mean_val:.2f}\n")
                print(f"     - Color range metadata: {os.path.basename(meta_file)}")
        
        if per_residue_sasa and output_pdb and os.path.exists(output_pdb):
            print(f"   ✓ PDB file created successfully\n")
            return output_pdb
        else:
            print(f"   ✗ Failed to create PDB\n")
            return None
