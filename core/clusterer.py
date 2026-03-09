"""
Main clustering class for PDB structures
"""

import os
import numpy as np
import MDAnalysis as mda
from sklearn.cluster import DBSCAN
from scipy.cluster.hierarchy import linkage, fcluster

from .metrics import DistanceMetrics
from .analysis import InterfaceAnalyzer


class PDBClusterer:
    """
    Main class for clustering PDB structures using various distance metrics
    """
    
    def __init__(self, pdb_files, selection='protein and name CA', align_selection=None):
        """
        Initialize the PDB clusterer
        
        Parameters:
        -----------
        pdb_files : list
            List of paths to PDB files
        selection : str
            MDAnalysis selection string for RMSD calculation
        align_selection : str
            MDAnalysis selection string for alignment
        """
        self.pdb_files = sorted(pdb_files)
        self.selection = selection
        self.align_selection = align_selection or selection
        self.n_structures = len(pdb_files)
        
        # State
        self.distance_matrix = None
        self.metric_name = "Unknown"
        self.metric_unit = ""
        self.labels = None
        self.universes = []
        self.representative_indices = {}
        
        # Jaccard-specific data
        self.contact_sets = None
        self.contact_residue_pairs = None
        
        # Residue offset correction
        self.reference_sequence = None  # Full-length reference sequence
        
        print(f"Found {self.n_structures} PDB files")
    
    def load_structures(self):
        """Load all PDB structures using MDAnalysis"""
        print("Loading structures...")
        self.universes = []
        
        for pdb_file in self.pdb_files:
            try:
                u = mda.Universe(pdb_file)
                self.universes.append(u)
            except Exception as e:
                print(f"Warning: Could not load {pdb_file}: {e}")
        
        print(f"Successfully loaded {len(self.universes)} structures")
    
    def compute_distance_matrix(self, reference_idx=0, align_structures=True, normalize=False, n_jobs=1):
        """
        Compute pairwise RMSD distance matrix
        
        Parameters:
        -----------
        reference_idx : int
            Index of reference structure for alignment
        align_structures : bool
            Whether to align structures before RMSD
        normalize : bool
            Whether to normalize to [0, 1]
        n_jobs : int
            Number of parallel jobs (-1 for all cores)
        """
        if len(self.universes) == 0:
            raise ValueError("No structures loaded. Call load_structures() first.")
        
        print(f"Computing pairwise RMSD distance matrix...")
        
        self.distance_matrix = DistanceMetrics.compute_rmsd_matrix(
            self.universes,
            selection=self.selection,
            align_structures=align_structures,
            reference_idx=reference_idx,
            normalize=normalize,
            n_jobs=n_jobs
        )
        
        self.metric_name = "RMSD"
        self.metric_unit = "Å"
        
        non_zero = self.distance_matrix[self.distance_matrix > 0]
        if len(non_zero) > 0:
            print(f"RMSD range: {non_zero.min():.2f} - {non_zero.max():.2f} {self.metric_unit}")
        print(f"Distance matrix shape: {self.distance_matrix.shape}")
    
    def compute_jaccard_contact_matrix(self, distance_cutoff=4.5,
                                       protein_selection='protein',
                                       nucleic_selection='nucleic',
                                       n_jobs=1, motif_residues=None,
                                       apply_offset=False):
        """
        Compute Jaccard distance matrix based on protein-DNA contacts
        
        Parameters:
        -----------
        distance_cutoff : float
            Contact distance cutoff in Angstroms
        protein_selection : str
            MDTraj selection for protein
        nucleic_selection : str
            MDTraj selection for nucleic acid
        n_jobs : int
            Number of parallel jobs (-1 for all cores)
        motif_residues : dict, optional
            Dict mapping chain_id -> list of residue numbers for motif screening
            Example: {'A': [10, 11, 12], 'B': [30, 31, 32]}
            Poses with no motif contacts are filtered out (marked as off-target)
        apply_offset : bool
            If True, apply residue number corrections for chopped AlphaFold models
            If False, use raw residue numbers from PDB (for full-length models)
        """
        self.distance_matrix, self.contact_sets, self.contact_residue_pairs, \
            pdb_files_filtered, filtered_indices = \
            DistanceMetrics.compute_jaccard_matrix(
                self.pdb_files,
                distance_cutoff=distance_cutoff,
                protein_selection=protein_selection,
                nucleic_selection=nucleic_selection,
                n_jobs=n_jobs,
                motif_residues=motif_residues,
                reference_sequence=self.reference_sequence if (apply_offset and self.reference_sequence) else None  # Only pass reference if offset enabled AND available
            )
        
        # Retrieve motif scores if motif was provided
        if motif_residues:
            self.motif_match_scores = getattr(DistanceMetrics, '_last_motif_match_scores', None)
        else:
            self.motif_match_scores = None
        
        # Update pdb_files to only include structures that passed filtering
        if len(filtered_indices) > 0:
            print(f"\n⚠️  {len(filtered_indices)} structure(s) filtered out (no contacts)")
            print(f"   Continuing analysis with {len(pdb_files_filtered)} valid structures")
            self.pdb_files = pdb_files_filtered
        
        # Store filtered indices regardless
        self.filtered_indices = filtered_indices
        
        # Store Jaccard matrix separately for later use
        self.jaccard_distance_matrix = self.distance_matrix.copy()
        
        self.metric_name = "Jaccard Contact Score"
        self.metric_unit = "(normalized)"
        
        print(f"✓ Contact lists stored (use export_contacts() to save)")
    
    def filter_duplicates(self, threshold=0.0001):
        """
        Filter out near-duplicate structures before clustering
        
        Identifies groups of structures with distance < threshold and keeps only one
        representative from each group. After clustering, duplicates are assigned
        to their representative's cluster.
        
        Parameters:
        -----------
        threshold : float
            Distance threshold below which structures are considered duplicates
            Default: 0.001 (very similar structures)
        
        Returns:
        --------
        reduced_matrix : np.ndarray
            Distance matrix with only representative structures
        representative_map : dict
            Maps original index -> representative index
        duplicate_groups : dict
            Maps representative index -> list of duplicate indices
        """
        if self.distance_matrix is None:
            raise ValueError("No distance matrix available. Run compute_distance_matrix() first.")
        
        n = self.distance_matrix.shape[0]
        
        # Track which structures are representatives and which are duplicates
        representative_indices = []  # Indices of structures to keep
        duplicate_groups = {}  # Maps representative -> list of duplicates
        assigned = set()  # Track which structures have been assigned
        
        print(f"\nFiltering duplicates (threshold={threshold})...")
        
        for i in range(n):
            if i in assigned:
                continue
            
            # This structure becomes a representative
            representative_indices.append(i)
            duplicates = [i]  # Include self
            
            # Find all duplicates of this structure
            for j in range(i + 1, n):
                if j in assigned:
                    continue
                
                if self.distance_matrix[i, j] < threshold:
                    duplicates.append(j)
                    assigned.add(j)
            
            duplicate_groups[i] = duplicates
            assigned.add(i)
        
        # Create mapping from original index to representative
        representative_map = {}
        for rep_idx, duplicates in duplicate_groups.items():
            for dup_idx in duplicates:
                representative_map[dup_idx] = rep_idx
        
        # Create reduced distance matrix with only representatives
        reduced_matrix = self.distance_matrix[np.ix_(representative_indices, representative_indices)]
        
        n_duplicates = n - len(representative_indices)
        print(f"Found {len(representative_indices)} representative structures")
        print(f"Filtered out {n_duplicates} near-duplicates ({n_duplicates/n*100:.1f}%)")
        
        # Show largest duplicate groups
        group_sizes = [(rep, len(dups)) for rep, dups in duplicate_groups.items()]
        group_sizes.sort(key=lambda x: x[1], reverse=True)
        print("\nLargest duplicate groups:")
        for rep, size in group_sizes[:5]:
            if size > 1:
                print(f"  Representative {rep}: {size} structures ({size-1} duplicates)")
        
        return reduced_matrix, representative_map, duplicate_groups, representative_indices
    
    def cluster_hdbscan(self, min_cluster_size=5, min_samples=2, filter_duplicates=True, duplicate_threshold=0.0001):
        """
        Perform HDBSCAN clustering (Hierarchical DBSCAN)
        
        HDBSCAN automatically finds optimal clustering without requiring eps parameter.
        Better for datasets with varying density clusters.
        
        Parameters:
        -----------
        min_cluster_size : int
            Minimum number of samples in a cluster (default: 5)
        min_samples : int
            Conservative noise threshold - higher values = more points labeled as noise (default: 2)
        filter_duplicates : bool
            If True, filter out near-duplicate structures before clustering (default: True)
        duplicate_threshold : float
            Distance threshold for identifying duplicates (default: 0.0001)
        """
        try:
            import hdbscan
        except ImportError:
            raise ImportError("hdbscan not installed. Install with: pip install hdbscan")
        
        print(f"\nPerforming HDBSCAN clustering (min_cluster_size={min_cluster_size}, min_samples={min_samples})...")
        
        if self.distance_matrix is None:
            raise ValueError("No distance matrix available. Run compute_distance_matrix() or compute_jaccard_contact_matrix() first.")
        
        # Optional: Filter duplicates
        if filter_duplicates:
            reduced_matrix, representative_map, duplicate_groups, representative_indices = self.filter_duplicates(duplicate_threshold)
            matrix_to_cluster = reduced_matrix
            self.representative_map = representative_map
            self.duplicate_groups = duplicate_groups
            self.representative_indices = representative_indices
        else:
            matrix_to_cluster = self.distance_matrix
            self.representative_map = None
            self.duplicate_groups = None
            self.representative_indices = None
        
        # Use the matrix as-is (duplicate filtering already handles near-zero distances)
        print(f"DEBUG: Distance range: {matrix_to_cluster[matrix_to_cluster > 0].min():.6e} - {matrix_to_cluster.max():.6f}", flush=True)
        
        # Convert to pure numpy array (in case it's a DataFrame or has pandas wrappers)
        # This ensures HDBSCAN's internal calculations work with clean numpy arrays
        import numpy as np
        if hasattr(matrix_to_cluster, 'to_numpy'):
            matrix_to_cluster = matrix_to_cluster.to_numpy()
        elif not isinstance(matrix_to_cluster, np.ndarray):
            matrix_to_cluster = np.asarray(matrix_to_cluster)
        print(f"DEBUG: Matrix type after conversion: {type(matrix_to_cluster)}", flush=True)
        
        # HDBSCAN clustering with precomputed distance matrix
        # gen_min_span_tree=True is required to compute DBCV (relative_validity_)
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric='precomputed',
            cluster_selection_method='eom',  # Excess of Mass - better for varying density
            gen_min_span_tree=True  # Required for DBCV calculation
        )
        reduced_labels = clusterer.fit_predict(matrix_to_cluster)
        labels_for_tree = reduced_labels  # This is what matches the condensed tree structure
        
        # Store the clusterer IMMEDIATELY after fitting, before any modifications
        self.hdbscan_clusterer = clusterer
        
        # If we filtered duplicates, map them back to full dataset
        if filter_duplicates:
            print(f"\nMapping {len(self.pdb_files) - len(representative_indices)} duplicates back to clusters...")
            full_labels = np.full(len(self.pdb_files), -1, dtype=int)
            
            # First pass: map based on representative's cluster
            for original_idx in range(len(self.pdb_files)):
                rep_idx = representative_map[original_idx]
                rep_position = representative_indices.index(rep_idx)
                full_labels[original_idx] = reduced_labels[rep_position]
            
            # Second pass: rescue duplicate groups whose representative is noise
            # These should form their own cluster instead of all being marked as noise
            next_cluster_id = max(reduced_labels) + 1 if len(reduced_labels) > 0 else 0
            rescued_clusters = 0
            
            for rep_idx, duplicates in duplicate_groups.items():
                rep_position = representative_indices.index(rep_idx)
                rep_label = reduced_labels[rep_position]
                
                # If representative is noise but has duplicates, form a new cluster
                if rep_label == -1 and len(duplicates) > 1:
                    # Assign all duplicates (including representative) to new cluster
                    for dup_idx in duplicates:
                        full_labels[dup_idx] = next_cluster_id
                    rescued_clusters += 1
                    next_cluster_id += 1
            
            if rescued_clusters > 0:
                print(f"✓ Rescued {rescued_clusters} duplicate groups from noise → formed new clusters")
            
            self.labels = full_labels  # Full labels for all structures (including duplicates)
            self.labels_reduced = reduced_labels  # Reduced labels matching the clustered matrix
            
            # Store which clusters are rescued for visualization
            self.rescued_cluster_ids = set(range(max(reduced_labels) + 1, next_cluster_id))
        else:
            self.labels = reduced_labels
            self.labels_reduced = None  # No reduction happened
            self.rescued_cluster_ids = set()  # No rescued clusters
        
        # Statistics
        n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise = list(self.labels).count(-1)
        
        print(f"Number of clusters: {n_clusters}")
        print(f"Number of noise points: {n_noise}")
        
        for i in range(n_clusters):
            cluster_size = list(self.labels).count(i)
            print(f"  Cluster {i}: {cluster_size} structures")
        
        # Debug: Print HDBSCAN internal metrics
        import sys
        print("\n=== HDBSCAN Clustering Quality ===", flush=True)
        sys.stdout.flush()
        
        # GLOSH (Global-Local Outlier Score from Hierarchies) is the quality metric
        # Lower scores = more inlier-like = better clustering quality
        if hasattr(clusterer, 'outlier_scores_'):
            outlier_scores = clusterer.outlier_scores_
            finite_scores = outlier_scores[np.isfinite(outlier_scores)]
            if len(finite_scores) > 0:
                median_glosh = np.median(finite_scores)
                print(f"✓ Median GLOSH Score: {median_glosh:.3f}", flush=True)
                print(f"  Range: {finite_scores.min():.3f} - {finite_scores.max():.3f}", flush=True)
                print(f"  (Lower = better clustering quality)", flush=True)
                sys.stdout.flush()
            else:
                print("⚠️  No finite GLOSH values", flush=True)
        else:
            print("⚠️  GLOSH outlier scores not available!", flush=True)
        
        print("=" * 40, flush=True)
        sys.stdout.flush()
        
        if hasattr(clusterer, 'condensed_tree_'):
            try:
                tree_df = clusterer.condensed_tree_.to_pandas()
                print(f"\nCondensed tree shape: {tree_df.shape}")
                print("Condensed tree sample (first 10 rows):")
                print(tree_df.head(10))
                
                # Check for infinite values
                n_infinite = np.isinf(tree_df['lambda_val']).sum()
                n_finite = np.isfinite(tree_df['lambda_val']).sum()
                print(f"\nLambda statistics:")
                print(f"  Finite values: {n_finite} ({n_finite/len(tree_df)*100:.1f}%)")
                print(f"  Infinite values: {n_infinite} ({n_infinite/len(tree_df)*100:.1f}%)")
                
                if n_finite > 0:
                    finite_lambdas = tree_df['lambda_val'][np.isfinite(tree_df['lambda_val'])]
                    print(f"  Finite range: min={finite_lambdas.min():.6f}, max={finite_lambdas.max():.6f}")
                    print(f"  Unique finite lambda values: {len(finite_lambdas.unique())}")
                
                # Show rows with infinite lambda
                if n_infinite > 0:
                    print(f"\nRows with infinite lambda (first 5):")
                    print(tree_df[np.isinf(tree_df['lambda_val'])].head())
            except Exception as e:
                print(f"Could not extract condensed tree data: {e}")
        
        if hasattr(clusterer, 'minimum_spanning_tree_'):
            try:
                mst = clusterer.minimum_spanning_tree_
                print(f"\nMinimum Spanning Tree edges: {mst.shape if hasattr(mst, 'shape') else 'N/A'}")
                mst_df = mst.to_pandas()
                print("MST sample (first 10 edges):")
                print(mst_df.head(10))
                print(f"MST distance range: min={mst_df['distance'].min():.6f}, max={mst_df['distance'].max():.6f}")
            except Exception as e:
                print(f"Could not extract MST data: {e}")
        
        print("=" * 40)
        
        return self.labels, labels_for_tree # was self.labels
    
    def cluster_dbscan(self, eps=2.0, min_samples=2):
        """
        Perform DBSCAN clustering (DEPRECATED - use cluster_hdbscan instead)
        
        Parameters:
        -----------
        eps : float
            Maximum distance between two samples
        min_samples : int
            Minimum number of samples in a neighborhood
        """
        print(f"\nPerforming DBSCAN clustering (eps={eps}, min_samples={min_samples})...")
        
        if self.distance_matrix is None:
            raise ValueError("No distance matrix available. Run compute_distance_matrix() or compute_jaccard_contact_matrix() first.")
        
        # DBSCAN clustering
        dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
        self.labels = dbscan.fit_predict(self.distance_matrix)
        
        # Statistics
        n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise = list(self.labels).count(-1)
        
        print(f"Number of clusters: {n_clusters}")
        print(f"Number of noise points: {n_noise}")
        
        for i in range(n_clusters):
            cluster_size = list(self.labels).count(i)
            print(f"  Cluster {i}: {cluster_size} structures")
        
        return self.labels
    
    def cluster_hierarchical(self, n_clusters=5, method='average', metric='precomputed'):
        """
        Perform hierarchical clustering
        
        Parameters:
        -----------
        n_clusters : int
            Number of clusters to form
        method : str
            Linkage method ('average', 'single', 'complete', 'ward')
        metric : str
            Distance metric ('precomputed' uses the distance matrix)
        """
        print(f"\nPerforming hierarchical clustering ({method} linkage, {n_clusters} clusters)...")
        
        if self.distance_matrix is None:
            raise ValueError("No distance matrix available")
        
        # Convert distance matrix to condensed form for linkage
        from scipy.spatial.distance import squareform
        condensed_dist = squareform(self.distance_matrix)
        
        # Perform linkage
        Z = linkage(condensed_dist, method=method)
        
        # Form flat clusters
        self.labels = fcluster(Z, n_clusters, criterion='maxclust') - 1  # 0-indexed
        
        print(f"Number of clusters: {n_clusters}")
        for i in range(n_clusters):
            cluster_size = list(self.labels).count(i)
            print(f"  Cluster {i}: {cluster_size} structures")
        
        return self.labels
    
    def get_interface_stats(self):
        """Get interface statistics (requires Jaccard computation)"""
        if self.contact_residue_pairs is None:
            raise ValueError("No contact data. Run compute_jaccard_contact_matrix() first.")
        
        return InterfaceAnalyzer.get_interface_stats(
            self.pdb_files,
            self.contact_residue_pairs,
            self.labels,
            self.representative_indices
        )
    
    def get_binding_hotspots(self, cluster_id=None):
        """Get binding hotspots (requires Jaccard computation)"""
        if self.contact_residue_pairs is None:
            raise ValueError("No contact data. Run compute_jaccard_contact_matrix() first.")
        
        return InterfaceAnalyzer.get_binding_hotspots(
            self.contact_residue_pairs,
            self.labels,
            cluster_id
        )
    
    def get_cluster_signature(self, cluster_id, threshold=50):
        """Get consensus interface for a cluster"""
        if self.contact_residue_pairs is None:
            raise ValueError("No contact data. Run compute_jaccard_contact_matrix() first.")
        
        return InterfaceAnalyzer.get_cluster_signature(
            self.contact_residue_pairs,
            self.labels,
            cluster_id,
            threshold
        )
    
    def get_cluster_summary(self, threshold=50):
        """Get summary table of all clusters with consensus residues"""
        if self.contact_residue_pairs is None:
            raise ValueError("No contact data. Run compute_jaccard_contact_matrix() first.")
        if self.labels is None:
            raise ValueError("No clustering results. Run cluster_dbscan() or cluster_hierarchical() first.")
        
        return InterfaceAnalyzer.get_cluster_summary(
            self.contact_residue_pairs,
            self.labels,
            threshold
        )
