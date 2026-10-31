"""
Main clustering class for PDB structures
"""

import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
from sklearn.cluster import DBSCAN
from scipy.cluster.hierarchy import linkage, fcluster
import tempfile

from .metrics import DistanceMetrics
from .analysis import InterfaceAnalyzer


def load_universe_with_cif_support(pdb_file):
    """
    Load MDAnalysis Universe with support for both PDB and CIF formats.
    CIF files are converted to PDB on-the-fly using BioPython.
    
    Parameters:
    -----------
    pdb_file : str or Path
        Path to structure file (.pdb or .cif)
    
    Returns:
    --------
    MDAnalysis.Universe
    """
    import MDAnalysis as mda
    
    file_ext = str(pdb_file).lower()
    if file_ext.endswith('.cif'):
        # AlphaFold 3 CIF files - convert to PDB using BioPython
        from Bio.PDB import MMCIFParser, PDBIO
        import atexit
        
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('structure', pdb_file)
        
        # Write to temporary PDB file - keep it persistent for MDAnalysis DSSP
        # Use delete=True so system cleans it up, but don't unlink immediately
        temp_pdb = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        temp_file_path = temp_pdb.name
        temp_pdb.close()
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(temp_file_path)
        
        # Load the converted PDB
        u = mda.Universe(temp_file_path)
        
        # Register cleanup to happen at program exit (not immediately)
        # This allows MDAnalysis DSSP to reopen the file when needed
        def cleanup_temp():
            try:
                if os.path.exists(temp_file_path):
                    os.unlink(temp_file_path)
            except:
                pass
        atexit.register(cleanup_temp)
        
        return u
    else:
        # Standard PDB files
        return mda.Universe(pdb_file)


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
        self.svd_coordinates = None
        self.svd_loadings = None
        self.svd_model_coordinates = None
        self.svd_feature_names = []
        self.svd_explained_variance_ratio = None
        self.svd_elbow_components = None
        self.svd_elbow_variance_ratio = None
        self.svd_elbow_data = None
        self.svd_cluster_elbow_data = None
        self.svd_grid_search_data = None
        self.svd_optimized_components = None
        self.svd_optimized_k = None
        self.svd_cluster_count_target = None
        self.svd_cluster_count_preferred_range = None
        self.svd_optimization_threshold_met = None
        self.svd_isolated_points = None
        self.snn_network = None
        
        # Residue offset correction
        self.reference_sequence = None  # Full-length reference sequence
        
        print(f"Found {self.n_structures} PDB files")
    
    def load_structures(self):
        """Load all PDB structures using MDAnalysis (supports PDB, CIF/mmCIF via BioPython conversion)"""
        print("Loading structures...")
        self.universes = []
        
        for pdb_file in self.pdb_files:
            try:
                u = load_universe_with_cif_support(pdb_file)
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

    def _project_distance_matrix_2d(self, distance_matrix, method='mds'):
        """
        Project a precomputed pairwise distance matrix into 2D coordinates.

        Parameters:
        -----------
        distance_matrix : np.ndarray
            Pairwise distance matrix (n x n)
        method : str
            Projection method: 'mds', 'tsne', 'umap', or 'pca'

        Returns:
        --------
        np.ndarray
            2D coordinates (n x 2)
        """
        method = str(method).strip().lower()
        if method not in {'mds', 'tsne', 'umap', 'pca'}:
            raise ValueError("Projection method must be one of: mds, tsne, umap, pca")

        matrix = np.asarray(distance_matrix, dtype=float)
        n_samples = matrix.shape[0]

        # Tiny datasets cannot be stably projected with nonlinear methods.
        if n_samples <= 2:
            return np.column_stack([np.arange(n_samples, dtype=float), np.zeros(n_samples, dtype=float)])

        if method == 'mds':
            from sklearn.manifold import MDS
            reducer = MDS(
                n_components=2,
                dissimilarity='precomputed',
                random_state=42,
                n_init=1,
                max_iter=300
            )
            return reducer.fit_transform(matrix)

        if method == 'tsne':
            from sklearn.manifold import TSNE
            perplexity = min(30, max(1, n_samples // 3))
            if perplexity >= n_samples:
                perplexity = max(1, n_samples - 1)
            reducer = TSNE(
                n_components=2,
                metric='precomputed',
                init='random',
                random_state=42,
                perplexity=perplexity,
                learning_rate='auto'
            )
            return reducer.fit_transform(matrix)

        if method == 'umap':
            try:
                from umap import UMAP
            except ImportError:
                raise ImportError("UMAP clustering requested, but umap-learn is not installed. Install with: pip install umap-learn")

            n_neighbors = min(15, max(2, n_samples - 1))
            reducer = UMAP(
                n_components=2,
                metric='precomputed',
                random_state=42,
                n_neighbors=n_neighbors
            )
            return reducer.fit_transform(matrix)

        # PCA on distance profiles (rows of the pairwise distance matrix).
        from sklearn.decomposition import PCA
        reducer = PCA(n_components=2)
        return reducer.fit_transform(matrix)
    
    def cluster_hdbscan(self, min_cluster_size=5, min_samples=2,
                        filter_duplicates=True, duplicate_threshold=0.0001,
                        cluster_selection_epsilon=0.0, cluster_selection_method='eom'):
        """
        Perform HDBSCAN clustering (Hierarchical DBSCAN)
        
        HDBSCAN automatically finds optimal clustering and supports
        optional epsilon-based cluster merging plus EOM/Leaf selection.
        
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
        cluster_selection_epsilon : float
            Merge clusters within this distance threshold (default: 0.0)
        cluster_selection_method : str
            Clustering strategy:
            - 'eom' or 'leaf': cluster on raw precomputed distance matrix
            - 'mds', 'tsne', 'umap', 'pca': project to 2D first, then cluster with EOM
        """
        try:
            import hdbscan
        except ImportError:
            raise ImportError("hdbscan not installed. Install with: pip install hdbscan")

        self.svd_elbow_components = None
        self.svd_elbow_variance_ratio = None
        self.svd_elbow_data = None
        self.svd_cluster_elbow_data = None
        self.svd_grid_search_data = None
        self.svd_optimized_components = None
        self.svd_optimized_k = None
        self.svd_optimization_threshold_met = None

        cluster_selection_method = str(cluster_selection_method).strip().lower()
        projection_methods = {'mds', 'tsne', 'umap', 'pca'}

        if cluster_selection_method in {'eom', 'leaf'}:
            projection_method = None
            hdbscan_selection_method = cluster_selection_method
        elif cluster_selection_method in projection_methods:
            projection_method = cluster_selection_method
            hdbscan_selection_method = 'eom'  # Default as requested for 2D methods
        else:
            raise ValueError(
                "cluster_selection_method must be one of: "
                "eom, leaf, mds, tsne, umap, pca"
            )

        if cluster_selection_epsilon is None:
            cluster_selection_epsilon = 0.0
        cluster_selection_epsilon = float(cluster_selection_epsilon)
        if cluster_selection_epsilon < 0:
            raise ValueError("cluster_selection_epsilon must be >= 0")
        
        print(
            f"\nPerforming HDBSCAN clustering "
            f"(min_cluster_size={min_cluster_size}, min_samples={min_samples}, "
            f"epsilon={cluster_selection_epsilon}, method={cluster_selection_method})..."
        )
        
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
        
        # Convert to pure numpy array (in case it's a DataFrame or has pandas wrappers)
        # This ensures HDBSCAN's internal calculations work with clean numpy arrays
        if hasattr(matrix_to_cluster, 'to_numpy'):
            matrix_to_cluster = matrix_to_cluster.to_numpy()
        elif not isinstance(matrix_to_cluster, np.ndarray):
            matrix_to_cluster = np.asarray(matrix_to_cluster)

        # Use the matrix as-is (duplicate filtering already handles near-zero distances)
        non_zero = matrix_to_cluster[matrix_to_cluster > 0]
        if len(non_zero) > 0:
            print(f"DEBUG: Distance range: {non_zero.min():.6e} - {matrix_to_cluster.max():.6f}", flush=True)
        else:
            print("DEBUG: Distance range: all pairwise distances are zero", flush=True)
        print(f"DEBUG: Matrix type after conversion: {type(matrix_to_cluster)}", flush=True)

        # Optional: project pairwise distances to 2D before clustering.
        data_to_cluster = matrix_to_cluster
        hdbscan_metric = 'precomputed'
        if projection_method is not None:
            print(
                f"Using {projection_method.upper()} 2D embedding for clustering "
                f"(HDBSCAN cluster_selection_method defaults to '{hdbscan_selection_method}')",
                flush=True
            )
            data_to_cluster = self._project_distance_matrix_2d(matrix_to_cluster, method=projection_method)
            hdbscan_metric = 'euclidean'
            self.clustering_input_mode = f"{projection_method}_2d"
            self.clustering_projection_2d = data_to_cluster
        else:
            self.clustering_input_mode = 'precomputed_distance'
            self.clustering_projection_2d = None
        
        # HDBSCAN clustering using either precomputed distances or 2D coordinates.
        # gen_min_span_tree=True is required to compute DBCV (relative_validity_)
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric=hdbscan_metric,
            cluster_selection_method=hdbscan_selection_method,
            cluster_selection_epsilon=cluster_selection_epsilon,
            gen_min_span_tree=True  # Required for DBCV calculation
        )
        reduced_labels = clusterer.fit_predict(data_to_cluster)
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

    def cluster_svd_leiden(self, max_components=100, n_neighbors=None,
                           n_components=None, snn_epsilon=0.0,
                           prune_edges=False, min_shared_neighbors=1,
                           optimize_grid=False, optimize_k=False,
                           max_average_intra_jaccard=0.6):
        """Cluster contact fingerprints with SVD coordinates and SNN-Leiden.

        The contact matrix is built from the protein residues observed in each
        structure. SVD provides the continuous representation, and a shared
        nearest-neighbor graph is partitioned with Leiden. Results needed for
        inspection and export are retained on this clusterer instance.
        """
        if self.contact_residue_pairs is None:
            raise ValueError("No contact data. Run compute_jaccard_contact_matrix() first.")

        n_structures = len(self.contact_residue_pairs)
        if n_structures != len(self.pdb_files):
            raise ValueError(
                f"Contact data ({n_structures}) does not match structures ({len(self.pdb_files)})."
            )
        if n_structures == 0:
            raise ValueError("No structures available for SVD-Leiden clustering.")
        self.svd_grid_search_data = None
        self.svd_optimized_components = None
        self.svd_optimized_k = None
        coherence_matrix = getattr(self, 'jaccard_distance_matrix', None)
        if coherence_matrix is None:
            coherence_matrix = self.distance_matrix
        if coherence_matrix is None:
            raise ValueError("No distance matrix available for coherence scoring.")
        if max_average_intra_jaccard is not None and max_average_intra_jaccard < 0:
            raise ValueError("max_average_intra_jaccard must be >= 0 or None.")

        from scipy.sparse import csr_matrix, triu
        from sklearn.decomposition import TruncatedSVD
        from sklearn.neighbors import kneighbors_graph
        import igraph as ig
        import leidenalg as la

        # Build one binary row per structure and one feature per contacted
        # protein residue, preserving the structure order used by labels.
        contact_sets = []
        for contacts in self.contact_residue_pairs:
            contact_sets.append({contact['protein_residue'] for contact in contacts})

        feature_names = sorted(set().union(*contact_sets)) if contact_sets else []
        self.svd_feature_names = feature_names
        if not feature_names:
            raise ValueError("No protein-residue contact features were found.")

        feature_index = {feature: index for index, feature in enumerate(feature_names)}
        feature_matrix = np.zeros((n_structures, len(feature_names)), dtype=float)
        for row_index, structure_contacts in enumerate(contact_sets):
            for feature in structure_contacts:
                feature_matrix[row_index, feature_index[feature]] = 1.0

        sparse_features = csr_matrix(feature_matrix)
        component_limit = max(1, min(
            int(max_components), n_structures - 1, len(feature_names)
        ))
        if n_components is None or str(n_components).lower() == 'auto':
            n_components = component_limit
            svd_probe = TruncatedSVD(n_components=component_limit, random_state=42)
            svd_probe.fit(sparse_features)
            variance_curve = svd_probe.explained_variance_ratio_
            if len(variance_curve) > 2:
                points = np.column_stack((
                    np.arange(len(variance_curve)), variance_curve
                ))
                line = points[-1] - points[0]
                line_norm = line / max(np.linalg.norm(line), 1e-12)
                offsets = points - points[0]
                projection = np.outer(offsets @ line_norm, line_norm)
                elbow_distances = np.linalg.norm(offsets - projection, axis=1)
                n_components = int(np.argmax(elbow_distances) + 1)
            # A one- or two-component solution is usually dominated by the
            # first variance spike; retain the practical three-component knee
            # when the dataset supports it.
            n_components = max(1, min(component_limit, max(3, n_components)))
            self.svd_elbow_components = n_components
            self.svd_elbow_variance_ratio = variance_curve
        else:
            n_components = max(1, min(int(n_components), component_limit))
            svd_probe = TruncatedSVD(n_components=component_limit, random_state=42)
            svd_probe.fit(sparse_features)
            self.svd_elbow_components = n_components
            self.svd_elbow_variance_ratio = svd_probe.explained_variance_ratio_
        self.svd_elbow_data = pd.DataFrame({
            'component': np.arange(1, len(self.svd_elbow_variance_ratio) + 1),
            'explained_variance_ratio': self.svd_elbow_variance_ratio,
            'selected': np.arange(1, len(self.svd_elbow_variance_ratio) + 1) == n_components
        })

        svd_model = TruncatedSVD(n_components=n_components, random_state=42)
        self.svd_coordinates = svd_model.fit_transform(sparse_features)
        self.svd_explained_variance_ratio = svd_model.explained_variance_ratio_
        coordinate_columns = [
            f'component_{index + 1}'
            for index in range(self.svd_coordinates.shape[1])
        ]
        self.svd_model_coordinates = pd.DataFrame(
            self.svd_coordinates,
            columns=coordinate_columns
        )
        self.svd_model_coordinates.insert(0, 'structure', [
            os.path.basename(path) for path in self.pdb_files
        ])
        self.svd_model_coordinates.insert(
            0, 'model_index', np.arange(n_structures)
        )

        loading_matrix = svd_model.components_.T * svd_model.singular_values_
        loading_columns = [
            f'component_{index + 1}_loading'
            for index in range(loading_matrix.shape[1])
        ]
        loading_df = pd.DataFrame(loading_matrix, columns=loading_columns)
        loading_df.insert(0, 'residue', feature_names)
        loading_df['max_absolute_loading'] = loading_df[loading_columns].abs().max(axis=1)
        loading_df['dominant_component'] = (
            loading_df[loading_columns].abs().idxmax(axis=1)
            .str.replace('_loading', '', regex=False)
        )
        self.svd_loadings = loading_df.sort_values(
            'max_absolute_loading', ascending=False
        ).reset_index(drop=True)

        # Estimate a useful cluster-count range before the grid search. This
        # keeps the elbow as a soft validation constraint while graph quality
        # still decides between otherwise comparable candidates.
        from sklearn.cluster import KMeans
        # Scale the cluster-count diagnostic with dataset size so larger
        # runs do not hide the selected cluster count beyond an arbitrary 20.
        cluster_probe_max = min(
            n_structures - 1,
            max(20, min(100, int(np.ceil(2 * np.sqrt(n_structures)))))
        )
        cluster_probe_rows = []
        if (optimize_grid or optimize_k) and n_structures > 3:
            for cluster_count in range(1, cluster_probe_max + 1):
                kmeans = KMeans(n_clusters=cluster_count, random_state=42, n_init=10)
                kmeans.fit(self.svd_coordinates)
                cluster_probe_rows.append((cluster_count, float(kmeans.inertia_)))
        if len(cluster_probe_rows) >= 3:
            probe_points = np.asarray(cluster_probe_rows, dtype=float)
            probe_line = probe_points[-1] - probe_points[0]
            probe_line_norm = probe_line / max(np.linalg.norm(probe_line), 1e-12)
            probe_offsets = probe_points - probe_points[0]
            probe_projection = np.outer(
                probe_offsets @ probe_line_norm, probe_line_norm
            )
            probe_distances = np.linalg.norm(
                probe_offsets - probe_projection, axis=1
            )
            cluster_count_target = max(2, int(
                cluster_probe_rows[int(np.argmax(probe_distances))][0]
            ))
        else:
            cluster_count_target = None
        cluster_count_min = (
            max(2, cluster_count_target - 1)
            if cluster_count_target is not None else None
        )
        cluster_count_max = (
            min(cluster_probe_max, cluster_count_target + 1)
            if cluster_count_target is not None else None
        )

        grid_optimized = False
        if (optimize_grid or optimize_k) and n_structures > 3:
            if optimize_grid:
                candidate_components = sorted(set(
                    max(1, min(component_limit, n_components + offset))
                    for offset in range(-5, 6)
                ))
            else:
                candidate_components = [n_components]
            coarse_k_min = 3
            coarse_k_max = min(30, n_structures - 1)
            # Use a contiguous k grid so the heatmap represents an actual
            # rectangular search surface rather than a sparse coarse/fine mix.
            coarse_k_values = list(range(coarse_k_min, coarse_k_max + 1))
            if not coarse_k_values:
                coarse_k_values = [min(3, n_structures - 1)]

            def evaluate_candidate(component_count, candidate_k):
                candidate_model = TruncatedSVD(
                    n_components=component_count, random_state=42
                )
                candidate_coordinates = candidate_model.fit_transform(sparse_features)
                effective_k = min(max(1, candidate_k), n_structures - 1)
                candidate_knn = kneighbors_graph(
                    candidate_coordinates,
                    n_neighbors=effective_k,
                    mode='connectivity',
                    include_self=True
                )
                candidate_snn = candidate_knn.dot(candidate_knn.T).tocsr()
                if float(snn_epsilon) < 0:
                    raise ValueError("snn_epsilon must be >= 0")
                if prune_edges:
                    minimum_weight = max(float(snn_epsilon), float(min_shared_neighbors))
                    candidate_snn.data[candidate_snn.data < minimum_weight] = 0
                    candidate_snn.eliminate_zeros()
                candidate_snn = candidate_snn.tolil()
                candidate_snn.setdiag(0)
                candidate_snn = candidate_snn.tocsr()
                upper_candidate_snn = triu(candidate_snn, k=1).tocsr()
                sources, targets = upper_candidate_snn.nonzero()
                graph = ig.Graph(n=n_structures, directed=False)
                if len(sources) == 0:
                    labels = np.arange(n_structures, dtype=int)
                    modularity = 0.0
                else:
                    graph.add_edges(
                        list(zip(sources.tolist(), targets.tolist())),
                        attributes={'weight': upper_candidate_snn.data.tolist()}
                    )
                    partition = la.find_partition(
                        graph, la.ModularityVertexPartition,
                        weights='weight', seed=42
                    )
                    labels = np.asarray(partition.membership, dtype=int)
                    modularity = float(partition.modularity)

                cluster_distances = []
                for cluster_id in np.unique(labels):
                    members = np.flatnonzero(labels == cluster_id)
                    if len(members) > 1:
                        distances = coherence_matrix[np.ix_(members, members)]
                        cluster_distances.extend(
                            distances[np.triu_indices(len(members), k=1)].tolist()
                        )
                average_intra_jaccard = (
                    float(np.mean(cluster_distances))
                    if cluster_distances else np.inf
                )
                cluster_sizes = np.bincount(labels)
                return {
                    'n_components': component_count,
                    'k': effective_k,
                    'n_clusters': int(len(np.unique(labels))),
                    'average_intra_jaccard': average_intra_jaccard,
                    'modularity': modularity,
                    'largest_cluster_fraction': float(cluster_sizes.max() / n_structures),
                    'valid': bool(
                        len(np.unique(labels)) >= 2 and len(cluster_distances) > 0
                    ),
                    'stability': np.nan,
                    'labels': labels,
                    'coordinates': candidate_coordinates
                }

            grid_rows = []
            for component_count in candidate_components:
                for candidate_k in coarse_k_values:
                    grid_rows.append(evaluate_candidate(component_count, candidate_k))
            valid_rows = [row for row in grid_rows if row['valid']]
            if valid_rows:
                from sklearn.metrics import adjusted_rand_score
                for row in grid_rows:
                    neighboring_rows = [
                        other for other in grid_rows
                        if other is not row and
                        abs(other['n_components'] - row['n_components']) +
                        abs(other['k'] - row['k']) <= 2
                    ]
                    row['stability'] = (
                        float(np.mean([
                            adjusted_rand_score(row['labels'], other['labels'])
                            for other in neighboring_rows
                        ])) if neighboring_rows else 0.0
                    )
                valid_rows = [row for row in grid_rows if row['valid']]
                threshold_rows = [
                    row for row in valid_rows
                    if max_average_intra_jaccard is None or
                    row['average_intra_jaccard'] <= max_average_intra_jaccard
                ]
                selection_pool = threshold_rows or valid_rows
                elbow_rows = [
                    row for row in selection_pool
                    if cluster_count_min is not None and
                    cluster_count_min <= row['n_clusters'] <= cluster_count_max
                ]
                if elbow_rows:
                    selection_pool = elbow_rows
                # The heatmap's primary metric must also drive selection.
                # Modularity and stability are useful tie-breakers, not a
                # reason to choose a visibly worse intra-cluster distance.
                best_row = min(
                    selection_pool,
                    key=lambda row: (
                        row['average_intra_jaccard'],
                        -row['stability'],
                        -row['modularity'],
                        row['largest_cluster_fraction'],
                        abs(row['n_clusters'] - (cluster_count_target or row['n_clusters']))
                    )
                )
                n_components = best_row['n_components']
                neighbor_count = best_row['k']
                self.svd_grid_search_data = pd.DataFrame([
                    {key: value for key, value in row.items()
                     if key not in {'labels', 'coordinates'}}
                    for row in grid_rows
                ]).drop_duplicates(subset=['n_components', 'k'])
                self.svd_grid_search_data['selected'] = (
                    (self.svd_grid_search_data['n_components'] == n_components) &
                    (self.svd_grid_search_data['k'] == neighbor_count)
                )
                self.svd_grid_search_data['cluster_count_target'] = cluster_count_target
                self.svd_grid_search_data['cluster_count_preferred'] = (
                    self.svd_grid_search_data['n_clusters'].between(
                        cluster_count_min, cluster_count_max
                    ) if cluster_count_target is not None else False
                )
                self.svd_cluster_count_target = cluster_count_target
                self.svd_cluster_count_preferred_range = (
                    (cluster_count_min, cluster_count_max)
                    if cluster_count_target is not None else None
                )
                self.svd_coordinates = best_row['coordinates']
                self.svd_optimized_components = n_components
                self.svd_optimized_k = neighbor_count
                self.svd_optimization_threshold_met = bool(threshold_rows)
                grid_optimized = True
                print(
                    f"Grid optimization selected components={n_components}, "
                    f"k={neighbor_count}, average intra-Jaccard="
                    f"{best_row['average_intra_jaccard']:.4f}"
                )

        if grid_optimized:
            # Refit the selected representation so coordinates and feature
            # loadings correspond to the optimized component count.
            svd_model = TruncatedSVD(n_components=n_components, random_state=42)
            self.svd_coordinates = svd_model.fit_transform(sparse_features)
            self.svd_explained_variance_ratio = svd_model.explained_variance_ratio_
            coordinate_columns = [
                f'component_{index + 1}'
                for index in range(self.svd_coordinates.shape[1])
            ]
            self.svd_model_coordinates = pd.DataFrame(
                self.svd_coordinates, columns=coordinate_columns
            )
            self.svd_model_coordinates.insert(0, 'structure', [
                os.path.basename(path) for path in self.pdb_files
            ])
            self.svd_model_coordinates.insert(0, 'model_index', np.arange(n_structures))
            loading_matrix = svd_model.components_.T * svd_model.singular_values_
            loading_columns = [
                f'component_{index + 1}_loading'
                for index in range(loading_matrix.shape[1])
            ]
            loading_df = pd.DataFrame(loading_matrix, columns=loading_columns)
            loading_df.insert(0, 'residue', feature_names)
            loading_df['max_absolute_loading'] = loading_df[loading_columns].abs().max(axis=1)
            loading_df['dominant_component'] = (
                loading_df[loading_columns].abs().idxmax(axis=1)
                .str.replace('_loading', '', regex=False)
            )
            self.svd_loadings = loading_df.sort_values(
                'max_absolute_loading', ascending=False
            ).reset_index(drop=True)

        if n_structures == 1:
            self.snn_network = csr_matrix((1, 1), dtype=float)
            self.svd_isolated_points = pd.DataFrame({
                'model_index': [0],
                'structure': [os.path.basename(self.pdb_files[0])],
                'network_degree': [0]
            })
            self.labels = np.array([0], dtype=int)
            self.metric_name = "SVD-SNN-Leiden"
            self.metric_unit = "contact features"
            return self.labels

        if grid_optimized:
            neighbor_count = self.svd_optimized_k
        elif n_neighbors is None or str(n_neighbors).lower() == 'auto':
            neighbor_count = min(max(2, int(np.sqrt(n_structures))), n_structures - 1)
        else:
            neighbor_count = min(max(1, int(n_neighbors)), n_structures - 1)
        # Conventional elbow diagnostic: vary the number of clusters and
        # measure KMeans inertia (WCSS) on the same SVD coordinates.
        from sklearn.cluster import KMeans
        cluster_sweep_max = cluster_probe_max
        cluster_elbow_rows = []
        for cluster_count in range(1, cluster_sweep_max + 1):
            kmeans = KMeans(n_clusters=cluster_count, random_state=42, n_init=10)
            kmeans.fit(self.svd_coordinates)
            cluster_elbow_rows.append({
                'n_clusters': cluster_count,
                'wcss': float(kmeans.inertia_),
                'wcss_per_point': float(kmeans.inertia_ / n_structures),
                'selected': False
            })
        self.svd_cluster_elbow_data = pd.DataFrame(cluster_elbow_rows)

        knn_graph = kneighbors_graph(
            self.svd_coordinates,
            n_neighbors=neighbor_count,
            mode='connectivity',
            include_self=True
        )
        snn_matrix_with_self = knn_graph.dot(knn_graph.T).tocsr()
        snn_matrix = snn_matrix_with_self.copy().tolil()
        snn_matrix.setdiag(0)
        snn_matrix = snn_matrix.tocsr()
        if float(snn_epsilon) < 0:
            raise ValueError("snn_epsilon must be >= 0")
        if prune_edges:
            minimum_weight = max(float(snn_epsilon), float(min_shared_neighbors))
            snn_network_data = snn_matrix.data
            snn_matrix.data[snn_network_data < minimum_weight] = 0
            snn_matrix.eliminate_zeros()
        network_degree = np.asarray((snn_matrix > 0).sum(axis=1)).ravel()
        isolated_indices = np.flatnonzero(network_degree == 0)
        self.snn_network = snn_matrix
        self.svd_isolated_points = pd.DataFrame({
            'model_index': isolated_indices,
            'structure': [os.path.basename(self.pdb_files[index]) for index in isolated_indices],
            'network_degree': network_degree[isolated_indices]
        })

        # Use one undirected edge per shared-neighbor pair for Leiden, matching
        # the graph representation used during grid evaluation.
        clean_snn = snn_matrix
        upper_snn = triu(clean_snn, k=1).tocsr()
        sources, targets = upper_snn.nonzero()
        weights = upper_snn.data

        if len(sources) == 0:
            cluster_labels = np.arange(n_structures, dtype=int)
        else:
            import igraph as ig
            import leidenalg as la

            graph = ig.Graph(n=n_structures, directed=False)
            graph.add_edges(
                list(zip(sources.tolist(), targets.tolist())),
                attributes={'weight': weights.tolist()}
            )
            partition = la.find_partition(
                graph,
                la.ModularityVertexPartition,
                weights='weight',
                seed=42
            )
            cluster_labels = np.asarray(partition.membership, dtype=int)

        self.labels = cluster_labels
        self.svd_cluster_elbow_data['selected'] = (
            self.svd_cluster_elbow_data['n_clusters'] == len(np.unique(cluster_labels))
        )
        self.labels_reduced = None
        self.hdbscan_clusterer = None
        self.representative_indices = None
        self.representative_map = None
        self.duplicate_groups = None
        self.rescued_cluster_ids = set()
        self.clustering_input_mode = "svd_snn"
        self.clustering_projection_2d = self.svd_coordinates
        self.svd_n_neighbors = neighbor_count
        self.svd_snn_epsilon = float(snn_epsilon)
        self.svd_prune_edges = bool(prune_edges)
        self.svd_min_shared_neighbors = int(min_shared_neighbors)
        self.metric_name = "SVD-SNN-Leiden"
        self.metric_unit = "contact features"

        print(
            f"SVD-SNN-Leiden detected {len(np.unique(cluster_labels))} clusters "
            f"from {n_structures} structures; {len(isolated_indices)} isolated points."
        )
        return self.labels
    
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
