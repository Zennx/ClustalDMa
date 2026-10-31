#Testing of new method of clustering

import pandas as pd
import numpy as np
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt
from sklearn.neighbors import kneighbors_graph
import umap
from scipy.sparse import csr_matrix
import igraph as ig
import leidenalg as la

# STEP 1: Load the data and data preprocessing
#load data (prototype reads CSV file)
df = pd.read_csv('/Users/Zen/Downloads/clustaldm_results_20260925_054426/cluster_statistics.csv', usecols=["structure", "protein_residues"])

print(df.head())

df['protein_residues'] = df['protein_residues'].astype(str).str.split(',').apply(lambda x: [i.strip() for i in x])  # Split the string into a list of residues

df_exploded = df.explode('protein_residues')  # Explode the list into separate rows

print(df_exploded.head())

#Cross-tabulation from two columns
ct = pd.crosstab(df_exploded['structure'], df_exploded['protein_residues'])
ct = (ct > 0).astype(int)

ct = ct.loc[:, (ct != 0).any(axis=0)]  # Remove columns with all zeros, theoretically this should not happen but just in case
ct = ct.loc[(ct != 0).any(axis=1), :]  # Remove rows with all zeros, theoretically this should not happen but just in case

print(ct)
ct.to_csv('/Users/Zen/Downloads/clustaldm_results_20260925_054426/crosstab_structure_protein_residues.csv', index=True)

residue_features = ct.columns.to_list()

ct = csr_matrix(ct.values)  # Convert to sparse matrix for efficiency


# STEP 2: Reduction of the data using SVD
# 2.1 Perform SVD on the cross-tabulation matrix and plot elbow plot
n_samples, n_features = ct.shape
theoretical_max_components = min(n_samples, n_features) -1

max_components = min(100, theoretical_max_components)  # Set the maximum number of components for SVD
svd_test = TruncatedSVD(n_components=max_components, random_state=42)
svd_test.fit(ct)

#automated n_components selection based on dynamic variance curve analysis
import numpy as np
from sklearn.decomposition import TruncatedSVD

def find_true_elbow(explained_variance_ratio):
    """Finds the maximum perpendicular distance from the line connecting start and end points."""
    n_points = len(explained_variance_ratio)
    coords = np.vstack((np.arange(n_points), explained_variance_ratio)).T
    
    # Line vector from first point to last point
    first_pt = coords[0]
    last_pt = coords[-1]
    line_vec = last_pt - first_pt
    line_vec_norm = line_vec / np.linalg.norm(line_vec)
    
    # Calculate perpendicular distance for all intermediate points
    vectors_to_first = coords - first_pt
    scalar_projections = np.dot(vectors_to_first, line_vec_norm)
    projected_vectors = np.outer(scalar_projections, line_vec_norm)
    perpendicular_vectors = vectors_to_first - projected_vectors
    distances = np.linalg.norm(perpendicular_vectors, axis=1)
    
    # The elbow is the point with the maximum distance
    return np.argmax(distances) + 1

def get_optimal_svd_coords_v2(sparse_contact_matrix):
    n_samples, n_features = sparse_contact_matrix.shape
    theoretical_max = min(n_samples, n_features) - 1
    
    if theoretical_max <= 2:
        return sparse_contact_matrix.toarray(), None
    
    # Fit the maximum reasonable components to capture the shape of the decay
    upper_bound = min(100, theoretical_max)
    svd_test = TruncatedSVD(n_components=upper_bound, random_state=42)
    svd_test.fit(sparse_contact_matrix)
    
    # Automatically locate the physical geometric elbow
    optimal_n = find_true_elbow(svd_test.explained_variance_ratio_)
    
    # Safety fallback: always retain at least 5 components if possible
    optimal_n = max(min(5, theoretical_max), optimal_n)
    
    print(f"Geometric Elbow Detector automatically selected: {optimal_n} components.")

    plt.figure(figsize=(10, 6))
    plt.plot(range(1, max_components + 1), svd_test.explained_variance_ratio_, marker='o')
    plt.axvline(x=optimal_n, color='r', linestyle='--', label=f'Selected: {optimal_n}')
    plt.legend()
    plt.title('Explained Variance Ratio by SVD Components')
    plt.xlabel('Number of Components')
    plt.ylabel('Explained Variance Ratio')
    plt.grid()
    plt.savefig('/Users/Zen/Downloads/clustaldm_results_20260925_054426/svd_explained_variance_ratio.png')
    plt.show()

    svd_final = TruncatedSVD(n_components=optimal_n, random_state=42)
    svd_coords = svd_final.fit_transform(sparse_contact_matrix)
    return svd_coords, svd_final

SVD, svd_model = get_optimal_svd_coords_v2(ct)

# Extract residue feature loadings from the fitted SVD model. The loading is
# scaled by the singular value so it reflects each residue's contribution in
# the retained coordinate space, rather than only its unit component weight.
if svd_model is not None:
    loading_matrix = svd_model.components_.T * svd_model.singular_values_
    loading_columns = [f'component_{i + 1}_loading' for i in range(loading_matrix.shape[1])]
    loading_df = pd.DataFrame(loading_matrix, columns=loading_columns)
    loading_df.insert(0, 'residue', residue_features)
    loading_df['max_absolute_loading'] = loading_df[loading_columns].abs().max(axis=1)
    loading_df['dominant_component'] = (
        loading_df[loading_columns].abs().idxmax(axis=1).str.replace('_loading', '', regex=False)
    )
    loading_df = loading_df.sort_values('max_absolute_loading', ascending=False)

    loading_df.to_csv(
        '/Users/Zen/Downloads/clustaldm_results_20260925_054426/svd_residue_feature_loadings.csv',
        index=False
    )
    print('\nTop predicted interacting residues by SVD loading:')
    for component_number, column in enumerate(loading_columns, start=1):
        top_residues = loading_df.loc[
            loading_df[column].abs().nlargest(10).index,
            ['residue', column]
        ]
        print(f'  Component {component_number}:')
        for _, row in top_residues.iterrows():
            print(f"    {row['residue']}: {row[column]:.6f}")
    print('Saved SVD residue loadings to svd_residue_feature_loadings.csv')
else:
    loading_df = pd.DataFrame()
    print('SVD residue loadings skipped: fewer than 3 samples/features available.')

# STEP 3: Clustering using Leiden algorithm
# 3.1. BUILD THE KNN GRAPH FROM SVD COORDINATES
# We use the clean, noise-free continuous dimensions from your geometric elbow
k_neighbors = 20
knn_graph = kneighbors_graph(SVD, n_neighbors=k_neighbors, mode='connectivity', include_self=True)

# 3.2. CONVERT KNN TO AN SNN WEIGHTED MATRIX (Jaccard-like neighborhood overlap)
snn_matrix = knn_graph.dot(knn_graph.T)

# The SNN matrix contains diagonal self-links. Remove them before testing for
# isolated points so a model counts as connected only if it shares neighbors
# with at least one other model.
snn_network = snn_matrix.copy().tolil()
snn_network.setdiag(0)
snn_network = snn_network.tocsr()
network_degree = np.asarray((snn_network > 0).sum(axis=1)).ravel()
isolated_indices = np.flatnonzero(network_degree == 0)

if len(isolated_indices) > 0:
    isolated_models = df.iloc[isolated_indices].copy()
    isolated_models.insert(0, 'model_index', isolated_indices)
    isolated_models['network_degree'] = network_degree[isolated_indices]
    isolated_models.to_csv(
        '/Users/Zen/Downloads/clustaldm_results_20260925_054426/isolated_network_points.csv',
        index=False
    )
    print(f'\nIsolated models unable to join an SNN network: {len(isolated_indices)}')
    print(isolated_models[['model_index', 'structure', 'network_degree']].to_string(index=False))
    print('Saved isolated model data to isolated_network_points.csv')
else:
    print('\nNo isolated models found in the off-diagonal SNN network.')

# 3.3. CONVERT SNN MATRIX TO AN IGRAPH OBJECT FOR LEIDEN
# We use the shared neighbor counts as edge weights
sources, targets = snn_network.nonzero()
weights = snn_network.data
g = ig.Graph(directed=False)
g.add_vertices(snn_matrix.shape[0])
g.add_edges(list(zip(sources, targets)), attributes={'weight': weights})

# 3.4. RUN LEIDEN CLUSTERING
# Adjust 'resolution_parameter' to get more (higher value) or fewer (lower value) clusters
partition = la.find_partition(g, la.ModularityVertexPartition, weights='weight')
cluster_labels = np.array(partition.membership)

print(f"SNN-Leiden pipeline successfully detected {len(np.unique(cluster_labels))} clusters!")

# Export one row per model so cluster assignments can be inspected alongside
# the contact-point data used to build the feature matrix.
if len(cluster_labels) != len(df):
    raise ValueError(
        f"Cluster label count ({len(cluster_labels)}) does not match model count ({len(df)})."
    )

model_cluster_stats = df.copy()
model_cluster_stats.insert(0, 'model_index', np.arange(len(model_cluster_stats)))
model_cluster_stats.insert(1, 'cluster', cluster_labels)

contact_columns = [
    column for column in ('protein_residues', 'nucleic_residues', 'contact_points')
    if column in model_cluster_stats.columns
]
preferred_columns = ['model_index', 'structure', 'cluster'] + contact_columns
remaining_columns = [
    column for column in model_cluster_stats.columns
    if column not in preferred_columns
]
model_cluster_stats = model_cluster_stats[
    [column for column in preferred_columns if column in model_cluster_stats.columns]
    + remaining_columns
].sort_values(['cluster', 'model_index'])

model_cluster_stats.to_csv(
    '/Users/Zen/Downloads/clustaldm_results_20260925_054426/model_cluster_contact_statistics.csv',
    index=False
)
print(
    "Saved per-model cluster assignments and contact points to "
    "model_cluster_contact_statistics.csv"
)


#UMAP dimensionality reduction visualisation and validation
umap_reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
umap_embedding = umap_reducer.fit_transform(SVD)

plt.figure(figsize=(9, 7))
scatter = plt.scatter(
    umap_embedding[:, 0], 
    umap_embedding[:, 1], 
    c=cluster_labels, 
    cmap='tab10',  # Distinct qualitative color palette for clusters
    s=12, 
    alpha=0.8
)

# Add a clean legend showing the cluster IDs
legend = plt.legend(*scatter.legend_elements(), title="SNN Clusters", loc="upper right")
plt.gca().add_artist(legend)

plt.title('UMAP Layout with Shared Nearest Neighbor (SNN) Clustering')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.grid(True, linestyle='--', alpha=0.3)
plt.savefig('/Users/Zen/Downloads/clustaldm_results_20260925_054426/umap_snn_leiden_clustering.png', dpi=300, bbox_inches='tight')
plt.show()