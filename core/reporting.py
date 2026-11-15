"""Shared report writers used by the CLI and Shiny exports."""

from pathlib import Path


def _save_plotly_html(figure, path):
    path.write_text(
        figure.to_html(
            full_html=True,
            include_plotlyjs='cdn',
            config={'displayModeBar': True, 'displaylogo': False}
        ),
        encoding='utf-8'
    )


def write_svd_diagnostics(clusterer, output_dir):
    """Write the SVD/Leiden CSV and HTML diagnostics for one analysis."""
    import plotly.graph_objects as go

    output_dir = Path(output_dir)
    csv_dir = output_dir / 'csv_exports'
    html_dir = output_dir / 'html_reports'
    csv_dir.mkdir(parents=True, exist_ok=True)
    html_dir.mkdir(parents=True, exist_ok=True)
    written = []

    elbow_df = getattr(clusterer, 'svd_elbow_data', None)
    if elbow_df is not None and not elbow_df.empty:
        elbow_df.to_csv(csv_dir / 'svd_elbow_curve.csv', index=False)
        figure = go.Figure(go.Scatter(
            x=elbow_df['component'],
            y=elbow_df['explained_variance_ratio'],
            mode='lines+markers',
            name='Explained variance ratio'
        ))
        selected = elbow_df[elbow_df['selected']]
        if not selected.empty:
            figure.add_vline(
                x=int(selected.iloc[0]['component']),
                line_dash='dash',
                line_color='red'
            )
        figure.update_layout(
            title='SVD Component Elbow',
            xaxis_title='SVD components',
            yaxis_title='Explained variance ratio',
            height=450
        )
        _save_plotly_html(figure, html_dir / 'svd_elbow.html')
        written.extend([
            csv_dir / 'svd_elbow_curve.csv',
            html_dir / 'svd_elbow.html'
        ])

    cluster_elbow_df = getattr(clusterer, 'svd_cluster_elbow_data', None)
    if cluster_elbow_df is not None and not cluster_elbow_df.empty:
        cluster_elbow_df.to_csv(
            csv_dir / 'svd_cluster_count_wcss_curve.csv', index=False
        )
        y_column = (
            'wcss_per_point' if 'wcss_per_point' in cluster_elbow_df.columns
            else 'wcss'
        )
        figure = go.Figure(go.Scatter(
            x=cluster_elbow_df['n_clusters'],
            y=cluster_elbow_df[y_column],
            customdata=cluster_elbow_df['wcss'],
            mode='lines+markers',
            name='KMeans WCSS per datapoint',
            hovertemplate=(
                'Clusters: %{x}<br>WCSS per datapoint: %{y:.4f}'
                '<br>Raw WCSS: %{customdata:.2f}<extra></extra>'
            )
        ))
        selected = cluster_elbow_df[cluster_elbow_df['selected']]
        if not selected.empty:
            figure.add_vline(
                x=int(selected.iloc[0]['n_clusters']),
                line_dash='dash',
                line_color='red'
            )
        figure.update_layout(
            title='Cluster-Count Elbow: WCSS per Datapoint',
            xaxis_title='Number of clusters',
            yaxis_title='WCSS per datapoint',
            height=450
        )
        _save_plotly_html(figure, html_dir / 'svd_cluster_count_wcss.html')
        written.extend([
            csv_dir / 'svd_cluster_count_wcss_curve.csv',
            html_dir / 'svd_cluster_count_wcss.html'
        ])

    grid_df = getattr(clusterer, 'svd_grid_search_data', None)
    if grid_df is not None and not grid_df.empty:
        grid_df.to_csv(csv_dir / 'svd_component_k_grid_search.csv', index=False)
        score_pivot = grid_df.pivot(
            index='n_components', columns='k', values='average_intra_jaccard'
        ).sort_index()
        figure = go.Figure(go.Heatmap(
            x=score_pivot.columns.tolist(),
            y=score_pivot.index.tolist(),
            z=score_pivot.to_numpy(),
            colorscale='RdYlGn_r',
            colorbar=dict(title='Avg intra-Jaccard distance'),
            hovertemplate=(
                'Components: %{y}<br>k: %{x}<br>'
                'Avg intra-Jaccard: %{z:.4f}<extra></extra>'
            )
        ))
        selected = grid_df[grid_df['selected']]
        if not selected.empty:
            figure.add_trace(go.Scatter(
                x=[int(selected.iloc[0]['k'])],
                y=[int(selected.iloc[0]['n_components'])],
                mode='markers',
                name='Selected candidate',
                marker=dict(
                    symbol='square-open', size=16, color='black',
                    line=dict(width=2)
                )
            ))
        figure.update_layout(
            title='SVD Component / SNN k Grid Search',
            xaxis_title='SNN neighborhood size (k)',
            yaxis_title='SVD components',
            height=600
        )
        _save_plotly_html(
            figure, html_dir / 'svd_component_k_grid_search.html'
        )
        written.extend([
            csv_dir / 'svd_component_k_grid_search.csv',
            html_dir / 'svd_component_k_grid_search.html'
        ])

    return written
