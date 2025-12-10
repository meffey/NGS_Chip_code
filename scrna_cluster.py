#!/usr/bin/env python3
"""
Single-cell RNA-seq Clustering Script
单细胞RNA-seq细胞聚类脚本

This script performs cell clustering using various algorithms:
- Louvain clustering
- Leiden clustering
- K-means clustering

Usage: python scrna_cluster.py <input_h5ad> <output_h5ad> [options]
"""

import argparse
import sys
import scanpy as sc
import matplotlib.pyplot as plt

def setup_scanpy():
    """Configure scanpy settings"""
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white')

def run_louvain(adata, resolution=1.0, plot=True, output_prefix='louvain'):
    """
    Run Louvain clustering
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    resolution : float
        Resolution parameter
    plot : bool
        Generate clustering plots
    output_prefix : str
        Prefix for output plots
    
    Returns:
    --------
    adata : AnnData
        Data with Louvain clustering results
    """
    print(f"\nRunning Louvain clustering (resolution={resolution})...")
    
    sc.tl.louvain(adata, resolution=resolution)
    
    n_clusters = len(adata.obs['louvain'].unique())
    print(f"Found {n_clusters} clusters")
    
    if plot:
        # Plot on UMAP if available
        if 'X_umap' in adata.obsm:
            sc.pl.umap(adata, color='louvain', save=f'_{output_prefix}_umap.png')
        
        # Plot on t-SNE if available
        if 'X_tsne' in adata.obsm:
            sc.pl.tsne(adata, color='louvain', save=f'_{output_prefix}_tsne.png')
    
    return adata

def run_leiden(adata, resolution=1.0, plot=True, output_prefix='leiden'):
    """
    Run Leiden clustering
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    resolution : float
        Resolution parameter
    plot : bool
        Generate clustering plots
    output_prefix : str
        Prefix for output plots
    
    Returns:
    --------
    adata : AnnData
        Data with Leiden clustering results
    """
    print(f"\nRunning Leiden clustering (resolution={resolution})...")
    
    sc.tl.leiden(adata, resolution=resolution)
    
    n_clusters = len(adata.obs['leiden'].unique())
    print(f"Found {n_clusters} clusters")
    
    if plot:
        # Plot on UMAP if available
        if 'X_umap' in adata.obsm:
            sc.pl.umap(adata, color='leiden', save=f'_{output_prefix}_umap.png')
        
        # Plot on t-SNE if available
        if 'X_tsne' in adata.obsm:
            sc.pl.tsne(adata, color='leiden', save=f'_{output_prefix}_tsne.png')
    
    return adata

def find_marker_genes(adata, groupby='leiden', method='wilcoxon', n_genes=25):
    """
    Find marker genes for each cluster
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    groupby : str
        Cluster annotation to use
    method : str
        Statistical test method
    n_genes : int
        Number of top genes per cluster
    
    Returns:
    --------
    adata : AnnData
        Data with marker gene results
    """
    print(f"\nFinding marker genes (method={method}, n_genes={n_genes})...")
    
    sc.tl.rank_genes_groups(adata, groupby, method=method, n_genes=n_genes)
    
    return adata

def plot_marker_genes(adata, n_genes=25, output_prefix='markers'):
    """
    Plot marker genes
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    n_genes : int
        Number of genes to plot
    output_prefix : str
        Prefix for output plots
    """
    print(f"\nPlotting marker genes...")
    
    # Ranking plot
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, 
                           save=f'_{output_prefix}_ranking.png')
    
    # Heatmap
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes, 
                                    save=f'_{output_prefix}_heatmap.png')
    
    # Dotplot
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, 
                                    save=f'_{output_prefix}_dotplot.png')
    
    # Violin plot for top genes
    sc.pl.rank_genes_groups_violin(adata, n_genes=min(5, n_genes), 
                                   save=f'_{output_prefix}_violin.png')

def export_marker_genes(adata, output_file, n_genes=100):
    """
    Export marker genes to file
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    output_file : str
        Output file path
    n_genes : int
        Number of genes to export
    """
    print(f"\nExporting top {n_genes} marker genes to {output_file}...")
    
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("cluster\tgene\tscore\tpval\tpval_adj\tlogfoldchange\n")
        
        # Write marker genes for each cluster
        for group in groups:
            for i in range(min(n_genes, len(result['names'][group]))):
                gene = result['names'][group][i]
                score = result['scores'][group][i]
                pval = result['pvals'][group][i]
                pval_adj = result['pvals_adj'][group][i]
                logfc = result['logfoldchanges'][group][i]
                
                f.write(f"{group}\t{gene}\t{score:.4f}\t{pval:.4e}\t{pval_adj:.4e}\t{logfc:.4f}\n")
    
    print(f"Marker genes exported successfully")

def main():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Clustering')
    parser.add_argument('input', help='Input h5ad file (with dimensionality reduction)')
    parser.add_argument('output', help='Output h5ad file')
    parser.add_argument('--method', default='leiden', choices=['louvain', 'leiden', 'both'],
                       help='Clustering method (default: leiden)')
    parser.add_argument('--resolution', type=float, default=1.0,
                       help='Clustering resolution (default: 1.0)')
    parser.add_argument('--find-markers', action='store_true', default=True,
                       help='Find marker genes for clusters')
    parser.add_argument('--marker-method', default='wilcoxon', 
                       choices=['wilcoxon', 't-test', 'logreg'],
                       help='Method for finding marker genes (default: wilcoxon)')
    parser.add_argument('--n-marker-genes', type=int, default=25,
                       help='Number of marker genes per cluster (default: 25)')
    parser.add_argument('--export-markers', type=str, default=None,
                       help='Export marker genes to file (e.g., markers.tsv)')
    parser.add_argument('--plot-prefix', default='cluster',
                       help='Prefix for output plots (default: cluster)')
    
    args = parser.parse_args()
    
    # Setup
    setup_scanpy()
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Check if neighbors have been computed
    if 'neighbors' not in adata.uns:
        print("\nWarning: Neighbor graph not found. Computing neighbors...")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Run clustering
    if args.method in ['louvain', 'both']:
        adata = run_louvain(adata, resolution=args.resolution, 
                           plot=True, output_prefix=args.plot_prefix)
    
    if args.method in ['leiden', 'both']:
        adata = run_leiden(adata, resolution=args.resolution, 
                          plot=True, output_prefix=args.plot_prefix)
    
    # Find marker genes
    if args.find_markers:
        # Use leiden if available, otherwise louvain
        groupby = 'leiden' if 'leiden' in adata.obs else 'louvain'
        
        adata = find_marker_genes(adata, groupby=groupby, 
                                 method=args.marker_method, 
                                 n_genes=args.n_marker_genes)
        
        # Plot marker genes
        plot_marker_genes(adata, n_genes=args.n_marker_genes, 
                         output_prefix=args.plot_prefix)
        
        # Export marker genes if requested
        if args.export_markers:
            export_marker_genes(adata, args.export_markers, 
                              n_genes=100)
    
    # Save results
    print(f"\nSaving results to {args.output}...")
    adata.write(args.output)
    
    print("\nClustering completed successfully!")
    print(f"Dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    if 'leiden' in adata.obs:
        print(f"Leiden clusters: {len(adata.obs['leiden'].unique())}")
    if 'louvain' in adata.obs:
        print(f"Louvain clusters: {len(adata.obs['louvain'].unique())}")

if __name__ == '__main__':
    main()
