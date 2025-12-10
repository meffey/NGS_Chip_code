#!/usr/bin/env python3
"""
Single-cell RNA-seq Dimensionality Reduction Script
单细胞RNA-seq降维分析脚本

This script performs dimensionality reduction including:
- Principal Component Analysis (PCA)
- t-SNE
- UMAP

Usage: python scrna_reduce_dims.py <input_h5ad> <output_h5ad> [options]
"""

import argparse
import sys
import scanpy as sc
import matplotlib.pyplot as plt

def setup_scanpy():
    """Configure scanpy settings"""
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white')

def run_pca(adata, n_comps=50, use_hvg=True, plot=True, output_prefix='pca'):
    """
    Run Principal Component Analysis
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    n_comps : int
        Number of principal components
    use_hvg : bool
        Use only highly variable genes
    plot : bool
        Generate PCA plots
    output_prefix : str
        Prefix for output plots
    
    Returns:
    --------
    adata : AnnData
        Data with PCA results
    """
    print(f"\nRunning PCA (n_comps={n_comps}, use_hvg={use_hvg})...")
    
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps, use_highly_variable=use_hvg)
    
    if plot:
        # Variance ratio plot
        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save=f'_{output_prefix}_variance.png')
        
        # PCA scatter plots
        sc.pl.pca(adata, color=['total_counts', 'n_genes_by_counts'], 
                 save=f'_{output_prefix}_scatter.png')
    
    print(f"PCA completed. Top PC explains {adata.uns['pca']['variance_ratio'][0]:.2%} variance")
    
    return adata

def compute_neighbors(adata, n_neighbors=10, n_pcs=40, metric='euclidean'):
    """
    Compute neighborhood graph
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    n_neighbors : int
        Number of neighbors
    n_pcs : int
        Number of PCs to use
    metric : str
        Distance metric
    
    Returns:
    --------
    adata : AnnData
        Data with neighbor graph
    """
    print(f"\nComputing neighbor graph (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")
    
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric=metric)
    
    return adata

def run_tsne(adata, n_pcs=40, plot=True, output_prefix='tsne'):
    """
    Run t-SNE
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    n_pcs : int
        Number of PCs to use
    plot : bool
        Generate t-SNE plots
    output_prefix : str
        Prefix for output plots
    
    Returns:
    --------
    adata : AnnData
        Data with t-SNE results
    """
    print(f"\nRunning t-SNE (n_pcs={n_pcs})...")
    
    sc.tl.tsne(adata, n_pcs=n_pcs)
    
    if plot:
        sc.pl.tsne(adata, color=['total_counts', 'n_genes_by_counts'], 
                  save=f'_{output_prefix}.png')
    
    print("t-SNE completed")
    
    return adata

def run_umap(adata, min_dist=0.5, spread=1.0, plot=True, output_prefix='umap'):
    """
    Run UMAP
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    min_dist : float
        Minimum distance parameter
    spread : float
        Spread parameter
    plot : bool
        Generate UMAP plots
    output_prefix : str
        Prefix for output plots
    
    Returns:
    --------
    adata : AnnData
        Data with UMAP results
    """
    print(f"\nRunning UMAP (min_dist={min_dist}, spread={spread})...")
    
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)
    
    if plot:
        sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts'], 
                  save=f'_{output_prefix}.png')
    
    print("UMAP completed")
    
    return adata

def run_diffmap(adata, n_comps=15):
    """
    Run Diffusion Map
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    n_comps : int
        Number of diffusion components
    
    Returns:
    --------
    adata : AnnData
        Data with diffusion map results
    """
    print(f"\nRunning Diffusion Map (n_comps={n_comps})...")
    
    sc.tl.diffmap(adata, n_comps=n_comps)
    
    print("Diffusion Map completed")
    
    return adata

def main():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Dimensionality Reduction')
    parser.add_argument('input', help='Input h5ad file (preprocessed data)')
    parser.add_argument('output', help='Output h5ad file')
    parser.add_argument('--n-pcs', type=int, default=50,
                       help='Number of principal components (default: 50)')
    parser.add_argument('--use-hvg', action='store_true', default=True,
                       help='Use only highly variable genes for PCA')
    parser.add_argument('--n-neighbors', type=int, default=10,
                       help='Number of neighbors for graph (default: 10)')
    parser.add_argument('--n-pcs-use', type=int, default=40,
                       help='Number of PCs to use for neighbor graph and t-SNE (default: 40)')
    parser.add_argument('--run-tsne', action='store_true',
                       help='Run t-SNE (may be slow for large datasets)')
    parser.add_argument('--run-umap', action='store_true', default=True,
                       help='Run UMAP (default: True)')
    parser.add_argument('--run-diffmap', action='store_true',
                       help='Run Diffusion Map')
    parser.add_argument('--umap-min-dist', type=float, default=0.5,
                       help='UMAP min_dist parameter (default: 0.5)')
    parser.add_argument('--umap-spread', type=float, default=1.0,
                       help='UMAP spread parameter (default: 1.0)')
    parser.add_argument('--plot-prefix', default='reduction',
                       help='Prefix for output plots (default: reduction)')
    
    args = parser.parse_args()
    
    # Setup
    setup_scanpy()
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Run PCA
    adata = run_pca(adata, n_comps=args.n_pcs, use_hvg=args.use_hvg, 
                   plot=True, output_prefix=args.plot_prefix)
    
    # Compute neighbors
    adata = compute_neighbors(adata, n_neighbors=args.n_neighbors, 
                             n_pcs=args.n_pcs_use)
    
    # Run t-SNE (optional)
    if args.run_tsne:
        adata = run_tsne(adata, n_pcs=args.n_pcs_use, plot=True, 
                        output_prefix=args.plot_prefix)
    
    # Run UMAP
    if args.run_umap:
        adata = run_umap(adata, min_dist=args.umap_min_dist, 
                        spread=args.umap_spread, plot=True, 
                        output_prefix=args.plot_prefix)
    
    # Run Diffusion Map (optional)
    if args.run_diffmap:
        adata = run_diffmap(adata)
    
    # Save results
    print(f"\nSaving results to {args.output}...")
    adata.write(args.output)
    
    print("\nDimensionality reduction completed successfully!")
    print(f"Dataset: {adata.n_obs} cells × {adata.n_vars} genes")

if __name__ == '__main__':
    main()
