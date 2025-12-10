#!/usr/bin/env python3
"""
Single-cell RNA-seq Visualization Script
单细胞RNA-seq数据可视化脚本

This script generates various visualizations for single-cell RNA-seq data

Usage: python scrna_visualize.py <input_h5ad> <output_dir> [options]
"""

import argparse
import os
import sys
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def setup_scanpy(output_dir):
    """Configure scanpy and matplotlib settings"""
    sc.settings.verbosity = 2
    
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    sc.settings.figdir = figures_dir
    sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)
    
    return figures_dir

def plot_clustering(adata, groupby='leiden', output_prefix='cluster'):
    """Plot clustering results"""
    print(f"\nGenerating clustering plots...")
    
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=groupby, legend_loc='on data', 
                  save=f'_{output_prefix}_umap.png')
        sc.pl.umap(adata, color=groupby, legend_loc='right margin',
                  save=f'_{output_prefix}_umap_legend.png')
    
    if 'X_tsne' in adata.obsm:
        sc.pl.tsne(adata, color=groupby, legend_loc='on data',
                  save=f'_{output_prefix}_tsne.png')

def plot_qc_metrics(adata, groupby='leiden', output_prefix='qc'):
    """Plot QC metrics by cluster"""
    print(f"\nGenerating QC metric plots...")
    
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                  save=f'_{output_prefix}_metrics_umap.png')
    
    # Violin plots
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                groupby=groupby, save=f'_{output_prefix}_metrics_violin.png')

def plot_gene_expression(adata, genes, output_prefix='genes'):
    """Plot gene expression patterns"""
    print(f"\nGenerating gene expression plots...")
    
    # Filter genes that exist in the data
    available_genes = [g for g in genes if g in adata.var_names]
    
    if not available_genes:
        print(f"Warning: None of the specified genes found in data")
        return
    
    print(f"Plotting {len(available_genes)} genes: {', '.join(available_genes)}")
    
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=available_genes, save=f'_{output_prefix}_umap.png')
    
    # Dot plot
    if 'leiden' in adata.obs.columns:
        sc.pl.dotplot(adata, available_genes, groupby='leiden', 
                     save=f'_{output_prefix}_dotplot.png')
    
    # Violin plot
    if 'leiden' in adata.obs.columns:
        sc.pl.violin(adata, available_genes, groupby='leiden',
                    save=f'_{output_prefix}_violin.png')

def plot_top_markers(adata, n_genes=10, groupby='leiden', output_prefix='markers'):
    """Plot top marker genes"""
    print(f"\nGenerating marker gene plots...")
    
    if 'rank_genes_groups' not in adata.uns:
        print("Warning: Marker genes not found. Run clustering with marker detection first.")
        return
    
    # Ranking plot
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False,
                           save=f'_{output_prefix}_ranking.png')
    
    # Heatmap
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes, groupby=groupby,
                                    save=f'_{output_prefix}_heatmap.png')
    
    # Dotplot
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, groupby=groupby,
                                    save=f'_{output_prefix}_dotplot.png')
    
    # Stacked violin plot
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=n_genes, groupby=groupby,
                                          save=f'_{output_prefix}_stacked_violin.png')
    
    # Matrix plot
    sc.pl.rank_genes_groups_matrixplot(adata, n_genes=n_genes, groupby=groupby,
                                      save=f'_{output_prefix}_matrix.png')

def plot_trajectory(adata, root_cluster='0', output_prefix='trajectory'):
    """Plot pseudotime trajectory (if applicable)"""
    print(f"\nGenerating trajectory plots...")
    
    if 'leiden' not in adata.obs.columns:
        print("Warning: Clustering information not found")
        return
    
    # Set root cell
    root_cells = adata.obs['leiden'] == root_cluster
    if not any(root_cells):
        print(f"Warning: Root cluster '{root_cluster}' not found")
        return
    
    try:
        adata.uns['iroot'] = root_cells.idxmax()
    except ValueError:
        print(f"Warning: Could not set root cell for cluster '{root_cluster}'")
        return
    
    # Run diffusion pseudotime
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    
    # Plot
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=['dpt_pseudotime'], save=f'_{output_prefix}_dpt_umap.png')

def create_summary_plot(adata, output_file='summary.png'):
    """Create a summary figure with multiple panels"""
    print(f"\nCreating summary plot...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    
    # UMAP with clusters
    if 'X_umap' in adata.obsm and 'leiden' in adata.obs.columns:
        sc.pl.umap(adata, color='leiden', ax=axes[0, 0], show=False, legend_loc='on data')
        axes[0, 0].set_title('Cell Clusters')
    
    # UMAP with gene counts
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[0, 1], show=False)
        axes[0, 1].set_title('Gene Counts per Cell')
    
    # QC violin plot
    if 'leiden' in adata.obs.columns:
        sc.pl.violin(adata, 'n_genes_by_counts', groupby='leiden', 
                    ax=axes[1, 0], show=False)
        axes[1, 0].set_title('Genes by Cluster')
    
    # Cluster size barplot
    if 'leiden' in adata.obs.columns:
        cluster_counts = adata.obs['leiden'].value_counts().sort_index()
        axes[1, 1].bar(range(len(cluster_counts)), cluster_counts.values)
        axes[1, 1].set_xlabel('Cluster')
        axes[1, 1].set_ylabel('Number of Cells')
        axes[1, 1].set_title('Cells per Cluster')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    
    print(f"Summary plot saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Visualization')
    parser.add_argument('input', help='Input h5ad file (processed data)')
    parser.add_argument('output_dir', help='Output directory for plots')
    parser.add_argument('--groupby', default='leiden',
                       help='Grouping variable for plots (default: leiden)')
    parser.add_argument('--genes', nargs='+', default=[],
                       help='Specific genes to plot')
    parser.add_argument('--n-marker-genes', type=int, default=10,
                       help='Number of marker genes to plot (default: 10)')
    parser.add_argument('--plot-trajectory', action='store_true',
                       help='Plot pseudotime trajectory')
    parser.add_argument('--root-cluster', default='0',
                       help='Root cluster for trajectory (default: 0)')
    parser.add_argument('--all', action='store_true',
                       help='Generate all available plots')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup
    figures_dir = setup_scanpy(args.output_dir)
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Check groupby variable
    if args.groupby not in adata.obs.columns:
        print(f"Warning: '{args.groupby}' not found. Available: {', '.join(adata.obs.columns)}")
        if 'leiden' in adata.obs.columns:
            args.groupby = 'leiden'
            print(f"Using 'leiden' instead")
        elif 'louvain' in adata.obs.columns:
            args.groupby = 'louvain'
            print(f"Using 'louvain' instead")
    
    # Generate plots
    plot_clustering(adata, args.groupby)
    plot_qc_metrics(adata, args.groupby)
    
    if args.genes:
        plot_gene_expression(adata, args.genes)
    
    if args.all or 'rank_genes_groups' in adata.uns:
        plot_top_markers(adata, args.n_marker_genes, args.groupby)
    
    if args.plot_trajectory:
        plot_trajectory(adata, args.root_cluster)
    
    # Create summary plot
    summary_file = os.path.join(args.output_dir, 'summary.png')
    create_summary_plot(adata, summary_file)
    
    print("\n" + "="*60)
    print("Visualization completed successfully!")
    print("="*60)
    print(f"Output directory: {args.output_dir}")
    print(f"Figures directory: {figures_dir}")

if __name__ == '__main__':
    main()
