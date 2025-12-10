#!/usr/bin/env python3
"""
Single-cell RNA-seq Analysis Pipeline
单细胞RNA-seq完整分析流程

This script runs a complete single-cell RNA-seq analysis pipeline including:
1. Quality control
2. Normalization and preprocessing
3. Dimensionality reduction
4. Clustering
5. Marker gene identification

Usage: python scrna_pipeline.py <input_data> <output_dir> [options]
"""

import argparse
import os
import sys
import scanpy as sc
import matplotlib.pyplot as plt

def setup_scanpy(output_dir):
    """Configure scanpy settings"""
    sc.settings.verbosity = 3
    
    # Create output directories
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    sc.settings.figdir = figures_dir
    sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)
    
    return figures_dir

def run_qc(adata, species='human', min_genes=200, max_mt_percent=20, min_cells=3):
    """Quality control step"""
    print("\n" + "="*60)
    print("STEP 1: Quality Control")
    print("="*60)
    
    # Identify mitochondrial genes
    # Note: Uses case-insensitive matching to catch MT-, Mt-, mt- variants
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Plot QC metrics before filtering
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save='_qc_before_filter.png')
    
    print(f"Before filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.pct_counts_mt < max_mt_percent, :]
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Plot QC metrics after filtering
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save='_qc_after_filter.png')
    
    return adata

def run_normalization(adata, target_sum=1e4, n_top_genes=2000):
    """Normalization and preprocessing step"""
    print("\n" + "="*60)
    print("STEP 2: Normalization and Preprocessing")
    print("="*60)
    
    # Store raw counts
    adata.raw = adata.copy()
    
    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    
    # Log transform
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print(f"Found {sum(adata.var['highly_variable'])} highly variable genes")
    
    sc.pl.highly_variable_genes(adata, save='_hvg.png')
    
    # Regress out effects
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # Scale
    sc.pp.scale(adata, max_value=10)
    
    return adata

def run_dimensionality_reduction(adata, n_pcs=50, n_neighbors=10):
    """Dimensionality reduction step"""
    print("\n" + "="*60)
    print("STEP 3: Dimensionality Reduction")
    print("="*60)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='_variance_ratio.png')
    
    # Compute neighbor graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=40)
    
    # UMAP
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'],
              save='_qc_metrics.png')
    
    return adata

def run_clustering(adata, resolution=1.0):
    """Clustering step"""
    print("\n" + "="*60)
    print("STEP 4: Clustering")
    print("="*60)
    
    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution)
    
    n_clusters = len(adata.obs['leiden'].unique())
    print(f"Found {n_clusters} clusters")
    
    # Plot clusters
    sc.pl.umap(adata, color='leiden', save='_leiden_clusters.png', legend_loc='on data')
    
    return adata

def find_markers(adata, n_genes=25):
    """Find marker genes"""
    print("\n" + "="*60)
    print("STEP 5: Marker Gene Identification")
    print("="*60)
    
    # Find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=n_genes)
    
    # Plot results
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, save='_marker_ranking.png')
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes, save='_marker_heatmap.png')
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, save='_marker_dotplot.png')
    
    return adata

def save_results(adata, output_dir):
    """Save analysis results"""
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)
    
    # Save processed data
    output_file = os.path.join(output_dir, 'processed_data.h5ad')
    adata.write(output_file)
    print(f"Processed data saved to: {output_file}")
    
    # Export marker genes
    marker_file = os.path.join(output_dir, 'marker_genes.tsv')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    with open(marker_file, 'w') as f:
        f.write("cluster\tgene\tscore\tpval\tpval_adj\tlogfoldchange\n")
        for group in groups:
            for i in range(len(result['names'][group])):
                gene = result['names'][group][i]
                score = result['scores'][group][i]
                pval = result['pvals'][group][i]
                pval_adj = result['pvals_adj'][group][i]
                logfc = result['logfoldchanges'][group][i]
                f.write(f"{group}\t{gene}\t{score:.4f}\t{pval:.4e}\t{pval_adj:.4e}\t{logfc:.4f}\n")
    
    print(f"Marker genes saved to: {marker_file}")
    
    # Export cluster information
    cluster_file = os.path.join(output_dir, 'cluster_assignments.tsv')
    adata.obs[['leiden', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt']].to_csv(
        cluster_file, sep='\t')
    print(f"Cluster assignments saved to: {cluster_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Single-cell RNA-seq Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Analyze 10x data
  python scrna_pipeline.py /path/to/10x_data output_dir --format 10x
  
  # Analyze h5ad file
  python scrna_pipeline.py data.h5ad output_dir --format h5ad
  
  # Customize parameters
  python scrna_pipeline.py data.h5ad output_dir --species mouse --min-genes 300 --resolution 0.5
        """
    )
    
    parser.add_argument('input', help='Input data (h5ad file, csv file, or 10x directory)')
    parser.add_argument('output_dir', help='Output directory for results')
    parser.add_argument('--format', default='auto', choices=['auto', 'h5ad', 'csv', '10x'],
                       help='Input data format (default: auto)')
    parser.add_argument('--species', default='human', choices=['human', 'mouse', 'other'],
                       help='Species (default: human)')
    parser.add_argument('--min-genes', type=int, default=200,
                       help='Minimum genes per cell (default: 200)')
    parser.add_argument('--max-mt-percent', type=float, default=20,
                       help='Maximum mitochondrial percentage (default: 20)')
    parser.add_argument('--min-cells', type=int, default=3,
                       help='Minimum cells per gene (default: 3)')
    parser.add_argument('--n-top-genes', type=int, default=2000,
                       help='Number of highly variable genes (default: 2000)')
    parser.add_argument('--n-pcs', type=int, default=50,
                       help='Number of principal components (default: 50)')
    parser.add_argument('--n-neighbors', type=int, default=10,
                       help='Number of neighbors (default: 10)')
    parser.add_argument('--resolution', type=float, default=1.0,
                       help='Clustering resolution (default: 1.0)')
    parser.add_argument('--n-marker-genes', type=int, default=25,
                       help='Number of marker genes per cluster (default: 25)')
    
    args = parser.parse_args()
    
    # Validate input parameters
    if args.min_genes <= 0:
        print("Error: --min-genes must be greater than 0")
        sys.exit(1)
    if args.max_mt_percent <= 0 or args.max_mt_percent >= 100:
        print("Error: --max-mt-percent must be between 0 and 100")
        sys.exit(1)
    if args.resolution <= 0:
        print("Error: --resolution must be greater than 0")
        sys.exit(1)
    if args.n_top_genes <= 0:
        print("Error: --n-top-genes must be greater than 0")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup
    figures_dir = setup_scanpy(args.output_dir)
    
    # Load data
    print("="*60)
    print("Loading Data")
    print("="*60)
    print(f"Input: {args.input}")
    print(f"Format: {args.format}")
    
    if args.format == 'auto':
        if args.input.endswith('.h5ad'):
            adata = sc.read_h5ad(args.input)
        elif args.input.endswith('.csv'):
            adata = sc.read_csv(args.input)
        elif os.path.isdir(args.input):
            adata = sc.read_10x_mtx(args.input, var_names='gene_symbols', cache=True)
        else:
            print("Error: Cannot determine input format. Please specify --format")
            sys.exit(1)
    elif args.format == 'h5ad':
        adata = sc.read_h5ad(args.input)
    elif args.format == 'csv':
        adata = sc.read_csv(args.input)
    elif args.format == '10x':
        adata = sc.read_10x_mtx(args.input, var_names='gene_symbols', cache=True)
    
    print(f"Loaded dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Run pipeline
    adata = run_qc(adata, args.species, args.min_genes, args.max_mt_percent, args.min_cells)
    adata = run_normalization(adata, n_top_genes=args.n_top_genes)
    adata = run_dimensionality_reduction(adata, args.n_pcs, args.n_neighbors)
    adata = run_clustering(adata, args.resolution)
    adata = find_markers(adata, args.n_marker_genes)
    
    # Save results
    save_results(adata, args.output_dir)
    
    print("\n" + "="*60)
    print("Pipeline Completed Successfully!")
    print("="*60)
    print(f"Output directory: {args.output_dir}")
    print(f"Final dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Number of clusters: {len(adata.obs['leiden'].unique())}")
    print(f"Figures saved in: {figures_dir}")

if __name__ == '__main__':
    main()
