#!/usr/bin/env python3
"""
Single-cell RNA-seq Differential Expression Analysis Script
单细胞RNA-seq差异表达基因分析脚本

This script performs differential expression analysis between groups of cells

Usage: python scrna_diff_expr.py <input_h5ad> <output_prefix> [options]
"""

import argparse
import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

def setup_scanpy():
    """Configure scanpy settings"""
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white')

def run_de_analysis(adata, groupby, group1, group2=None, method='wilcoxon'):
    """
    Run differential expression analysis
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    groupby : str
        Column in adata.obs to group by
    group1 : str
        First group for comparison
    group2 : str or None
        Second group for comparison (None for one-vs-rest)
    method : str
        Statistical test method
    
    Returns:
    --------
    result : DataFrame
        Differential expression results
    """
    if group2 is None:
        print(f"\nRunning DE analysis: {group1} vs rest (method={method})...")
        sc.tl.rank_genes_groups(adata, groupby, groups=[group1], 
                               reference='rest', method=method)
    else:
        print(f"\nRunning DE analysis: {group1} vs {group2} (method={method})...")
        sc.tl.rank_genes_groups(adata, groupby, groups=[group1], 
                               reference=group2, method=method)
    
    # Extract results
    result = sc.get.rank_genes_groups_df(adata, group=group1)
    
    return result

def filter_de_genes(result, pval_cutoff=0.05, logfc_cutoff=1.0):
    """
    Filter differentially expressed genes
    
    Parameters:
    -----------
    result : DataFrame
        DE analysis results
    pval_cutoff : float
        P-value cutoff
    logfc_cutoff : float
        Log fold change cutoff
    
    Returns:
    --------
    filtered : DataFrame
        Filtered DE genes
    """
    print(f"\nFiltering DE genes (pval < {pval_cutoff}, |logFC| > {logfc_cutoff})...")
    
    filtered = result[
        (result['pvals_adj'] < pval_cutoff) & 
        (abs(result['logfoldchanges']) > logfc_cutoff)
    ]
    
    print(f"Found {len(filtered)} significant DE genes")
    
    return filtered

def plot_de_results(adata, result, top_n=20, output_prefix='de'):
    """
    Plot differential expression results
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    result : DataFrame
        DE analysis results
    top_n : int
        Number of top genes to plot
    output_prefix : str
        Prefix for output plots
    """
    print(f"\nGenerating DE plots...")
    
    # Get top genes
    top_genes = result.nsmallest(top_n, 'pvals_adj')['names'].tolist()
    
    # Volcano plot-like scatter
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(result['logfoldchanges'], -np.log10(result['pvals_adj']), 
              alpha=0.5, s=10)
    ax.axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    ax.axvline(-1, color='blue', linestyle='--', alpha=0.5)
    ax.axvline(1, color='blue', linestyle='--', alpha=0.5, label='|logFC|=1')
    ax.set_xlabel('Log Fold Change')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title('Volcano Plot')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'figures/{output_prefix}_volcano.png', dpi=150)
    plt.close()
    
    # Heatmap of top genes
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=top_genes[:10], 
                  save=f'_{output_prefix}_top_genes_umap.png')
    
    print(f"DE plots saved with prefix: {output_prefix}")

def export_de_results(result, output_file):
    """
    Export DE results to file
    
    Parameters:
    -----------
    result : DataFrame
        DE analysis results
    output_file : str
        Output file path
    """
    print(f"\nExporting DE results to {output_file}...")
    
    result.to_csv(output_file, sep='\t', index=False)
    
    print(f"DE results exported successfully")

def main():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Differential Expression Analysis')
    parser.add_argument('input', help='Input h5ad file (clustered data)')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('--groupby', default='leiden',
                       help='Column in adata.obs to group by (default: leiden)')
    parser.add_argument('--group1', required=True,
                       help='First group for comparison (required)')
    parser.add_argument('--group2', default=None,
                       help='Second group for comparison (default: None for one-vs-rest)')
    parser.add_argument('--method', default='wilcoxon', 
                       choices=['wilcoxon', 't-test', 'logreg', 't-test_overestim_var'],
                       help='Statistical test method (default: wilcoxon)')
    parser.add_argument('--pval-cutoff', type=float, default=0.05,
                       help='Adjusted p-value cutoff (default: 0.05)')
    parser.add_argument('--logfc-cutoff', type=float, default=1.0,
                       help='Log fold change cutoff (default: 1.0)')
    parser.add_argument('--top-n', type=int, default=20,
                       help='Number of top genes to plot (default: 20)')
    
    args = parser.parse_args()
    
    # Setup
    setup_scanpy()
    
    # Import numpy for plotting
    import numpy as np
    globals()['np'] = np
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Check if groupby column exists
    if args.groupby not in adata.obs.columns:
        print(f"\nError: Column '{args.groupby}' not found in data")
        print(f"Available columns: {', '.join(adata.obs.columns)}")
        sys.exit(1)
    
    # Check if groups exist
    available_groups = adata.obs[args.groupby].unique().tolist()
    if args.group1 not in available_groups:
        print(f"\nError: Group '{args.group1}' not found")
        print(f"Available groups: {', '.join(map(str, available_groups))}")
        sys.exit(1)
    
    if args.group2 is not None and args.group2 not in available_groups:
        print(f"\nError: Group '{args.group2}' not found")
        print(f"Available groups: {', '.join(map(str, available_groups))}")
        sys.exit(1)
    
    # Run DE analysis
    result = run_de_analysis(adata, args.groupby, args.group1, args.group2, args.method)
    
    # Export all results
    all_results_file = f"{args.output_prefix}_all_genes.tsv"
    export_de_results(result, all_results_file)
    
    # Filter significant genes
    filtered_result = filter_de_genes(result, args.pval_cutoff, args.logfc_cutoff)
    
    # Export filtered results
    if len(filtered_result) > 0:
        filtered_file = f"{args.output_prefix}_significant_genes.tsv"
        export_de_results(filtered_result, filtered_file)
    
    # Plot results
    plot_de_results(adata, result, args.top_n, args.output_prefix)
    
    print("\nDifferential expression analysis completed successfully!")
    print(f"Total genes tested: {len(result)}")
    print(f"Significant genes: {len(filtered_result)}")

if __name__ == '__main__':
    main()
