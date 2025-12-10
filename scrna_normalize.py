#!/usr/bin/env python3
"""
Single-cell RNA-seq Normalization and Preprocessing Script
单细胞RNA-seq数据标准化和预处理脚本

This script performs normalization and preprocessing including:
- Normalization of expression data
- Log transformation
- Identification of highly variable genes
- Scaling and centering

Usage: python scrna_normalize.py <input_h5ad> <output_h5ad> [options]
"""

import argparse
import sys
import scanpy as sc
import matplotlib.pyplot as plt

def setup_scanpy():
    """Configure scanpy settings"""
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white')

def normalize_data(adata, target_sum=1e4):
    """
    Normalize counts per cell
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    target_sum : float
        Target sum for normalization (default: 1e4)
    
    Returns:
    --------
    adata : AnnData
        Normalized data
    """
    print(f"\nNormalizing data to {target_sum} counts per cell...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    return adata

def log_transform(adata):
    """
    Log-transform the data
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    
    Returns:
    --------
    adata : AnnData
        Log-transformed data
    """
    print("\nLog-transforming data...")
    sc.pp.log1p(adata)
    return adata

def identify_hvgs(adata, n_top_genes=2000, method='seurat', plot=True, output_prefix='hvg'):
    """
    Identify highly variable genes
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    n_top_genes : int
        Number of highly variable genes to identify
    method : str
        Method for HVG identification ('seurat', 'cell_ranger', 'seurat_v3')
    plot : bool
        Whether to generate HVG plot
    output_prefix : str
        Prefix for output plot
    
    Returns:
    --------
    adata : AnnData
        Data with HVG annotations
    """
    print(f"\nIdentifying {n_top_genes} highly variable genes using {method} method...")
    
    # Identify highly variable genes
    # Parameters explanation:
    # - min_mean: minimum average expression level (0.0125)
    # - max_mean: maximum average expression level (3)
    # - min_disp: minimum dispersion (variability) (0.5)
    # These thresholds help identify genes with biologically meaningful variation
    if method == 'seurat':
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    elif method == 'cell_ranger':
        sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=n_top_genes)
    elif method == 'seurat_v3':
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=n_top_genes)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    print(f"Found {sum(adata.var['highly_variable'])} highly variable genes")
    
    if plot:
        sc.pl.highly_variable_genes(adata, save=f'_{output_prefix}.png')
    
    return adata

def regress_out_effects(adata, vars_to_regress=['total_counts', 'pct_counts_mt']):
    """
    Regress out unwanted sources of variation
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    vars_to_regress : list
        Variables to regress out
    
    Returns:
    --------
    adata : AnnData
        Data with regressed effects
    """
    print(f"\nRegressing out effects of: {', '.join(vars_to_regress)}...")
    
    # Check which variables exist
    available_vars = [var for var in vars_to_regress if var in adata.obs.columns]
    
    if available_vars:
        sc.pp.regress_out(adata, available_vars)
        print(f"Regressed out: {', '.join(available_vars)}")
    else:
        print("Warning: None of the specified variables found in data")
    
    return adata

def scale_data(adata, max_value=10):
    """
    Scale data to unit variance and zero mean
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    max_value : float
        Clip values to this maximum
    
    Returns:
    --------
    adata : AnnData
        Scaled data
    """
    print(f"\nScaling data (max value: {max_value})...")
    sc.pp.scale(adata, max_value=max_value)
    return adata

def main():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Normalization and Preprocessing')
    parser.add_argument('input', help='Input h5ad file (quality-controlled data)')
    parser.add_argument('output', help='Output h5ad file')
    parser.add_argument('--target-sum', type=float, default=1e4,
                       help='Target sum for normalization (default: 1e4)')
    parser.add_argument('--n-top-genes', type=int, default=2000,
                       help='Number of highly variable genes (default: 2000)')
    parser.add_argument('--hvg-method', default='seurat', choices=['seurat', 'cell_ranger', 'seurat_v3'],
                       help='Method for identifying highly variable genes (default: seurat)')
    parser.add_argument('--regress', nargs='+', default=['total_counts', 'pct_counts_mt'],
                       help='Variables to regress out (default: total_counts pct_counts_mt)')
    parser.add_argument('--no-regress', action='store_true',
                       help='Skip regression step')
    parser.add_argument('--max-value', type=float, default=10,
                       help='Maximum value for scaling (default: 10)')
    parser.add_argument('--plot-prefix', default='preprocess',
                       help='Prefix for output plots (default: preprocess)')
    
    args = parser.parse_args()
    
    # Setup
    setup_scanpy()
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Store raw counts
    adata.raw = adata.copy()
    
    # Normalize
    adata = normalize_data(adata, args.target_sum)
    
    # Log transform
    adata = log_transform(adata)
    
    # Identify highly variable genes
    adata = identify_hvgs(adata, args.n_top_genes, args.hvg_method, 
                         plot=True, output_prefix=args.plot_prefix)
    
    # Regress out unwanted variation (optional)
    if not args.no_regress:
        adata = regress_out_effects(adata, args.regress)
    
    # Scale data
    adata = scale_data(adata, args.max_value)
    
    # Save processed data
    print(f"\nSaving processed data to {args.output}...")
    adata.write(args.output)
    
    print("\nPreprocessing completed successfully!")
    print(f"Processed dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Highly variable genes: {sum(adata.var['highly_variable'])}")

if __name__ == '__main__':
    main()
