#!/usr/bin/env python3
"""
Single-cell RNA-seq Quality Control Script
单细胞RNA-seq数据质控脚本

This script performs quality control on single-cell RNA-seq data including:
- Filtering cells based on gene counts
- Filtering cells based on mitochondrial gene percentage
- Filtering genes based on cell expression
- Generating QC plots

Usage: python scrna_qc.py <input_matrix> <output_h5ad> [options]
"""

import argparse
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def setup_scanpy():
    """Configure scanpy settings"""
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_data(input_file, data_format='auto'):
    """
    Load single-cell RNA-seq data
    
    Parameters:
    -----------
    input_file : str
        Path to input file (h5ad, csv, tsv, or 10x format)
    data_format : str
        Format of input data ('auto', 'h5ad', 'csv', '10x')
    
    Returns:
    --------
    adata : AnnData
        Annotated data matrix
    """
    if data_format == 'auto':
        if input_file.endswith('.h5ad'):
            adata = sc.read_h5ad(input_file)
        elif input_file.endswith('.csv'):
            adata = sc.read_csv(input_file)
        elif input_file.endswith('.tsv') or input_file.endswith('.txt'):
            adata = sc.read_csv(input_file, delimiter='\t')
        else:
            # Try 10x format
            adata = sc.read_10x_mtx(input_file, var_names='gene_symbols', cache=True)
    elif data_format == 'h5ad':
        adata = sc.read_h5ad(input_file)
    elif data_format == 'csv':
        adata = sc.read_csv(input_file)
    elif data_format == '10x':
        adata = sc.read_10x_mtx(input_file, var_names='gene_symbols', cache=True)
    else:
        raise ValueError(f"Unsupported data format: {data_format}")
    
    return adata

def calculate_qc_metrics(adata, species='human'):
    """
    Calculate QC metrics for cells and genes
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    species : str
        Species ('human' or 'mouse') for mitochondrial gene detection
    
    Returns:
    --------
    adata : AnnData
        Updated annotated data with QC metrics
    """
    # Identify mitochondrial genes
    # Note: This uses case-insensitive matching to catch MT-, Mt-, mt- variants
    if species.lower() == 'human':
        adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    elif species.lower() == 'mouse':
        adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    else:
        # For other species, try both common patterns
        adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    return adata

def filter_cells_and_genes(adata, min_genes=200, max_genes=None, min_counts=None, 
                          max_counts=None, max_mt_percent=20, min_cells=3):
    """
    Filter cells and genes based on QC metrics
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    min_genes : int
        Minimum number of genes expressed per cell
    max_genes : int
        Maximum number of genes expressed per cell (None for no limit)
    min_counts : int
        Minimum number of counts per cell (None for no limit)
    max_counts : int
        Maximum number of counts per cell (None for no limit)
    max_mt_percent : float
        Maximum percentage of mitochondrial genes
    min_cells : int
        Minimum number of cells expressing a gene
    
    Returns:
    --------
    adata : AnnData
        Filtered annotated data
    """
    print(f"\nBefore filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes is not None:
        adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    if min_counts is not None:
        adata = adata[adata.obs.n_counts >= min_counts, :]
    if max_counts is not None:
        adata = adata[adata.obs.n_counts < max_counts, :]
    
    # Filter by mitochondrial percentage
    adata = adata[adata.obs.pct_counts_mt < max_mt_percent, :]
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    return adata

def plot_qc_metrics(adata, output_prefix):
    """
    Generate QC plots
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    output_prefix : str
        Prefix for output plot files
    """
    # Violin plots for QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save=f'_{output_prefix}_violin.png')
    
    # Scatter plots
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save=f'_{output_prefix}_counts_vs_mt.png')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save=f'_{output_prefix}_counts_vs_genes.png')
    
    print(f"\nQC plots saved with prefix: {output_prefix}")

def main():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Quality Control')
    parser.add_argument('input', help='Input data file (h5ad, csv, tsv, or 10x directory)')
    parser.add_argument('output', help='Output h5ad file')
    parser.add_argument('--format', default='auto', choices=['auto', 'h5ad', 'csv', '10x'],
                       help='Input data format (default: auto)')
    parser.add_argument('--species', default='human', choices=['human', 'mouse', 'other'],
                       help='Species for mitochondrial gene detection (default: human)')
    parser.add_argument('--min-genes', type=int, default=200,
                       help='Minimum number of genes per cell (default: 200)')
    parser.add_argument('--max-genes', type=int, default=None,
                       help='Maximum number of genes per cell (default: None)')
    parser.add_argument('--min-counts', type=int, default=None,
                       help='Minimum number of counts per cell (default: None)')
    parser.add_argument('--max-counts', type=int, default=None,
                       help='Maximum number of counts per cell (default: None)')
    parser.add_argument('--max-mt-percent', type=float, default=20,
                       help='Maximum mitochondrial gene percentage (default: 20)')
    parser.add_argument('--min-cells', type=int, default=3,
                       help='Minimum number of cells per gene (default: 3)')
    parser.add_argument('--plot-prefix', default='qc',
                       help='Prefix for QC plot files (default: qc)')
    
    args = parser.parse_args()
    
    # Setup
    setup_scanpy()
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = load_data(args.input, args.format)
    
    # Calculate QC metrics
    print("\nCalculating QC metrics...")
    adata = calculate_qc_metrics(adata, args.species)
    
    # Generate QC plots before filtering
    print("\nGenerating QC plots...")
    plot_qc_metrics(adata, f"{args.plot_prefix}_before")
    
    # Filter cells and genes
    print("\nFiltering cells and genes...")
    adata = filter_cells_and_genes(
        adata,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        min_counts=args.min_counts,
        max_counts=args.max_counts,
        max_mt_percent=args.max_mt_percent,
        min_cells=args.min_cells
    )
    
    # Generate QC plots after filtering
    plot_qc_metrics(adata, f"{args.plot_prefix}_after")
    
    # Save filtered data
    print(f"\nSaving filtered data to {args.output}...")
    adata.write(args.output)
    
    print("\nQuality control completed successfully!")
    print(f"Final dataset: {adata.n_obs} cells × {adata.n_vars} genes")

if __name__ == '__main__':
    main()
