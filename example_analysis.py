#!/usr/bin/env python3
"""
Example script demonstrating single-cell RNA-seq analysis
单细胞RNA-seq分析示例脚本

This script shows how to use the scRNA-seq analysis tools step by step
"""

import scanpy as sc
import numpy as np
import os

def main():
    """Run the example analysis"""
    # This example uses the built-in PBMC dataset from scanpy
    print("="*60)
    print("Single-cell RNA-seq Analysis Example")
    print("单细胞RNA-seq分析示例")
    print("="*60)

    # Create output directory
    output_dir = "example_output"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)

    sc.settings.figdir = f"{output_dir}/figures"
    sc.settings.set_figure_params(dpi=150)

    print("\nStep 1: Loading example data...")
    print("步骤1：加载示例数据...")

    # Load example data (PBMC 3k dataset)
    # This will download the data if not already cached
    adata = sc.datasets.pbmc3k()
    print(f"Loaded {adata.n_obs} cells × {adata.n_vars} genes")

    print("\nStep 2: Quality control...")
    print("步骤2：质量控制...")

    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Plot QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save='_qc_before.png')

    # Filter cells and genes
    print(f"Before filtering: {adata.n_obs} cells")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    print(f"After filtering: {adata.n_obs} cells")

    print("\nStep 3: Normalization...")
    print("步骤3：标准化...")

    # Store raw counts
    adata.raw = adata.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print(f"Found {sum(adata.var['highly_variable'])} highly variable genes")
    sc.pl.highly_variable_genes(adata, save='_hvg.png')

    # Regress out effects and scale
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    print("\nStep 4: Dimensionality reduction...")
    print("步骤4：降维分析...")

    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True, save='_variance.png')

    # Compute neighbors
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    # UMAP
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['n_genes_by_counts', 'pct_counts_mt'], save='_qc.png')

    print("\nStep 5: Clustering...")
    print("步骤5：聚类分析...")

    # Leiden clustering
    sc.tl.leiden(adata, resolution=0.9)
    print(f"Found {len(adata.obs['leiden'].unique())} clusters")

    # Plot clusters
    sc.pl.umap(adata, color='leiden', save='_clusters.png', legend_loc='on data')

    print("\nStep 6: Finding marker genes...")
    print("步骤6：寻找标记基因...")

    # Find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_markers.png')

    # Plot marker genes
    marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                    'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                    'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']

    # Filter to genes that exist in dataset
    available_markers = [g for g in marker_genes if g in adata.var_names]
    print(f"Plotting {len(available_markers)} marker genes")

    sc.pl.dotplot(adata, available_markers, groupby='leiden', save='_known_markers.png')
    sc.pl.umap(adata, color=available_markers[:9], save='_marker_expression.png')

    print("\nStep 7: Saving results...")
    print("步骤7：保存结果...")

    # Save processed data
    adata.write(f"{output_dir}/pbmc3k_processed.h5ad")

    # Export marker genes
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    with open(f"{output_dir}/marker_genes.tsv", 'w') as f:
        f.write("cluster\tgene\tscore\tpval\tpval_adj\tlogfoldchange\n")
        for group in groups:
            for i in range(min(25, len(result['names'][group]))):
                gene = result['names'][group][i]
                score = result['scores'][group][i]
                pval = result['pvals'][group][i]
                pval_adj = result['pvals_adj'][group][i]
                logfc = result['logfoldchanges'][group][i]
                f.write(f"{group}\t{gene}\t{score:.4f}\t{pval:.4e}\t{pval_adj:.4e}\t{logfc:.4f}\n")

    print("\n" + "="*60)
    print("Analysis completed successfully!")
    print("分析完成！")
    print("="*60)
    print(f"\nResults saved to: {output_dir}/")
    print(f"Processed data: {output_dir}/pbmc3k_processed.h5ad")
    print(f"Marker genes: {output_dir}/marker_genes.tsv")
    print(f"Figures: {output_dir}/figures/")
    print(f"\nFinal dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Number of clusters: {len(adata.obs['leiden'].unique())}")


if __name__ == '__main__':
    main()
