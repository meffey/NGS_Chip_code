# NGS Analysis Tools / NGS分析工具集

This repository contains tools for analyzing both ChIP-seq and single-cell RNA-seq data.

本仓库包含ChIP-seq和单细胞RNA-seq数据分析工具。

---

## ChIP-seq Analysis Tools (Perl Scripts)

### bed_overlap.pl:
   Description: This script find overlap peaks between two peaks. bed format is input type.
#### Usage:
```bash
perl bed_overlap.pl input_a.bed input_b.bed out.bed
```

### extract_seq.pl：
   Description: This script extract equally length sequence according to peaks.
#### Usage:
```bash
perl extract_seq.pl genome.fa peaks.bed length out.fa
```
Example:
```bash
perl extract_seq.pl A.fa A_peaks.bed 1000 A_1000.fa
```

### peaks_annotation.pl:
   Description: This script annotated peaks according to relationship between peaks and gene(include exon,intron,3'UTR,5'UTR and promote).
#### Usage:
```bash
perl peaks_annotation.pl genes.gff peaks.bed
```

### peaks_annotation_tss.pl:
   Description: This script annotated peaks according to relationship between peaks and transcription start sites(TSS-400bp).
#### Usage:
```bash
perl peaks_annotation_tss.pl genes_transposons.gff peaks.bed
```

---

## Single-cell RNA-seq Analysis Tools (Python Scripts)
## 单细胞RNA-seq分析工具（Python脚本）

### Installation / 安装依赖

First, install the required Python packages:
首先安装所需的Python包：

```bash
pip install -r requirements.txt
```

### Complete Analysis Pipeline / 完整分析流程

**scrna_pipeline.py** - Complete single-cell RNA-seq analysis pipeline
完整的单细胞RNA-seq分析流程

This script performs the complete analysis workflow including QC, normalization, dimensionality reduction, clustering, and marker gene identification.
此脚本执行完整的分析工作流程，包括质控、标准化、降维、聚类和标记基因识别。

#### Usage / 使用方法:
```bash
# Analyze 10x Genomics data
python scrna_pipeline.py /path/to/10x_data output_dir --format 10x

# Analyze h5ad file
python scrna_pipeline.py data.h5ad output_dir --format h5ad

# Customize parameters for mouse data
python scrna_pipeline.py data.h5ad output_dir --species mouse --min-genes 300 --resolution 0.5
```

#### Key Parameters / 主要参数:
- `--species`: Species (human/mouse/other) - 物种
- `--min-genes`: Minimum genes per cell (default: 200) - 每个细胞最少基因数
- `--max-mt-percent`: Maximum mitochondrial percentage (default: 20) - 最大线粒体基因百分比
- `--resolution`: Clustering resolution (default: 1.0) - 聚类分辨率
- `--n-top-genes`: Number of highly variable genes (default: 2000) - 高变基因数量

#### Outputs / 输出:
- `processed_data.h5ad`: Processed data with all analysis results
- `marker_genes.tsv`: Marker genes for each cluster
- `cluster_assignments.tsv`: Cell cluster assignments
- `figures/`: All generated plots and visualizations

---

### Individual Analysis Steps / 独立分析步骤

#### 1. Quality Control / 质量控制
**scrna_qc.py** - Perform quality control and filtering

```bash
python scrna_qc.py input_data.h5ad qc_filtered.h5ad \
    --species human \
    --min-genes 200 \
    --max-mt-percent 20
```

#### 2. Normalization and Preprocessing / 标准化和预处理
**scrna_normalize.py** - Normalize and preprocess data

```bash
python scrna_normalize.py qc_filtered.h5ad normalized.h5ad \
    --n-top-genes 2000 \
    --hvg-method seurat
```

#### 3. Dimensionality Reduction / 降维分析
**scrna_reduce_dims.py** - Perform PCA, t-SNE, and UMAP

```bash
python scrna_reduce_dims.py normalized.h5ad reduced.h5ad \
    --n-pcs 50 \
    --run-umap \
    --run-tsne
```

#### 4. Clustering / 聚类分析
**scrna_cluster.py** - Cluster cells and identify marker genes

```bash
python scrna_cluster.py reduced.h5ad clustered.h5ad \
    --method leiden \
    --resolution 1.0 \
    --export-markers markers.tsv
```

#### 5. Differential Expression / 差异表达分析
**scrna_diff_expr.py** - Find differentially expressed genes between groups

```bash
python scrna_diff_expr.py clustered.h5ad de_results \
    --groupby leiden \
    --group1 0 \
    --group2 1
```

#### 6. Visualization / 可视化
**scrna_visualize.py** - Generate comprehensive visualizations

```bash
python scrna_visualize.py clustered.h5ad visualization_output \
    --genes CD3D CD8A CD4 \
    --all
```

---

## Example Workflow / 示例工作流程

### Complete pipeline in one command / 一键完整流程:
```bash
python scrna_pipeline.py my_data.h5ad results/ \
    --species human \
    --min-genes 200 \
    --max-mt-percent 20 \
    --resolution 1.0
```

### Step-by-step analysis / 分步分析:
```bash
# Step 1: Quality control
python scrna_qc.py raw_data.h5ad qc_data.h5ad

# Step 2: Normalization
python scrna_normalize.py qc_data.h5ad norm_data.h5ad

# Step 3: Dimensionality reduction
python scrna_reduce_dims.py norm_data.h5ad reduced_data.h5ad --run-umap

# Step 4: Clustering
python scrna_cluster.py reduced_data.h5ad clustered_data.h5ad

# Step 5: Visualization
python scrna_visualize.py clustered_data.h5ad plots/ --all
```

---

## Input Data Formats / 输入数据格式

The scripts support multiple input formats:
脚本支持多种输入格式：

- **10x Genomics**: Directory containing `matrix.mtx`, `genes.tsv`, `barcodes.tsv`
- **H5AD**: AnnData HDF5 format
- **CSV/TSV**: Gene expression matrix (genes × cells)

---

## Output Files / 输出文件

- **H5AD files**: Processed data at each step (可用于下一步分析)
- **TSV files**: Tables with marker genes, DE results, cluster assignments
- **PNG files**: High-resolution plots and visualizations
- **figures/**: Directory containing all generated figures

---

## Citation / 引用

If you use these tools in your research, please cite the underlying packages:
如果您在研究中使用这些工具，请引用以下软件包：

- **Scanpy**: Wolf, F. A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology.
- **AnnData**: https://anndata.readthedocs.io/

---

## Contact / 联系方式

For questions and issues, please open an issue on GitHub.
如有问题，请在GitHub上提交issue。
