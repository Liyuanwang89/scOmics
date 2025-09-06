# scOmics

Single-Cell Multi-Omics Data Analysis Pipeline

0. Requirements
The following environments/packages are required to run the scripts.

Single cell data analysis (scRNA_processing.r and scATAC_processing.r)
R v4.1.2

ArchR v1.0.1

Seurat v4.1.0

dplyr v1.0.4

patchwork v1.1.1

tidyverse v1.3.1+

clusterProfiler v4.0+

TxDb.Sscrofa.UCSC.susScr11.refGene

org.Ss.eg.db

BSgenome.Sscrofa.UCSC.susScr11

scran v1.20.0+

RColorBrewer v1.1-2

ggsci v2.9

pheatmap v1.0.12

Operating System
Linux (Ubuntu 18.04+ or CentOS 7+ recommended)

Memory: Minimum 32GB, Recommended 64GB+

Storage: At least 500GB available space

1. Single nucleus ATAC-seq (snATAC-seq) data analysis
All the scripts for the snATAC-seq data analysis are included in the scATAC_processing.r file.

The snATAC-seq data were generated for this study and pre-processed with 10x Genomics Cell Ranger ATAC pipeline (v1.2.0, with default parameters). The output files from the Cell Ranger ATAC pipeline were used as input for our analysis. The script contains the following parts:

Quality control and filtering: Cells were filtered based on TSS enrichment score (>4) and number of unique fragments (>1000)

Dimensionality reduction: Iterative LSI dimensionality reduction with 50,000 variable features and 25 dimensions

Batch correction: Harmony integration to correct for batch effects across samples

Clustering: Seurat-based clustering with resolution 0.9

Visualization: UMAP and t-SNE embeddings

Gene scoring: Calculation of gene activity scores from chromatin accessibility

Marker identification: Wilcoxon rank sum test for identifying cell-type-specific features

Motif enrichment: JASPAR 2016 database annotation and ChromVAR analysis

Integration with scRNA-seq: Label transfer using FindTransferAnchors

2. Single cell RNA-seq (scRNA-seq) data analysis
All the scripts for the scRNA-seq data analysis are included in the scRNA_processing.r file.

The scRNA-seq data span 21 porcine tissues across four developmental stages (E65, D1, Y1, Y3). The analysis provides a comprehensive reference atlas for integration with snATAC-seq data. The script contains the following parts:

Quality control: Filtering cells with <500 genes, <1000 UMIs, or >20% mitochondrial genes

Ambient RNA removal: DecontX implementation for removing contamination

Doublet detection: DoubletFinder with default parameters

Normalization and scaling: Log2 normalization and data scaling

Integration: Harmony batch correction across developmental stages

Clustering: FindNeighbors and FindClusters with resolution 0.5

Marker gene identification: FindAllMarkers with min.pct = 0.25 and logfc.threshold = 0.25

Cell type annotation: SCSA and SingleR with manual curation based on marker genes

3. Multi-omics integration analysis
The integration analysis scripts are included in the integration_analysis.r file, performing joint analysis of scRNA-seq and snATAC-seq data using both Seurat and ArchR frameworks.

The pipeline includes:

Data integration: Using Seurat's FindTransferAnchors function with canonical correlation analysis

Label transfer: Predicting cell types from scRNA-seq to snATAC-seq data

Co-embedding: Merging scRNA-seq and snATAC-seq datasets in shared embedding space

Multi-omics visualization: UMAP and t-SNE plots showing both modalities

Confusion matrix analysis: Evaluating clustering consistency between modalities
