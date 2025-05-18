# Bulk RNA-seq Analysis

**Tools Used:** DESeq2, edgeR  
**Goal:** Identify differentially expressed genes in brown adipose tissue between wild-type and Hoxa5 knockout mice.

## Methods
- Normalized count data using DESeq2.
- Performed differential expression analysis.
- Conducted downstream GO enrichment analysis.

> 🔐 Note: Due to data confidentiality, results are based on synthetic data resembling real distribution patterns.

## Script
See `DEseq2.Rmd` and `edgeR.Rmd` for the workflow.

## Sample Output
- Volcano plot of DEGs
- Heatmap of top 50 DEGs

## Folder structure
- `Bulk_RNAseq_Analysis/`
  - `README.md`: Description of the analysis
  - `DEseq2.Rmd` or `edgeR.Rmd`: Scripts
  - `result_writeup.pdf`: A short paper on the results and methods used
  - `requirements.txt`: Environment setup




