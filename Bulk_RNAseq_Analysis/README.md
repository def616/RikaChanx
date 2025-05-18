# Bulk RNA-seq Analysis

**Tools Used:** DESeq2, edgeR  
**Goal:** Identify differentially expressed genes in brown adipose tissue between wild-type and Hoxa5 knockout mice.

## Methods
- Preprocessed count data.
- Performed differential expression analysis.
- Conducted pathways analysis.

## Folder structure
- `Bulk_RNAseq_Analysis/`
  - `README.md`: Description of the analysis
  - `deseq2.Rmd` and `edgeR.Rmd`: Scripts
  - `y_chrom_removed.R`: Analysis after the removal of Y-linked genes
  - `result_writeup.pdf`: A short paper on the results and methods used
  - `requirements.txt`: Environment setup

## Description
This project investigates differential gene expression in brown adipose tissue (BAT) between wild-type and Hoxa5 knockout mice using both DESeq2 and edgeR. While edgeR identified more differentially expressed genes (DEGs) than DESeq2, many of these were driven by Y-linked genes such as Ddx3y—highlighting sex-specific biases in the data. After excluding Y-linked genes, only the control Hoxa5 remained significantly differentially expressed, suggesting their strong influence on initial results.

The analysis reflects known limitations of DEG methods in small-sample RNA-seq studies and raises questions about the interplay between Hoxa5 and Y-linked gene expression. This work points toward a possible sex-specific regulation of BAT development, with implications for understanding metabolic differences between male and female mice.




