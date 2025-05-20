# Single cell RNA-seq Analysis

**Tools Used:** Scanpy, Velocyto, Cellrank \
**Goal:** Identify populations of cells and characterize signaling pathways that correspond to Brown Adipose Tissue (BAT) progenitor cells population

## Methods
- Preprocessed count data including normalization, logarithmization, mitochondrial genes detection, and doublet detection
- Performed PCA and UMAP for linear and non-linear dimensionality reduction, respectively
- Performed Leiden clustering algorithm for cell types clustering
- Enriching for BAT progenitor cells by selecting genes proposed to be responsible for the cell population
- Conducted pathways analysis using differentially expressed genes from the enriched cells
- Performed trajectory analysis using experimental times and RNA velocity

## Folder structure
`single_cell_analysis/` \
`single_ceLL_analysis.py`: a snippet of codes used to preprocess, analyze, and generate important figures \
`trajectory.py`: a snippet of codes used to perform trajectory inference using experimental time via Cellrank

## Description
Brown adipose tissue (BAT) plays a crucial role in thermogenesis and metabolism, yet the lineage of brown adipocyte progenitors (BAP) and brown adipocyte (BA) during embryonic development remains incompletely understood. Here, we investigate BAT development in mouse embryos from embryonic day E11.5 to E13.5 using single cell RNA sequencing (scRNA-seq). Single cell analysis of Ebf2+/Sox9+/Col2a1- cells reveal a distinct spatial arrangement of BAP and BA. Differentially expressed genes analysis reveal key adipogenic regulators such as Pparg and Cebpa. Trajectory inference of Pdgfra+ mesenchymal cells across the embryonic time points confirmed a lineage trajectory toward BAT, with Cxcl12 emerging as a potential driver via ERK signaling. Additional candidates including Fst and Id3 (BMP pathway) further support the involvement of multiple signaling pathways in BAT differentiation. We uncovered evidence of ERK and BMP signaling in BAT as well as transcriptional and spatial signatures of embryonic BAT development and propose regulatory genes for future validation. These insights lay the groundwork for mapping BAT development and understanding its implications in therapeutics for metabolic-related diseases.

**As this is an ongoing and unpublished work, dataset and results are not shared. Snippets of codes are shared.**





