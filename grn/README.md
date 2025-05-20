# Gene Regulatory Network (GRN) analysis

**Tools Used:** SCENIC, GRNBoost2, CellOracle, SCopeLoomR
**Goal:** Identify potential Brown Adipose Tissue (BAT) progenitor cells population using single cell data

## Methods
- Preprocessed single cell dataset: filter low-expressing genes, mitochondrial genes, and doublets
- SCENIC with GRNBoost2 was used to infer gene regulatory networks; outputs were processed in R using SCopeLoomR
- Filtered dataset was analyzed using CellOracle to infer cluster-specific GRNs
- 
## Folder structure
- `GRN/`
    - `grn_code_cleaned`: code used perform gene regulatory network inference
    - `results_writeup.pdf`: a short paper on the results obtained

## Description
Characterizing progenitors of Brown Adipose Tissue (BAT) is crucial to studying BAT development. SCENIC and CellOracle were used to inferred gene regulatory networks (GRNs), and identify possible
BAT progenitor cell types in single cell data of mouse tissue taken at E12.5. Cxcl14 was identified as a marker for a potential BAT progenitor cell population.





