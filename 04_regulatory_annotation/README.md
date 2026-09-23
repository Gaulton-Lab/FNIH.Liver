# 04 Regulatory annotation

### `Liver_analysis.ipynb`
Main analysis of the cCRE catalog:
- ABC links and marker gene plots
- comparison of cCREs to ENCODE SCREEN and HEA, and genomic feature annotation (Figure 2A)
- cell type chromatin states from ChromHMM, and chromatin state of each cCRE (Figure 2B)
- cell type- and sub-type-specific cCREs by Jensen-Shannon specificity, using the SnapATAC2 LR test results (Figures 2F and 6E)
- NMF cCRE modules and their chromatin states (Figure 2D)
- SCENIC+ cell type-specific links, target gene expression, and APA of links in hepatocytes vs. myeloid cells (Figure 2C)
- motif and GO enrichment of specific cCREs (Figures 2G and 2H)
- k-means clustering of hepatocyte cCREs changing in Fib+ MASH, with their chromatin states (Figures 3D and 3E)

### `241217_WE_Prep_HOMER_For_TSCC.ipynb`
Prepares and reads back HOMER `findMotifsGenome.pl` runs (JASPAR 2022) for cell type markers, cell type-specific cCREs, sub-type markers, differential peaks, and clusters. The commands were run in batches on TSCC.

### `241224_WE_Homer_fGSEA_on_specific_links.ipynb`
fGSEA (GO:BP) of genes linked to cell type-specific cCREs, and HOMER on the linked cCREs (Figures 2G and 2H).
