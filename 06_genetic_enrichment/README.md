# 06 Genetic enrichment

LDSC partitioned heritability and FINRICH enrichment of fine-mapped cALT (MVP) credible set variants, plus HOMER on the same site sets. For FINRICH, each enrichment test is written as a separate command, batch-submitted to the HPC, and the results are read back in.

## File descriptions

### `240924_WE_FINRICH_intput_prep.ipynb`
Lifts the GWAS credible sets to hg38 and writes feature and background sets in bed format.

### `241025_WE_Plot_MVP_GWAS_Clean.ipynb`
Re-plots the MVP cALT GWAS Manhattan plot.

### `241217_WE_Prep_FINRICH_For_TSCC.ipynb`
FINRICH on differential sites from DESeq2.

### `241221_WE_Prep_FINRICH_For_TSCC_states_GRNs_H3K27me3_bins.ipynb`
FINRICH on chromatin state clusters, SCENIC+ GRNs, and differential H3K27me3 15 kb bins (Figure 3I).

### `250219_WE_FINRICH_and_Homer_For_TSCC_DAR_Clusters_High_Fibrosis.ipynb`
FINRICH and HOMER on clustered differentially accessible cCREs in high fibrosis (Figures 3D and 6I).

### `250220_WE_Prep_FINRICH_For_TSCC_GRNs.ipynb`
FINRICH on clustered GRN cCREs (Figure 3I).

### `250408_WE_FINRICH_and_Homer_For_TSCC_High_Fibrosis_and_cRE_States.ipynb`
FINRICH and HOMER on cCREs by chromatin state and on sites changing in accessibility in high fibrosis (Figures 2E and 6I).

### `250428_WE_FINRICH_and_Homer_For_Clustered_Marker_Peaks.ipynb`
FINRICH and HOMER on NMF cCRE modules (Figure 2D).

### `LDSC.sh`
Example of how LDSC was run. Peaks were lifted to hg19. The steps are: munge the GWAS, make annotations, compute LD scores, then run partitioned heritability with `--h2-cts`.

### `241221_WE_LDSC_Plotting.ipynb`
Plots LDSC enrichment for cell types, modalities, sub-types, and chromatin states (Figure 2E).
