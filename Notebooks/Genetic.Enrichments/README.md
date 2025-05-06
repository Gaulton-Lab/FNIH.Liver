# Multiome Processing Notebooks

This directory contains Jupyter notebooks to run DESeq, fGSEA, and homer.

## File Descriptions

### `240924_WE_FINRICH_intput_prep.ipynb`
**Description**: Lift the GWAS credible set to hg38 and format different features into bed format

### `241025_WE_Plot_MVP_GWAS_Clean.ipynb`
**Description**: The Manhattan plot from the MVP GWAS used in this analysis is not powerpoint friendly; replotting here.

### `241217_WE_Prep_FINRICH_For_TSCC.ipynb`
**Description**: For FINRICH we are writing each enrichment test as a distinct command and batch submitting to the HPC before reading in the results. 

### `241221_WE_LDSC_Plotting.ipynb`
**Description**: Plotting LDSC enirchment for cell types and modalities.

### `241221_WE_Prep_FINRICH_For_TSCC_states_GRNs_H3K27me3_bins.ipynb`
**Description**: For FINRICH we are writing each enrichment test as a distinct command and batch submitting to the HPC before reading in the results. This is enrichment on the differentially methylated 15kb bins from DESeq analysis.

### `250124_WE_Hep_cREs_Intersect.ipynb`
**Description**: In this notebook we examine the overlap of NAFLD variants with cREs by cell type and active regions by cell type.

### `250219_WE_FINRICH_and_Homer_For_TSCC_DAR_Clusters_High_Fibrosis.ipynb`
**Description**: For FINRICH we are writing each enrichment test as a distinct command and batch submitting to the HPC before reading in the results. This enrichment is looking at clustered differentially accessible cREs.

### `250220_WE_Prep_FINRICH_For_TSCC_GRNs.ipynb`
**Description**: For FINRICH we are writing each enrichment test as a distinct command and batch submitting to the HPC before reading in the results. This enrichment is looking at clustered GRN cREs.

### `250303_WE_Prep_ChrmoBP_NarrowPeak.ipynb`
**Description**: Based on the credible sets this makes a narrow peak like format of regions for running chromBPNet.

### `250408_WE_FINRICH_and_Homer_For_TSCC_High_Fibrosis_and_cRE_States.ipynb`
**Description**: For FINRICH we are writing each enrichment test as a distinct command and batch submitting to the HPC before reading in the results. This enrichment is looking at cREs by states and sites changing in accessibility (ATAC) in high fibrosis.

### `250428_WE_FINRICH_and_Homer_For_Clustered_Marker_Peaks.ipynb`
**Description**: For FINRICH we are writing each enrichment test as a distinct command and batch submitting to the HPC before reading in the results. This enrichment is looking at the clustering of cREs by accessibility in various cell types into modules.

### `chrombpnet.variant.effects.sh`
**Description**: This script is how we made bigwigs for each allele of a target variant. First we inject the variant into the reference genome, then run chromBPNet on both the reference and alternate reference genome.

### `LDSC.sh`
**Description**: This scripts is an example of how LDSC was run for this project. LDSC was run using hg19 reference files so peak annotations were lifted back to hg19. The GWAS was munged, then annotations made, LD score regression run, and finally heritability was partitioned using the h2-cts flag.