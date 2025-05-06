# Multiome Processing Notebooks

This directory contains code to run tensorQTL and downstream testing based on QTL analysis.

## File Descriptions


### `240829_WE_Genotype_PCA_Visualization.ipynb`
**Description**: Genotype PCA was run using plink with the 1000 Genomes project data, which is used as a covariate in the QTL model to account for ancestry.

### `240902_WE_Normalize_Transform_PCA.ipynb`
**Description**: Normalize count data across modalities and generate feature PCs which are used as additional covariates to control for latent factors.

### `240911_WE_TensorQTL_Inputs.ipynb`
**Description**: Final reformatting for tensorQTL and preparing code to run on the HPC.

### `240912_WE_*QTL_Summary.ipynb`
**Description**: Aggregate variant level summary stats and add per feature FDR cutoffs and flags to each row. Also, generate summary of number of features that are significant by cell type. acQTL=H3K27ac, eQTL=RNA, caQTL=ATAC, me3QTL=H3K27me3

### `240930_WE_Coloc_GWAS.ipynb`
**Description**: Code to colocalize QTLs with NAFLD GWAS signal. The actual coloc was run on HPC but this sets up the process and was copied over.

### `240930_WE_Coloc_QTLs.ipynb`
**Description**: We colocalized QTLs to identify any effects spanning modalities and create QTL modules. Actual colocalization was run on HPC.

### `241122_WE_Specificity_*QTLs.ipynb`
**Description**: Using mashr we restimated the QTL effects across cell types for each modality and identified QTL cell type specificity for all lead variants. acQTL=H3K27ac, eQTL=RNA, caQTL=ATAC, me3QTL=H3K27me3

### `241122_WE_Plotting_Summaries.ipynb`
**Description**: Plot summary of significant QTL features across modalities and cell types.

### `241206_WE_MotifbreakR.ipynb`
**Description**: MotifbreakR was used to identify motifs disrupted by QTL variants, with an aggregate enrichment performed to identify motifs enriched for each cell type and modality.

### `241206_WE_Quick_Interaction_Test.ipynb`
**Description**: For FINRICH we are writing each enrichment test as a distinct command and batch submitting to the HPC before reading in the results. This enrichment is looking at cREs by states and sites changing in accessibility (ATAC) in high fibrosis.

### `241219_QTL_Summary_Stats_Intersect_Finemapping.ipynb`
**Description**: Identify fine-mapped GWAS varaints that are also QTL variants.

### `250303_WE_Plot_QTLs.ipynb`
**Description**: This script is for making plots for QTLs, such as boxplots by genotype and locus zoom plots and mashr effect plots.

### `250325_WE_annotating_gwas_v2.ipynb`
**Description**: Annotation of NAFLD fine-mapped variants based on QTLs, cRE overlaps, and cRE-gene links. 