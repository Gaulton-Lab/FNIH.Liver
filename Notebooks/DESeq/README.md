# Multiome Processing Notebooks

This directory contains Jupyter notebooks to run DESeq, fGSEA, and homer.

## File Descriptions

### `1_Pseudobulk_TPM_CPM.ipynb`
**Description**: This notebook psuedobulks counts by cell type and donor and calculates TPM for RNA and CPM for ATAC and DPT histone modalities.

### `240829_WE_Liver_RNA_DEseq.ipynb`
**Description**: Performes DESeq and fGSEA analysis for the RNA modality for each cell type

### `240918_MK_Liver_RNA_DEseq_Hepatocytes_Subtypes.ipynb`
**Description**: Performes DESeq and fGSEA analysis for the RNA modality for each Hepatocyte subtype

### `241022_WE_Liver_ATAC_DEseq.ipynb`
**Description**: Performes DESeq analysis for the ATAC modality for each cell type

### `2241024_WE_Liver_H3K27ac_DEseq_ATAC_Peaks.ipynb`
**Description**: Performes DESeq analysis for the H3K27ac modality quantified in ATAC peaks for each cell type

### `241217_WE_Prep_HOMER_For_TSCC.ipynb`
**Description**: Prepare file inputs for running homer on the UCSD HPC. The commands were then run in batches and read back in and results compiled.

### `241219_WE_Liver_H3K27me3_DEseq_Bins.ipynb`
**Description**: Performes DESeq analysis for the H3K27me3 modality quantified in 15kb bins for each cell type

### `241227_WE_Homer_DAC_cluster_sites.ipynb`
**Description**: Homer run on the clustered sites changing in disease. 

### `250117_WE_Liver_H3K27me3_DEseq_Bins_Hepatocyte_Subtypes.ipynb`
**Description**: Performes DESeq analysis for the H3K27me3 modality quantified in 15kb bins for each Hepatocyte subtype

### `250212_WE_Liver_ATAC_DEseq_Subtypes_Celltype_Peaks.ipynb`
**Description**: Performes DESeq analysis for the ATAC modality for each Hepatocyte subtype

### `250224_WE_Liver_H3K27ac_DEseq_Subtypes.ipynb`
**Description**: Performes DESeq analysis for the H3K27ac modality quantified in ATAC peaks for each cell type