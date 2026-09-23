# 05 Differential analysis

Pseudobulk DESeq2 by cell type (and hepatocyte sub-type) for each modality, comparing MASL, MASH, and MetALD to normal livers, fibrosis stage (quantitative and low/high), and steatosis grade. Covariates are scaled age, scaled BMI, sex, and pool. fGSEA uses KEGG, Reactome, and GO.

## File descriptions

### `1_Pseudobulk_TPM_CPM.ipynb`
Pseudobulks counts by cell type and donor, and calculates TPM for RNA and CPM for ATAC, H3K27ac, and H3K27me3. Used for marker plots (Figure 1C) and QTL inputs.

### `240829_WE_Liver_RNA_DEseq.ipynb`
DESeq2 and fGSEA for RNA in each cell type (Figures 3A and 3B).

### `240918_MK_Liver_RNA_DEseq_Hepatocytes_Subtypes.ipynb`
DESeq2 and fGSEA for RNA in each hepatocyte sub-type (Figures 6G and 6H).

### `240925_MK_Calculate_per_Hepatocyte_Cellsubtype_Entropy2_Clean.ipynb`
Finds hepatocyte sub-type marker genes from pseudobulk TPM. It converts TPM to proportions across the five sub-types and computes entropy with DescTools `Entropy`. A gene counts as a marker for a sub-type if it has low entropy, its maximum TPM is in that sub-type, and TPM > 1. It exports the unique marker genes per sub-type and the top 500 lowest-entropy genes (Figure 6B).

### `241022_WE_Liver_ATAC_DEseq.ipynb`
DESeq2 for ATAC in each cell type (Figure 3A).

### `241024_WE_Liver_H3K27ac_DEseq_ATAC_Peaks.ipynb`
DESeq2 for H3K27ac, quantified in ATAC peaks, for each cell type (Figure 3A).

### `241219_WE_Liver_H3K27me3_DEseq_Bins.ipynb`
DESeq2 for H3K27me3, quantified in 15 kb bins, for each cell type (Figure 3A).

### `241227_WE_Homer_DAC_cluster_sites.ipynb`
HOMER on clusters of sites changing in disease (Figure 3D).

### `250117_WE_Liver_H3K27me3_DEseq_Bins_Hepatocyte_Subtypes.ipynb`
DESeq2 for H3K27me3 bins in each hepatocyte sub-type, plus the overview of sub-type contrasts (Figures 6G and 6H).

### `250212_WE_Liver_ATAC_DEseq_Subtypes_Celltype_Peaks.ipynb`
DESeq2 for ATAC in each hepatocyte sub-type (Figures 6G and 6I).

### `250224_WE_Liver_H3K27ac_DEseq_Subtypes.ipynb`
DESeq2 for H3K27ac in each hepatocyte sub-type (Figure 6G).
