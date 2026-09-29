# 07 QTL mapping

Code to run tensorQTL for RNA (eQTL), ATAC (caQTL), H3K27ac (acQTL), and H3K27me3 (me3QTL) in each cell type, plus downstream specificity, colocalization, and motif analyses.

## File descriptions

### `240829_WE_Genotype_PCA_Visualization.ipynb`
Genotype PCA with 1000 Genomes (run in plink). The PCs are used as ancestry covariates.

### `240902_WE_Normalize_Transform_PCA.ipynb`
TMM normalization and inverse-normal transformation of the count data for each modality. Also computes feature PCs, used as covariates for latent factors.

### `240911_WE_TensorQTL_Inputs.ipynb`
Writes phenotype bed files, covariates, and the VCF, and prepares the tensorQTL SLURM commands.

### `240912_WE_*QTL_Summary.ipynb`
Aggregate variant-level summary stats, add per-feature FDR cutoffs and significance flags, and count significant features per cell type.

### `241122_WE_Plotting_Summaries.ipynb`
Plots the number of significant QTL features across modalities and cell types (Figure 5A).

### `241122_WE_Specificity_*QTLs.ipynb`
mashr re-estimation of QTL effects across cell types and QTL cell type specificity for all lead variants (Figure 5C, bottom).

### `240930_WE_Coloc_QTLs.ipynb`
Colocalization of QTLs across modalities and construction of cross-modality QTL modules with igraph/Louvain (Figure 5C, top). The coloc itself was run on the HPC.

### `241206_WE_MotifbreakR.ipynb`
motifbreakR on QTL variants, with an aggregate enrichment of disrupted motifs per cell type and modality (Figure 5B).
