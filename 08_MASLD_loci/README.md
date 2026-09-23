# 08 MASLD risk loci

Annotation of fine-mapped MASLD (cALT) loci using hepatocyte QTLs, cCREs, chromatin states, chromatin loops, GRNs, and chromBPNet.

## File descriptions

### `240930_WE_Coloc_GWAS.ipynb`
Converts the Vujkovic et al. GWAS to hg38 and colocalizes QTLs with the GWAS (coloc.abf, PPH4>0.8) (Figure 5D).

### `241219_QTL_Summary_Stats_Intersect_Finemapping.ipynb`
Finds fine-mapped GWAS variants that are also QTL variants.

### `250124_WE_Hep_cREs_Intersect.ipynb`
Overlap of fine-mapped variants with cCREs and active states in each cell type (Figure 2E, bottom right).

### `250303_WE_Prep_ChrmoBP_NarrowPeak.ipynb`
Writes narrowPeak-style regions around credible set variants for chromBPNet.

### `chrombpnet.variant.effects.sh`
Makes predicted count and contribution bigwigs for each allele of a target variant. It inserts the variant into the reference genome, then runs chromBPNet on both the reference and alternate genomes (Figures 5G and 5I).

### `250325_WE_annotating_gwas_v2.ipynb`
Builds the locus annotation table from summary stats and credible sets, coloc results, cCRE and active-state overlap, and chromBPNet effects (Figure 5E).

### `250303_WE_Plot_QTLs.ipynb`
QTL boxplots by genotype, locus zoom plots, and mashr effect plots for KRT8, PPP1R3B, EFHD1, CEBPA, and XBP1 (Figures 5F–I).

### `241206_WE_Quick_Interaction_Test.ipynb`
Test of genotype by disease interaction at the KRT8 locus.
