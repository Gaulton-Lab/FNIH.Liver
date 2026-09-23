# FNIH Liver  <img src="./images/FNIH_liver_color.png" align="right" width="200"/>
FNIH Liver is a single-nuclei, multi-modal epigenomic atlas of metabolic dysfunction-associated steatotic liver disease (MASLD).
<br/><br/>

This repository contains the code used for the manuscript: **Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis** ([medRxiv 10.1101/2025.05.09.25327043](https://www.medrxiv.org/content/10.1101/2025.05.09.25327043v1)).

## Abstract
Metabolic dysfunction-associated steatotic liver disease (MASLD) has limited treatments, and cell type-specific regulatory networks driving MASLD represent therapeutic avenues. We assayed five transcriptomic and epigenomic modalities in 2.4M cells from 86 livers across MASLD stages. Integrating modalities increased annotation of the genome in liver cell types several-fold over previous catalogs. We identified cell type regulatory networks of MASLD progression, including distinct hepatocyte networks driving MASL and mild and severe fibrosis MASH. Our single cell atlas annotated 88% of MASH-associated loci, including a third affecting hepatocyte regulation which we linked to distal target genes. Finally, we characterized hepatocyte heterogeneity, including MASH-enriched populations with altered repression, localization, and signaling. Overall, our results provide high-resolution maps of liver cell types and revealed novel targets for anti-MASH therapy.

## Repository layout

Directories are numbered in the order the analyses were run. Each directory has its own README describing every file.

| Directory | Contents | Main figures |
|---|---|---|
| [`01_preprocessing/`](01_preprocessing) | Genotyping/imputation, 10x multiome, Droplet Paired-Tag, Droplet Hi-C and Visium HD processing (incl. Tangram) | Methods |
| [`02_integration/`](02_integration) | Reference mapping of Paired-Tag and Hi-C onto multiome, joint UMAP, donor/QC summaries | 1B, 1C, 1F |
| [`03_composition/`](03_composition) | Cell type proportion changes (scCODA, miloR) | 1G, 6D |
| [`04_regulatory_annotation/`](04_regulatory_annotation) | cCREs, ChromHMM states, cell type-specific cCREs, NMF modules, motifs, SCENIC+ links | 2A–D, 2F–H, 3D–E |
| [`05_differential_analysis/`](05_differential_analysis) | Pseudobulk, DESeq2 and fGSEA across MASLD stages for cell types and hepatocyte sub-types | 3A–B, 6G–I |
| [`06_genetic_enrichment/`](06_genetic_enrichment) | LDSC partitioned heritability, FINRICH fine-mapped variant enrichment, HOMER on differential/clustered sites | 2E, 3I |
| [`07_QTL/`](07_QTL) | tensorQTL mapping, mashr specificity, cross-modality QTL coloc, motifbreakR | 5A–C |
| [`08_MASLD_loci/`](08_MASLD_loci) | GWAS–QTL colocalization, annotation of MASLD loci, chromBPNet variant effects, locus plots | 5D–I |
| [`09_hepatocyte_trajectory/`](09_hepatocyte_trajectory) | Hepatocyte union peaks and Monocle3 pseudotime | 6J |
| [`10_spatial_analysis/`](10_spatial_analysis) | Visium HD marker heatmaps, inflammatory hepatocyte module scores and neighborhood enrichment | 1I, 6K–L |
| [`envs/`](envs) | Conda environments for the Jupyter kernels used in the notebooks | – |

Related repositories:
- Peak calling pipeline: https://github.com/Gaulton-Lab/peak-call-pipeline

## Figure → code index

| Figure | Panel(s) | Code |
|---|---|---|
| 1 | B, C, F | `02_integration/Liver_integration.ipynb` (C also uses `05_differential_analysis/1_Pseudobulk_TPM_CPM.ipynb`) |
| 1 | G | `03_composition/260302_WE_scCODA*.ipynb` |
| 1 | H | `01_preprocessing/spatial/` (see missing code below) |
| 1 | I | `10_spatial_analysis/fig1I_marker_heatmap/` |
| 2 | A, B, D, F, H | `04_regulatory_annotation/Liver_analysis.ipynb` |
| 2 | C | `04_regulatory_annotation/Liver_analysis.ipynb` (SCENIC+ links, APA plots) |
| 2 | E | `06_genetic_enrichment/241221_WE_LDSC_Plotting.ipynb`, `250408_WE_FINRICH_and_Homer_For_TSCC_High_Fibrosis_and_cRE_States.ipynb`, `08_MASLD_loci/250124_WE_Hep_cREs_Intersect.ipynb` |
| 2 | G | `04_regulatory_annotation/241217_WE_Prep_HOMER_For_TSCC.ipynb`, `241224_WE_Homer_fGSEA_on_specific_links.ipynb` |
| 3 | A | `05_differential_analysis/*DEseq*.ipynb` |
| 3 | B | `05_differential_analysis/240829_WE_Liver_RNA_DEseq.ipynb` |
| 3 | D, E | `04_regulatory_annotation/Liver_analysis.ipynb`, `05_differential_analysis/241227_WE_Homer_DAC_cluster_sites.ipynb`, `06_genetic_enrichment/250219_*` |
| 3 | I | `06_genetic_enrichment/250220_WE_Prep_FINRICH_For_TSCC_GRNs.ipynb`, `241221_WE_Prep_FINRICH_For_TSCC_states_GRNs_H3K27me3_bins.ipynb` |
| 4 | A–E | Wet-lab experiments (qRT-PCR, ImageJ quantification); no code |
| 5 | A | `07_QTL/241122_WE_Plotting_Summaries.ipynb` |
| 5 | B | `07_QTL/241206_WE_MotifbreakR.ipynb` |
| 5 | C | `07_QTL/240930_WE_Coloc_QTLs.ipynb`, `07_QTL/241122_WE_Specificity_*QTLs.ipynb` |
| 5 | D | `08_MASLD_loci/240930_WE_Coloc_GWAS.ipynb` |
| 5 | E | `08_MASLD_loci/250325_WE_annotating_gwas_v2.ipynb` |
| 5 | F–I | `08_MASLD_loci/250303_WE_Plot_QTLs.ipynb`, `08_MASLD_loci/chrombpnet.variant.effects.sh` |
| 6 | A | `01_preprocessing/multiome/241018_WE_Liver_Peaks_Add_New_Peak_Mat_Our_Pipeline.ipynb` |
| 6 | C | `01_preprocessing/spatial/03_tangram/` |
| 6 | D | `03_composition/` (scCODA hepatocytes, miloR) |
| 6 | E | `04_regulatory_annotation/Liver_analysis.ipynb`, `241217_WE_Prep_HOMER_For_TSCC.ipynb` |
| 6 | G, H | `05_differential_analysis/*Subtypes*.ipynb` |
| 6 | I | `06_genetic_enrichment/250219_*`, `250408_*` |
| 6 | J | `09_hepatocyte_trajectory/241113_WE_Liver_Hepatocyte_Union_Peaks.ipynb` |
| 6 | K, L | `10_spatial_analysis/fig6K-L_inflammatory_hepatocytes/` |

Panels not listed (1A, 1D, 1E, 2C loops, 3C, 3F–H, 6B, 6F, 6M) are either schematics, genome browser views, or produced by code listed in [`MISSING_CODE.md`](MISSING_CODE.md).

## Data
Processed data are available at https://epigenome.wustl.edu/MASLD/. Raw and supplementary data will be available soon.

## Notes on running the code
The notebooks were run on the UCSD TSCC HPC and contain absolute paths to lab storage (e.g. `/tscc/projects/ps-gaultonlab/...`, `/nfs/lab/...`). To re-run, replace these with the locations of the processed data downloaded above. Notebook outputs are kept so results can be inspected without re-running. Notebooks larger than ~10 MB do not render on GitHub; download them or view them via [nbviewer](https://nbviewer.org/).

## Citation
Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis. *medRxiv* (2025). https://doi.org/10.1101/2025.05.09.25327043

## License
MIT, see [LICENSE](LICENSE).
