# FNIH Liver  <img src="./images/FNIH_liver_color.png" align="right" width="200"/>
FNIH Liver is a single-nuclei, multi-modal epigenomic atlas of metabolic dysfunction-associated steatotic liver disease (MASLD).
<br/><br/>

This repository contains the code used for the manuscript: **Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis** ([medRxiv 10.1101/2025.05.09.25327043](https://www.medrxiv.org/content/10.1101/2025.05.09.25327043v1)).

## Abstract
Metabolic dysfunction-associated steatotic liver disease (MASLD) has limited treatments, and cell type-specific regulatory networks driving MASLD represent therapeutic avenues. We assayed five transcriptomic and epigenomic modalities in 2.4M cells from 86 livers across MASLD stages. Integrating modalities increased annotation of the genome in liver cell types several-fold over previous catalogs. We identified cell type regulatory networks of MASLD progression, including distinct hepatocyte networks driving MASL and mild and severe fibrosis MASH. Our single cell atlas annotated 88% of MASH-associated loci, including a third affecting hepatocyte regulation which we linked to distal target genes. Finally, we characterized hepatocyte heterogeneity, including MASH-enriched populations with altered repression, localization, and signaling. Overall, our results provide high-resolution maps of liver cell types and revealed novel targets for anti-MASH therapy.

## Repository layout

Directories are numbered in the order the analyses were run. Each directory has its own README describing every file.

| Directory | Contents |
|---|---|
| [`01_preprocessing/`](01_preprocessing) | Genotyping/imputation and demultiplexing validation, 10x multiome, Droplet Paired-Tag, Droplet Hi-C and Visium HD processing (incl. Tangram) |
| [`02_integration/`](02_integration) | Reference mapping of Paired-Tag and Hi-C onto multiome, joint UMAP, donor/QC summaries |
| [`03_composition/`](03_composition) | Cell type proportion changes (scCODA, miloR) |
| [`04_regulatory_annotation/`](04_regulatory_annotation) | cCREs, ChromHMM states, cell type-specific cCREs, NMF modules, motifs, SCENIC+ GRNs, Hi-C chromatin loops, gnomAD-SV overlap |
| [`05_differential_analysis/`](05_differential_analysis) | Pseudobulk, DESeq2 and fGSEA across MASLD stages for cell types and hepatocyte sub-types |
| [`06_genetic_enrichment/`](06_genetic_enrichment) | LDSC partitioned heritability, FINRICH fine-mapped variant enrichment, HOMER on differential/clustered sites |
| [`07_QTL/`](07_QTL) | tensorQTL mapping, mashr specificity, cross-modality QTL coloc, motifbreakR |
| [`08_MASLD_loci/`](08_MASLD_loci) | GWAS–QTL colocalization, annotation of MASLD loci, chromBPNet variant effects, locus plots |
| [`09_hepatocyte_trajectory/`](09_hepatocyte_trajectory) | Hepatocyte union peaks and Monocle3 pseudotime |
| [`10_spatial_analysis/`](10_spatial_analysis) | Visium HD cell type plots, marker heatmaps, pathway and inflammatory hepatocyte module scores, neighborhood enrichment |
| [`envs/`](envs) | Conda environments for the Jupyter kernels used in the notebooks |

Related repositories:
- Peak calling pipeline: https://github.com/Gaulton-Lab/peak-call-pipeline

## Data
Processed data are available at https://epigenome.wustl.edu/MASLD/. Raw and supplementary data will be available soon.

## Notes on running the code
The notebooks were run on the UCSD TSCC HPC and contain absolute paths to lab storage (e.g. `/tscc/projects/ps-gaultonlab/...`, `/nfs/lab/...`). To re-run, replace these with the locations of the processed data downloaded above. Notebook outputs are kept so results can be inspected without re-running. Notebooks larger than ~10 MB do not render on GitHub; download them or view them via [nbviewer](https://nbviewer.org/).

## Citation
Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis. *medRxiv* (2025). https://doi.org/10.1101/2025.05.09.25327043

## License
MIT, see [LICENSE](LICENSE).
