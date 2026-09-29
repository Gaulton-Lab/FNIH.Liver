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

## chromhmm/
- `chmm.md`: ChromHMM commands for splitting fragments by cell type, `BinarizeBed`, `LearnModel` (5–8 states tested, 5 used), `MakeSegmentation`, and annotating cCREs with states (Figure 2B).
- `extract10X_from_fragment.py`: Splits 10x fragment files into per-cell type/mark bed files for ChromHMM.

## nmf_modules/
NMF of the cell type-by-cCRE accessibility matrix into cis-regulatory modules (Figure 2D). Run in order: `proc_npz.sh` (build the input matrix), `runNMF.sh` (100 runs per rank), `stable_rank_select.sh` (sparseness/entropy rank selection). These scripts call the `nmfATAC.*` helpers from `snATACutils` (not included).

## chromatin_loops/
Droplet Hi-C loop calling at 10 kb (Figures 1E, 2C, 3F, 5F, 5H):
- `hicluster_setup_percondition_loops.sh`, `hicluster_make_sbatch_percondition_loops.sh`: scHiCluster loop calling per cell type and condition.
- `peakachu_run.sh`, `mustache_run.sh`, `hiccups_run.sh`: Peakachu, Mustache (p < 0.1), and HiCCUPS (juicer_tools 1.22.01, CPU mode).
- `fithic2_run.sh`: Two-stage hepatocyte loop filter (Peakachu t=0.7, then FitHiC2 FDR < 0.01).
- `pileup_celltype.py`: coolpuppy observed/expected pileups (APA) of cell type loop sets across cell types.

## evolutionary_analysis/
- `extract_cnv.sh`, `run_gnomad_sv_overlap.sh`: Extract PASS gnomAD-SV v4.1 CNVs (DEL/DUP/CNV), stratify them by allele frequency, and test cCRE overlap (Supplementary).

## scenicplus/
SCENIC+ gene regulatory networks (Figures 2C, 3G–H, 6F):
- `config.yaml`: SCENIC+ Snakemake configuration for the hepatocyte run.
- `reports_otsu.ipynb`: eRegulon reports, specificity scores, and heatmaps using Otsu-thresholded region sets.
- `heatmap_utils.py`: Plotting helpers used by `reports_otsu.ipynb`.
- `Liver_GRNs_fgsea_plots.ipynb`: clusterProfiler GSEA of TF GRN target genes across MASLD stages, NES heatmaps, and GRN network plots (Figures 3G and 3H).
- `QTL_Overlap.ipynb`: Overlap of GRN regions and target genes with QTLs and colocalized signals.
