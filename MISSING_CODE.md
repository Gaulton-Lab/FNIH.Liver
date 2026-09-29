# Missing code (to add before publication)

This list was built by comparing the Methods and figure legends against the code in this repository. Delete this file once each item is added or marked as not applicable.

## Missing analyses (no code in the repo)

| Methods section | What is missing | Figures | Likely owner | Updated |
|---|---|---|---|---|
| Visualization of ATAC and histone tracks | Scripts that split bams/tagAligns by cell type/condition and run `bamCoverage` (RPGC, no chrY) | 1D, 3F, 5F, 5H | Ren lab / Gaulton lab | Partial: `01_preprocessing/paired_tag/12.bamCoverage.sh`, `04_regulatory_annotation/chromhmm/extract10X_from_fragment.py` |
| Identifying chromatin loops | scHiCluster loop calling, Peakachu, Mustache, HiCCUPS; consensus loops (NMS / union-find); O/E loop strength (cooltools); eulerr plots | 1E, 2C, 3F, 5F, 5H | Ren lab (Yang Xie) | Partial: `04_regulatory_annotation/chromatin_loops/` (scHiCluster, Peakachu, HiCCUPS, Mustache, FitHiC2). Still missing: consensus loops, O/E loop strength, eulerr |
| Compartment analysis | Compartment switch / monotonic analysis and dcHiC differential compartments (only `hicluster compartment` is present) | Supp | Ren lab | No |
| APA | coolpuppy pileups and APA score calculation (only plotting of `APA.txt` is in `Liver_analysis.ipynb`) | 2C | Ren lab | Partial: `04_regulatory_annotation/chromatin_loops/pileup_celltype.py` |
| Loop annotations | Anchor annotation (TSS/cCRE), hepatocyte anchor k-means, clusterProfiler `enrichGO` | Supp | Ren lab | No |
| Genome annotation using ChromHMM | Splitting fragments by cell type, `BinarizeBed`, `LearnModel` (5 states). The notebook only reads the results | 2B | Ren lab | Yes: `04_regulatory_annotation/chromhmm/` |
| Evolutionary analysis of cCREs | RepeatMasker/TE age overlap, phastCons profiles (deepTools), gnomAD-SV CNV Fisher tests | Supp | – | Partial: `04_regulatory_annotation/evolutionary_analysis/` (gnomAD-SV). Still missing: RepeatMasker/TE, phastCons |
| Cell type marker genes | Entropy-based marker selection at the **cell type** level. The hepatocyte sub-type version is in `05_differential_analysis/240925_MK_Calculate_per_Hepatocyte_Cellsubtype_Entropy2_Clean.ipynb` (Fig 6B) | Supp | Gaulton lab | No |
| Cell type-specific cCREs | SnapATAC2 regression test (the notebook reads `snapatac2_cellsubtype_LR_test.csv`) | 2F, 6E | Ren lab | No |
| cis-regulatory modules | scikit-learn NMF runs and rank selection (the notebook reads the `NMF/res/*.r9n10` outputs) | 2D | Ren lab | Yes: `04_regulatory_annotation/nmf_modules/` (needs the external `snATACutils` helpers, see below) |
| Gene regulatory network modeling | SketchData downsampling, pycisTopic, SCENIC+ run | 2C, 3G–I, 6F | – | Partial: `04_regulatory_annotation/scenicplus/config.yaml` (SCENIC+ run config), `reports_otsu.ipynb` (eRegulon heatmaps). Still missing: SketchData downsampling, pycisTopic |
| GRN activity in MASLD | clusterProfiler `GSEA` on TF GRN targets; Cytoscape network export | 3G, 3H | – | Yes: `04_regulatory_annotation/scenicplus/Liver_GRNs_fgsea_plots.ipynb` (clusterProfiler GSEA and igraph network plot; confirm whether Cytoscape was still used for 3H) |
| GRN validation | regioneR `permTest` with ENCODE HepG2 ChIP-seq | Supp | – | No |
| GRN connectivity | igraph graph of TF GRNs linked by fine-mapped variants | 3I (bottom) | Gaulton lab | No |
| Effect size correlations | Spearman correlation across MASLD stages; cross-modality Pearson correlation | Supp | Gaulton lab | No |
| ChromBPNet | Bias and bias-factorized model training for the four hepatocyte conditions; variant scoring (`log_counts_diff`). Only the peak prep and allele bigwig script are present | 5E, 5G, 5I | – | No |
| Comparison to other studies | AUCell scoring of Gribben/Karpova markers; DESeq2 on Gribben et al.; fGSEA of GRNs on Gribben results | Supp | Gaulton lab | No |
| Pseudotime robustness | Palantir and Slingshot (SeuratExtend), per-donor analyses | Supp | – | No |
| Spatial: broad cell types | Tangram mapping of broad cell types and per-disease-group reference markers (only hepatocyte sub-type Tangram is present) | 1H | Spatial team | Yes: `01_preprocessing/spatial/03_tangram/broad_celltypes/`, `reference_markers/2411*`, `10_spatial_analysis/fig1H_broad_celltypes/` |
| Spatial: pathway scores | Module scores for fatty acid and ECM interaction genes | 3C | Spatial team | Yes: `10_spatial_analysis/supp_fig21_ecm_fatty_acid/` (KEGG ECM and Reactome fatty acid scores; confirm these also make Fig 3C) |
| Spatial: hepatocyte sub-type plots | Plotting of Tangram hepatocyte sub-type labels | 6C | Spatial team | No |
| Cell-cell interaction | CellChat on Visium HD | 6M | Spatial team | No |

## Helper scripts called but not included

These scripts are called by code in the repository but live elsewhere on TSCC. Add them to the repo or link a public repo that has them.

| Called from | Missing script(s) | Updated |
|---|---|---|
| `01_preprocessing/spatial/03_tangram/*/tans_tangram_train_sp_gpu_arr.sh` | `tans_tangram_train_sp_gpu_arr.py` | No. `Split_h5ad.py` was added, but the `.sh` wrappers still call the copies in `/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/` |
| `01_preprocessing/spatial/03_tangram/*/Tangram_postPrediction_arr.sh` | `Tangram_postPrediction_arr.py` | No |
| `01_preprocessing/spatial/03_tangram/*/merge_patches.sh` | `merge_patches.py` | No |
| `01_preprocessing/droplet_hic/*.sh` | `phc.count_pairs_sc.py`, `phc.plot_fragment.R`, `phc.summarize_pairs_lec.R`, `phc.batch_splitPairs_single.sh`, `scifi.CB_to_BB.py`, `10XcountFrag.py`, `getSize.py` | No |
| `04_regulatory_annotation/nmf_modules/*.sh` | `snATACutils` helpers: `nmfATAC.cluster2peak.lite.py`, `nmfATAC.stat.py`, `nmfATAC.plotH.R`, `nmfATAC.plotW.R`, `nmfATAC.statBox.R`, `nmfATAC.plotBox.R` (from `~/scripts/git/snATACutils/bin`); link the public repo if there is one | No |
| `04_regulatory_annotation/chromatin_loops/fithic2_run.sh` | `summarize_intra_contacts.py`, `mcool_to_fithic2_bias.py` | No |
| `07_QTL/`, `08_MASLD_loci/` | The HPC coloc and tensorQTL batch scripts; the notebooks say these were run on the HPC and the results copied back | No |

## Things to check

- `01_preprocessing/paired_tag/08.proc_10Xarc_DNA.sh` points to an **mm10** reference (`refdata-cellranger-arc-mm10-2020-A-2.0.0`). Check whether this is the version run on the human liver data.
- `04_regulatory_annotation/scenicplus/config.yaml` has the SCENIC+ inputs for **Hepatocytes**, but `dem_db_fname` points to the **Endothelial** cistarget scores database. Check whether this is intended.
- `04_regulatory_annotation/chromatin_loops/hicluster_make_sbatch_percondition_loops.sh` and `hicluster_setup_percondition_loops.sh` write paths for a "mediator" mount (`WRITE_BASE`/`PATH_BASE`). Check that they still run as documented.
- `01_preprocessing/multiome/240827_WE_Liver_Peaks_Add_New_Peak_Mat.ipynb` was superseded by `241018_...Our_Pipeline.ipynb`. Keep it for transparency or remove it.
- `03_merge_and_call_genotypes.sb` has the per-donor GATK `GenotypeConcordance` step commented out, and discordance was computed with `bcftools gtcheck` instead. Consider removing the commented block or noting this in the script.
- `06_genetic_enrichment/LDSC.sh` computes cALT–cirrhosis and cALT–NAFLD `--rg`. The Methods say correlations were calculated between all three GWAS, so NAFLD–cirrhosis may be missing. The script also uses FinnGen **R11** NAFLD summary stats, but the Methods cite FinnGen **R9**.
- The spatial SLURM scripts contain a personal email in `--mail-user`. Consider replacing it with a placeholder.
- The `Seurat5.0 DecontX` kernel (`seurat5.0.1.decontx`), used by most R notebooks, is not one of the exported conda environments. See `envs/README.md`.
