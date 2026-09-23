# Missing code (to add before publication)

This list was built by comparing the Methods and figure legends against the code in this repository. Delete this file once each item is added or marked as not applicable.

## Missing analyses (no code in the repo)

| Methods section | What is missing | Figures | Likely owner |
|---|---|---|---|
| Single cell multiome clustering | Demultiplexing benchmark: sinto split by donor, GATK HaplotypeCaller, `bcftools gtcheck` discordance | Supp | Gaulton lab |
| Visualization of ATAC and histone tracks | Scripts that split bams/tagAligns by cell type/condition and run `bamCoverage` (RPGC, no chrY) | 1D, 3F, 5F, 5H | Ren lab / Gaulton lab |
| Identifying chromatin loops | scHiCluster loop calling, Peakachu, Mustache, HiCCUPS; consensus loops (NMS / union-find); O/E loop strength (cooltools); eulerr plots | 1E, 2C, 3F, 5F, 5H | Ren lab (Yang Xie) |
| Compartment analysis | Compartment switch / monotonic analysis and dcHiC differential compartments (only `hicluster compartment` is present) | Supp | Ren lab |
| APA | coolpuppy pileups and APA score calculation (only plotting of `APA.txt` is in `Liver_analysis.ipynb`) | 2C | Ren lab |
| Loop annotations | Anchor annotation (TSS/cCRE), hepatocyte anchor k-means, clusterProfiler `enrichGO` | Supp | Ren lab |
| Genome annotation using ChromHMM | Splitting fragments by cell type, `BinarizeBed`, `LearnModel` (5 states). The notebook only reads the results | 2B | Ren lab |
| Evolutionary analysis of cCREs | RepeatMasker/TE age overlap, phastCons profiles (deepTools), gnomAD-SV CNV Fisher tests | Supp | – |
| Cell type marker genes | Entropy-based marker selection (DescTools `Entropy`) | 6B, Supp | Gaulton lab |
| Cell type-specific cCREs | SnapATAC2 regression test (the notebook reads `snapatac2_cellsubtype_LR_test.csv`) | 2F, 6E | Ren lab |
| cis-regulatory modules | scikit-learn NMF runs and rank selection (the notebook reads the `NMF/res/*.r9n10` outputs) | 2D | Ren lab |
| Gene regulatory network modeling | SketchData downsampling, pycisTopic, SCENIC+ run | 2C, 3G–I, 6F | – |
| GRN activity in MASLD | clusterProfiler `GSEA` on TF GRN targets; Cytoscape network export | 3G, 3H | – |
| GRN validation | regioneR `permTest` with ENCODE HepG2 ChIP-seq | Supp | – |
| GRN connectivity | igraph graph of TF GRNs linked by fine-mapped variants | 3I (bottom) | Gaulton lab |
| Effect size correlations | Spearman correlation across MASLD stages; cross-modality Pearson correlation | Supp | Gaulton lab |
| LDSC | Genetic correlation (`--rg`) between cALT, NAFLD and cirrhosis GWAS | Supp | Gaulton lab |
| ChromBPNet | Bias and bias-factorized model training for the four hepatocyte conditions; variant scoring (`log_counts_diff`). Only the peak prep and allele bigwig script are present | 5E, 5G, 5I | – |
| Comparison to other studies | AUCell scoring of Gribben/Karpova markers; DESeq2 on Gribben et al.; fGSEA of GRNs on Gribben results | Supp | Gaulton lab |
| Pseudotime robustness | Palantir and Slingshot (SeuratExtend), per-donor analyses | Supp | – |
| Spatial: broad cell types | Tangram mapping of broad cell types and per-disease-group reference markers (only hepatocyte sub-type Tangram is present) | 1H | Spatial team |
| Spatial: pathway scores | Module scores for fatty acid and ECM interaction genes | 3C | Spatial team |
| Spatial: hepatocyte sub-type plots | Plotting of Tangram hepatocyte sub-type labels | 6C | Spatial team |
| Cell-cell interaction | CellChat on Visium HD | 6M | Spatial team |

## Helper scripts called but not included

These scripts are called by code in the repository but live elsewhere on TSCC. Add them to the repo or link a public repo that has them.

| Called from | Missing script(s) |
|---|---|
| `01_preprocessing/spatial/03_tangram/*/tans_tangram_train_sp_gpu_arr.sh` | `tans_tangram_train_sp_gpu_arr.py` |
| `01_preprocessing/spatial/03_tangram/*/Tangram_postPrediction_arr.sh` | `Tangram_postPrediction_arr.py` |
| `01_preprocessing/spatial/03_tangram/*/merge_patches.sh` | `merge_patches.py` |
| `01_preprocessing/droplet_hic/*.sh` | `phc.count_pairs_sc.py`, `phc.plot_fragment.R`, `phc.summarize_pairs_lec.R`, `phc.batch_splitPairs_single.sh`, `scifi.CB_to_BB.py`, `10XcountFrag.py`, `getSize.py` |
| `07_QTL/`, `08_MASLD_loci/` | The HPC coloc and tensorQTL batch scripts; the notebooks say these were run on the HPC and the results copied back |

## Things to check

- `01_preprocessing/paired_tag/08.proc_10Xarc_DNA.sh` points to an **mm10** reference (`refdata-cellranger-arc-mm10-2020-A-2.0.0`). Check whether this is the version run on the human liver data.
- `01_preprocessing/multiome/240827_WE_Liver_Peaks_Add_New_Peak_Mat.ipynb` was superseded by `241018_...Our_Pipeline.ipynb`. Keep it for transparency or remove it.
- The spatial SLURM scripts contain a personal email in `--mail-user`. Consider replacing it with a placeholder.
- The `Seurat5.0 DecontX` kernel (`seurat5.0.1.decontx`), used by most R notebooks, is not one of the exported conda environments. See `envs/README.md`.
