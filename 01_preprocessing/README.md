# 01 Preprocessing

Processing of each assay from raw data to annotated single-cell objects.

| Directory | Assay | Summary |
|---|---|---|
| [`genotypes/`](genotypes) | Illumina InfiniumCoreExome arrays | GenomeStudio export, QC with plink, McCarthy tools, TOPMed r3 imputation, VCF prep for demuxlet |
| [`multiome/`](multiome) | 10x Multiome (ATAC + RNA) | Per-library filtering, demuxlet, SoupX, merging, doublet removal, clustering, peak calling and cell type annotation |
| [`paired_tag/`](paired_tag) | Droplet Paired-Tag (H3K27ac, H3K27me3 + RNA) | cellranger-arc processing, QC/FRiP filtering, SoupX, clustering |
| [`droplet_hic/`](droplet_hic) | Droplet Hi-C | Barcode extraction, alignment, pairs processing, scHiCluster imputation, gene scores, compartments, domains, embedding |
| [`spatial/`](spatial) | Visium HD | SpaceRanger, bin2cell segmentation, Tangram label transfer |

## genotypes/
- `Genotype_Processing.md`: Notes and commands for array export, imputation on the TOPMed server, filtering (R2>0.9, MAF>0.01), and subsetting VCFs per pool for demuxlet.

## paired_tag/
- `07.proc_10Xarc_RNA.sh`, `08.proc_10Xarc_DNA.sh`: cellranger-arc/cellranger-atac processing of the RNA and histone (DNA) libraries.
- `DPT_help.R`: Helper functions for Droplet Paired-Tag analysis in R.
- `DPT_preprocess.ipynb`: QC of the DNA modality, FRiP, and manual DNA + RNA barcode filtering per library.
- `DPT_clustering.ipynb`: Reads demultiplexing results, merges libraries, applies quality filtering and runs SoupX per library.

Reference mapping of Paired-Tag onto the multiome is in [`02_integration/`](../02_integration).

## droplet_hic/
Run in numeric order:
- `00.split_single_cell.sh`: Split pairs into single-cell contact files.
- `01.filter-contact.sh`: Filter contacts (blacklist) with scHiCluster.
- `02.gene-score.sh`: scHiCluster gene scores (scGAD) at 10 kb for reference mapping.
- `03.preproc_paired_hic_tscc2.sh`, `04.proc_paired_hic_v4.sh`: Barcode extraction, trimming, BWA-MEM alignment and pairtools parsing/deduplication.
- `05.compartment.sh`: A/B compartment scores (100 kb).
- `06.insulation.sh`: Domain calling (25 kb).
- `hicluster_embedding.ipynb`: Builds the scGAD matrix and maps Hi-C cells onto the multiome RNA reference (PyNNDescent label transfer).

## spatial/
- `01_spaceranger/spaceranger_reruns_liver.sh`: `spaceranger count` with manual Loupe alignment.
- `02_bin2cell/`: StarDist/bin2cell segmentation of 2 µm bins into cells (`bin2cell_script.py`, array submission script and sample sheet).
- `03_tangram/`: Tangram mapping of hepatocyte sub-types onto Visium HD hepatocytes (Figure 6C).
  - `reference_markers/`: Notebooks defining the top-50 `rank_genes_groups` markers per disease group (normal, MASL, MASH).
  - `hepatocyte_subtypes_ds10k_ref/`, `hepatocyte_subtypes_highfib_lowfib/`: Split h5ad into 10k-cell patches, train Tangram on GPU, post-process and merge patches. The `.tsv` files list the input datasets.
