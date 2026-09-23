# Multiome processing

Jupyter notebooks and scripts for processing the 10x multiome data from Cell Ranger ARC output to clustering and cell type annotation. Files are listed in run order.

## pipeline/
Per-library processing, run via `multiome_sample_processing_v6_Sierra.sh`.

### `multiome_sample_processing_v6_Sierra.sh`
Wrapper that runs initial barcode filtering, quantifies ATAC in genomic windows, runs demuxlet, and clusters.

### `1sample_metrics_filtering_v4.R`
First sub-script. Filters barcodes by number of genes and fragments and makes an initial clustering.

### `2ATAC_processing_v4.py`
Second sub-script. Creates counts per barcode across windows tiling the genome.

### `Popscle_demuxletAllinOneScript_final_v1.sh`
Runs demuxlet on the merged RNA and ATAC bams for each library using the pool's genotype VCF.

### `3sample_postfilter_SoupX_ATACwindows_demuxlet_v5_fNIHmod.R`
Third sub-script. Runs ambient RNA removal, adds the windows matrix, filters barcodes using demuxlet, and adds donor assignments.

### `write_pool_3_4_commands.sh`, `parallel_pool_3_4_commands.sh`
Write and run the pipeline commands for the libraries in pools 3 and 4.

## Notebooks

### `240510_Corban_fNIH_liver_Using_rScript2output_demuxlet_merge.ipynb`
Merges the first four libraries, with batch correction and clustering checks.

### `240516_Corban_fNIH_liver_1.MergingAllPoolsTogether.ipynb`
Merges data from all libraries, plots additional QC metrics such as TSS enrichment, filters cells, and adds metadata.

### `240620_WE_Liver_Windows_Subclustering_Doublets.ipynb`
Clusters at high resolution to find doublet-enriched clusters (Amulet) and remove them.

### `240704_WE_Liver_Windows_Subclustering_CellTypes_Peak_Call.ipynb`
Explores cell sub-types and calls initial peaks per cell type.

### `240812_WE_Liver_Peaks_Clustering.ipynb`
Generates counts on peaks and re-clusters.

### `240827_WE_Liver_Peaks_Add_New_Peak_Mat.ipynb`
Explores cell sub-types with peaks called in Seurat. Superseded by `241018_WE_Liver_Peaks_Add_New_Peak_Mat_Our_Pipeline.ipynb`, which redoes almost all of it.

### `241018_WE_Liver_Peaks_Add_New_Peak_Mat_Our_Pipeline.ipynb`
Adds the peak matrix from the updated peak calling ([peak-call-pipeline](https://github.com/Gaulton-Lab/peak-call-pipeline)), annotates cell types and sub-types, and makes core annotation plots (Figure 6A).

### `241023_WE_Make_Barcode_List_Subtypes_and_Disease.ipynb`
Writes barcode lists for sub-type and disease groups, used to make bigwigs and call sub-type peaks.
