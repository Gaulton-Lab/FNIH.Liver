# Multiome Processing Notebooks

This directory contains Jupyter notebooks and scripts used for processing and analyzing multiome data from cellranger to clustering and cell annotation in the FNIH Liver project. Below is an overview of the files in this directory in rough order:

## File Descriptions

### `Genotype_Processing.md`
**Description**: Notes on imputation and preparing genotypes for demultiplexing

### `multiome_sample_processing_v6_Sierra.sh`
**Description**: Wrapper script to run initial barcode filtering, quantifying ATAC data in genomic windows, demuxlet for demultiplexing, and clustering.

### `1sample_metrics_filtering_v4.R`
**Description**: Wrapped by multiome_sample_processing_v6_Sierra.sh. First subscript called in multiome_sample_processing_v6_Sierra.sh. This script filters barcodes by # genes and # fragments and makes an initial clustering.

### `2ATAC_processing_v4.py`
**Description**: Wrapped by multiome_sample_processing_v6_Sierra.sh. Create counts per barcode across 50kb windows tiling the genome.

### `Popscle_demuxletAllinOneScript_final_v1.sh`
**Description**: Runs on combined RNA and ATAC bams for the 10X multiome libraries using a reference vcf.

### `3sample_postfilter_SoupX_ATACwindows_demuxlet_v5_fNIHmod.R`
**Description**: Wrapped by multiome_sample_processing_v6_Sierra.sh. This code runs ambient RNA removal, adds the windows based matrix, and filters barcodes based on demuxlet and adds donor assignments.

### `240516_Corban_fNIH_liver_1.MergingAllPoolsTogether.ipynb`
**Description**: Merge data from individual libraries. Plot additional QC metrics like TSS enrichment and filter. Add in meta data.

### `240620_WE_Liver_Windows_Subclustering_Doublets.ipynb`
**Description**: Cluster at high resolution to examine doublet enriched clusters and remove them.

### `240704_WE_Liver_Windows_Subclustering_CellTypes_Peak_Call.ipynb`
**Description**: Explore cell subtypes and call peaks

### `240812_WE_Liver_Peaks_Clustering.ipynb`
**Description**: Generate counts on peaks and recluster. 

### `240827_WE_Liver_Peaks_Add_New_Peak_Mat.ipynb`
**Description**: Explored cell subtypes and called peaks with Seurat. After this decided to call peaks via a different pipeline so this is redone almost completely in 241018_WE_Liver_Peaks_Add_New_Peak_Mat_Our_Pipeline.ipynb

### `241018_WE_Liver_Peaks_Add_New_Peak_Mat_Our_Pipeline.ipynb`
**Description**: Added in the peaks based matrix using the updated peak calling process and annotated cell types and subtypes as well as some core plots of the annotations. Peaks called based on this repo: https://github.com/Gaulton-Lab/peak-call-pipeline

### `241023_WE_Make_Barcode_List_Subtypes_and_Disease.ipynb`
**Description**: Wrote barcode lists for subtype and disease groups to make bigwigs and call subtype peaks.

### `241113_WE_Liver_Hepatocyte_Union_Peaks.ipynb`
**Description**: Add in a hepatocyte specific peak matrix, merge with the Hepatocyte paired tag, and explore the trajectory analysis.