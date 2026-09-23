# 02 Integration

### `Liver_integration.ipynb`
Maps Droplet Paired-Tag (SCTransform, `FindTransferAnchors`/`MapQuery`) and Droplet Hi-C (scGAD embedding from `01_preprocessing/droplet_hic/hicluster_embedding.ipynb`, `ProjectUMAP`) onto the multiome RNA reference. Filters cells by prediction score and writes the final object. Produces:
- joint UMAP of all modalities (Figure 1B)
- marker gene activity across modalities (Figure 1C)
- donor characteristics (Figure 1F)
- cell type proportions by disease and fibrosis (Figure 1G bars)
- QC by modality and label-transfer prediction scores (Supplementary)
