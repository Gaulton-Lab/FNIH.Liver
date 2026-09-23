# Environments

Conda environments exported with `conda env export --no-builds` from the environments used on TSCC. Create one with:

```bash
conda env create -f envs/<name>.yml
```

These exports reflect the environments as of 2026-09. Some packages may have been updated since the analyses were run. Versions reported in the Methods take precedence.

| Notebook kernel | Environment file | Used in |
|---|---|---|
| `Seurat5.0 DecontX` | `seurat5.yml` (closest match, see note) | Most R notebooks in 01, 05, 06, 07, 08, 09 |
| `seurat5` | `seurat5.yml` | `03_composition/260302_WE_scCODA_R_Plotting.ipynb` |
| `R 4.2` | `multiome_pipeline.yml` | `01_preprocessing/multiome/` pipeline and merge notebooks |
| `FINRICH` | `FINRICH.yml` (+ `homer.yml` for HOMER) | `04_regulatory_annotation/`, `06_genetic_enrichment/` |
| `coloc.susie` | `coloc.susie.yml` | `07_QTL/240930_WE_Coloc_QTLs.ipynb`, `08_MASLD_loci/` |
| `motifbreakR` | `motifbreakR.yml` | `07_QTL/241206_WE_MotifbreakR.ipynb` |
| `scCODA` | `scCODA.yml` | `03_composition/260302_WE_scCODA.ipynb` |
| `milor` | `milor.yml` | `03_composition/260203_miloR.ipynb` |
| (command line) | `a100.tensorQTL.yml` | tensorQTL runs from `07_QTL/240911_WE_TensorQTL_Inputs.ipynb` |
| (command line) | `ldsc.yml` | `06_genetic_enrichment/LDSC.sh` |

**Not yet exported:**
- The `Seurat5.0 DecontX` kernel (`seurat5.0.1.decontx`), which is not one of the exported conda environments.
- The `R 4.1` kernel (mashr notebooks in `07_QTL/`).
- The Ren lab environments (`seurat`, `schicluster`, and the R kernel used in `02_integration/` and `04_regulatory_annotation/`).
- The spatial Python environment (scanpy, squidpy, tangram, bin2cell).
- chromBPNet.
