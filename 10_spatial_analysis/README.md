# 10 Spatial analysis

Downstream analysis of Visium HD data after segmentation and Tangram mapping (see [`01_preprocessing/spatial/`](../01_preprocessing/spatial)).

### `fig1I_marker_heatmap/`
- `250321_01_imputed_heatmap_all_downsample_using_nash_fibstage2_rm_thy1_aspn.ipynb`: Heatmap of Tangram-imputed marker gene expression per cell type, downsampled to 100k cells per dataset (Figure 1I).

### `fig6K-L_inflammatory_hepatocytes/`
Run in this order:
1. `250303_01_module_score_top50_deseq2_fc_hep_sen.ipynb`: Module score of the top 50 genes up-regulated in inflammatory hepatocytes in Fib+ MASH (from DESeq2) in the Fib+ MASH donor. Labels the top 1% as inflammatory+ hepatocytes.
2. `250303_03_neighborhood_analysis_pt02_cells_nash_donor_stage4_highfib_sen_label_subset_bylogfc.ipynb`: Squidpy `nhood_enrichment` for cell types, hepatocyte sub-types, and inflammatory+ hepatocytes (Figure 6L).
3. `250312_01_neighborhood_analysis_recolor.ipynb`: Final colouring of the neighborhood enrichment heatmaps.
4. `250401_spatial_plots_injured_pop_fibstage4_figure5.ipynb`: Spatial plots of inflammatory+ hepatocytes (Figure 6K).

### `supp_fig32_all_donors/`
Parameterized versions of the above, looped over all Visium HD donors (Supplementary Figure 32):
- `260918_01_module_score_top50_deseq2_fc_hep_sen_mash_stage2_clean_parametrized.ipynb`
- `260918_02_neighborhood_analysis_combined_looped.ipynb`
