#!/usr/bin/env python3
"""
bin2cell_script.py

Per-sample Bin2cell processing pipeline for Visium HD spatial data.

Given a single sample's Space Ranger 2 micron output plus its H&E image,
this script runs the full Bin2cell (https://github.com/Teichlab/bin2cell)
workflow to go from raw 2 micron bins to single-cell-resolution "bin2cell"
(b2c) objects, with QC plots and intermediate `.h5ad` checkpoints saved
along the way. Broadly, the pipeline:

  1. Reads the raw Visium HD object and does light gene/bin filtering.
  2. Saves a scaled H&E image and destripes the counts.
  3. Segments cells from the H&E image with StarDist ("HE" segmentation),
     inserts those labels back onto the bins, and expands them.
  4. Builds a gene-expression-density image and segments *that* with
     StarDist as well ("GEX" segmentation), to catch cells the H&E
     segmentation misses.
  5. "Salvages" GEX-only labels into the expanded H&E labels to get a
     single joint label set, then aggregates bins into per-cell ("b2c")
     pseudo-cells with `b2c.bin_to_cell`.
  6. Normalizes, clusters (Leiden) and UMAPs the resulting b2c matrix,
     saving spatial/UMAP QC plots for both the whole tissue and a small
     cropped region (for a closer look at segmentation quality).

Throughout, `mask` / `cell_mask` define a small rectangular region of the
tissue (in bin/array coordinates) that gets re-plotted after every major
step as a quick visual sanity check of the segmentation, since the whole
image is too large/slow to inspect that way at every step.

Usage
-----
    python bin2cell_script.py \\
        <patient_id> <library_id> <condition> \\
        <source_img_path> <mic2_path> <spatial_dir_path> \\
        <outdir_samp_path> <mpp>

Positional arguments
---------------------
patient_id        Short sample/patient identifier; used as the prefix for
                   every output file and as the per-sample subdirectory
                   name under `outdir_samp_path`.
library_id        Sequencing library ID; stored as sample metadata only
                   (does not affect any file paths).
condition         Condition/group label (e.g. disease stage); stored as
                   sample metadata only.
source_img_path   Path to the full-resolution source H&E image.
mic2_path         Path to the Space Ranger 2 micron (`square_002um`)
                   output directory, passed to `b2c.read_visium`.
spatial_dir_path   Path to the matching Space Ranger `spatial` directory
                   (fiducials/scale factors etc.), passed to
                   `b2c.read_visium`.
outdir_samp_path   Base output directory; a `<patient_id>` subdirectory is
                   created under it for all of this sample's outputs.
mpp               Microns-per-pixel to render the H&E/GEX images at
                   (passed to Bin2cell's image-rendering/segmentation
                   steps); parsed as a float.
"""

import sys
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
import datetime
import bin2cell as b2c


def main():
    # Print scanpy/anndata/numpy version info, useful for debugging environment issues.
    sc.logging.print_header()

    # ---- Parse positional CLI arguments ----
    # Expected invocation (see module docstring for the full description of each):
    #   python bin2cell_script.py $PATIENT_ID $LIBRARY_ID $CONDITION \
    #       $SOURCE_IMG_PATH $MIC2_PATH $SPATIAL_DIR $OUTDIR_SAMP $MPP
    if len(sys.argv) != 9:
        print("Usage: bin2cell.py patient_id library_id condition source_img_path mic2_path spatial_dir_path outdir_samp_path mpp")
        sys.exit(1)
    else:
        print('Yay, all good in the hood!')

    patient_id = sys.argv[1]
    library_id = sys.argv[2]
    condition = sys.argv[3]
    img_path = sys.argv[4]
    two_micro_p = sys.argv[5]
    spatial_dir = sys.argv[6]
    outdir = sys.argv[7]
    mpp_inp = sys.argv[8]

    print('Reading in data')
    print(patient_id)
    print(library_id)
    print(condition)
    print(img_path)
    print(two_micro_p)
    print(spatial_dir)
    print(outdir)
    print(mpp_inp)

    # Per-sample output directory: <outdir>/<patient_id>/
    odir_pat = outdir + patient_id
    print(odir_pat)
    # Create directory for StarDist input/output files and all other per-sample outputs.
    os.makedirs(odir_pat, exist_ok=True)

    # ---- Load the raw Visium HD (2 micron bin) object ----
    adata = b2c.read_visium(
        two_micro_p,
        source_image_path=img_path,
        spaceranger_image_path=spatial_dir,
    )
    adata.var_names_make_unique()
    # NOTE: bare `adata` here is a no-op when run as a script (it only has an
    # effect in a notebook/REPL, where it would display the object's repr).
    # Kept as-is; it doesn't change behavior.
    adata

    # Light filtering: drop genes detected in fewer than 3 bins, and bins with
    # zero counts.
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)
    adata  # NOTE: also a no-op in script context; see above.

    # ---- Attach sample-level metadata (does not affect processing) ----
    adata.obs['patient_id'] = patient_id
    adata.obs['library_id'] = library_id
    adata.obs['condition'] = condition

    # ---- Render a scaled H&E image at the requested resolution ----
    # This rendered image is what StarDist segments in the "HE" pass below.
    print("Saving scaled image")
    mpp = float(mpp_inp)
    img_key_str = str(mpp) + "_mpp"
    b2c.scaled_he_image(adata, mpp=mpp, save_path=odir_pat + "/" + patient_id + "_he.tiff")

    # Destripe the raw bin counts (Bin2cell's correction for the systematic
    # row/column banding artifacts seen in Visium HD 2 micron bins).
    b2c.destripe(adata)

    # Define a small rectangular crop (in bin array coordinates) used purely
    # for QC plotting throughout the rest of the pipeline -- re-plotting the
    # full tissue after every step would be slow and hard to visually inspect.
    mask = ((adata.obs['array_row'] >= 1450) &
            (adata.obs['array_row'] <= 1550) &
            (adata.obs['array_col'] >= 250) &
            (adata.obs['array_col'] <= 450)
            )

    bdata = adata[mask]
    # QC plot: raw vs. destriped counts within the cropped region.
    sc.pl.spatial(bdata, color=[None, "n_counts", "n_counts_adjusted"], img_key=img_key_str,
                  basis="spatial_cropped", show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_mask_ncounts.pdf', format="pdf")

    # ---- StarDist segmentation pass 1: on the H&E image ----
    print("Running segmentation on HE")
    current_time = datetime.datetime.now()
    print(current_time)
    b2c.stardist(
        image_path=odir_pat + "/" + patient_id + "_he.tiff",
        labels_npz_path=odir_pat + "/" + patient_id + "_he.npz",
        stardist_model="2D_versatile_he",
        prob_thresh=0.01,
    )
    current_time = datetime.datetime.now()
    print(current_time)

    # Map the H&E-derived StarDist labels back onto the bins.
    b2c.insert_labels(
        adata,
        labels_npz_path=odir_pat + "/" + patient_id + "_he.npz",
        basis="spatial",
        spatial_key="spatial_cropped",
        mpp=mpp,
        labels_key="labels_he",
    )

    # QC plot: H&E-derived cell labels within the cropped region (label 0 = unassigned bin, dropped).
    bdata = adata[mask]
    bdata = bdata[bdata.obs['labels_he'] > 0]
    bdata.obs['labels_he'] = bdata.obs['labels_he'].astype(str)

    labels_he_fn = odir_pat + "/" + patient_id + '_labels_he.pdf'

    sc.pl.spatial(bdata, color=[None, "labels_he"], img_key=img_key_str, basis="spatial_cropped",
                  show=False)
    plt.savefig(labels_he_fn, format="pdf")

    # QC render: overlay the raw StarDist label mask on the cropped H&E image
    # itself (as opposed to the spatial scatter plots above), for a more
    # direct look at segmentation quality.
    print("Saving HE render")
    # The label viewer wants a crop of the processed image; get the
    # corresponding coordinates spanning the subset object.
    crop = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial_cropped", mpp=mpp)

    # If this errors about missing get_cmap(), downgrade matplotlib below 3.9.0.
    rendered = b2c.view_stardist_labels(
        image_path=odir_pat + "/" + patient_id + "_he.tiff",
        labels_npz_path=odir_pat + "/" + patient_id + "_he.npz",
        crop=crop,
    )
    plt.imshow(rendered)
    plt.savefig(odir_pat + "/" + patient_id + '_rendered_he.pdf')

    # Dilate the H&E-derived labels outward slightly so they cover bins just
    # outside the detected nucleus/cell boundary (nuclei-based segmentation
    # tends to undercall the true cell footprint).
    b2c.expand_labels(
        adata,
        labels_key='labels_he',
        expanded_labels_key="labels_he_expanded",
    )

    print("Saving expanded labels")
    bdata = adata[mask]

    bdata = bdata[bdata.obs['labels_he_expanded'] > 0]
    bdata.obs['labels_he_expanded'] = bdata.obs['labels_he_expanded'].astype(str)

    # QC plot: expanded H&E labels within the cropped region.
    sc.pl.spatial(bdata, color=[None, "labels_he_expanded"], img_key=img_key_str,
                  basis="spatial_cropped", show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_he_expanded.pdf', format="pdf")

    # Build a smoothed gene-expression-density image (used as the input for
    # the second, GEX-based StarDist segmentation pass below). This picks up
    # cells that are hard to see in the H&E image but still show a distinct
    # local expression signal.
    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5,
                    save_path=odir_pat + "/" + patient_id + "_gex.tiff")

    # Checkpoint: save the object as it stands after H&E-only segmentation,
    # before running the (slower) GEX segmentation pass.
    print("Saving adata with HE")
    adata.write_h5ad(odir_pat + "/" + patient_id + "_adata_he_only.h5ad")
    adata  # NOTE: no-op in script context; see note above.

    # ---- StarDist segmentation pass 2: on the gene-expression-density image ----
    print("Running segmentation on gene expression")
    current_time = datetime.datetime.now()
    print(current_time)
    b2c.stardist(
        image_path=odir_pat + "/" + patient_id + "_gex.tiff",
        labels_npz_path=odir_pat + "/" + patient_id + "_gex.npz",
        stardist_model="2D_versatile_fluo",
        prob_thresh=0.05,
        nms_thresh=0.5,
    )
    current_time = datetime.datetime.now()
    print(current_time)

    # Map the GEX-derived StarDist labels back onto the bins. Note this uses
    # basis="array" (bin grid coordinates) rather than "spatial", since the
    # GEX image was rendered directly from the bin grid.
    b2c.insert_labels(
        adata,
        labels_npz_path=odir_pat + "/" + patient_id + "_gex.npz",
        basis="array",
        mpp=mpp,
        labels_key="labels_gex",
    )
    print("Saving labels GEX")
    bdata = adata[mask]

    # 0 means unassigned.
    bdata = bdata[bdata.obs['labels_gex'] > 0]
    bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(str)

    # QC plot: GEX-derived cell labels within the cropped region.
    sc.pl.spatial(bdata, color=[None, "labels_gex"], img_key=img_key_str, basis="spatial_cropped",
                  show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_gex.pdf', format="pdf")

    print("Saving rendered GEX")
    # The label viewer wants a crop of the processed image; get the
    # corresponding coordinates spanning the subset object.
    crop = b2c.get_crop(bdata, basis="array", mpp=mpp)

    # If this errors about missing get_cmap(), downgrade matplotlib below 3.9.0.
    rendered = b2c.view_stardist_labels(
        image_path=odir_pat + "/" + patient_id + "_gex.tiff",
        labels_npz_path=odir_pat + "/" + patient_id + "_gex.npz",
        crop=crop,
    )
    plt.imshow(rendered)
    plt.savefig(odir_pat + "/" + patient_id + '_rendered_gex.pdf')

    # ---- Combine the two segmentations into one joint label set ----
    # Start from the (dilated) H&E labels as the primary source of truth, and
    # "salvage" in any additional cells found only by the GEX segmentation
    # (i.e. bins not covered by any H&E label) using the GEX labels instead.
    b2c.salvage_secondary_labels(
        adata,
        primary_label="labels_he_expanded",
        secondary_label="labels_gex",
        labels_key="labels_joint",
    )

    print("Saving joint labels")
    bdata = adata[mask]

    # 0 means unassigned.
    bdata = bdata[bdata.obs['labels_joint'] > 0]
    bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(str)

    # QC plot: joint labels within the cropped region, plus `labels_joint_source`
    # (which segmentation -- HE or GEX -- each label ultimately came from).
    sc.pl.spatial(bdata, color=[None, "labels_joint_source", "labels_joint"], img_key=img_key_str,
                  basis="spatial_cropped", show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_joint.pdf', format="pdf")

    print("Saving bdata mask gex .h5ad")
    # Save the small cropped/masked object (post joint-label assignment) as
    # its own checkpoint, separate from the full-tissue objects saved below.
    bdata.write_h5ad(odir_pat + "/" + patient_id + "_bdata_gex.h5ad")

    # ---- Aggregate bins into per-cell ("b2c") pseudo-cells ----
    # Sums all bins sharing the same joint label into one pseudo-cell,
    # carrying along both the original and cropped spatial coordinates.
    cdata = b2c.bin_to_cell(adata, labels_key="labels_joint",
                             spatial_keys=["spatial", "spatial_cropped"])

    # Checkpoint: whole-tissue b2c object before normalization/clustering.
    cdata.write_h5ad(odir_pat + "/" + patient_id + "_cdata_b2c.h5ad")

    print("Saving joint labels on whole tissue")
    sc.pl.spatial(cdata, color=["bin_count", "labels_joint_source"], img_key=img_key_str, basis="spatial_cropped",
                  show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_joint_bins_all.pdf', format="pdf")

    print("Saving cdata whole object .h5ad")
    # NOTE: this re-saves cdata to the exact same path as the write two lines
    # above `_cdata_b2c.h5ad`; cdata hasn't changed in between (only a plot
    # was made), so this second write is redundant. Left as-is per the
    # original script.
    cdata.write_h5ad(odir_pat + "/" + patient_id + "_cdata_b2c.h5ad")

    # ---- Normalize, find HVGs, and cluster the b2c matrix ----
    print("Normalizing and clustering bin2cell matrix")
    sc.pp.normalize_total(cdata, inplace=True)
    sc.pp.log1p(cdata)
    sc.pp.highly_variable_genes(cdata, flavor="seurat", n_top_genes=2000)

    sc.pp.pca(cdata)
    sc.pp.neighbors(cdata)
    sc.tl.umap(cdata)
    sc.tl.leiden(
        cdata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2
    )

    # QC plots: UMAP and whole-tissue spatial, colored by label source and cluster.
    sc.pl.umap(cdata, color=["labels_joint_source", "clusters"], wspace=0.4, show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_joint_clusters_umap_all.pdf', format="pdf")

    sc.pl.spatial(cdata, img_key="hires", color=["labels_joint_source", "clusters"], show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_joint_clusters_spatial_all.pdf', format="pdf")

    print("Saving cdata whole object .h5ad after clustering")
    cdata.write_h5ad(odir_pat + "/" + patient_id + "_cdata_b2c_clust.h5ad")

    # ---- Repeat the cropped-region QC plots, now on the clustered b2c object ----
    print("Plotting trimmed image of bin2cell called cells")
    cell_mask = ((cdata.obs['array_row'] >= 1450) &
                 (cdata.obs['array_row'] <= 1550) &
                 (cdata.obs['array_col'] >= 250) &
                 (cdata.obs['array_col'] <= 450)
                 )

    ddata = cdata[cell_mask]
    # NOTE: plotted twice in a row -- once interactively (show defaults to
    # True) and once with show=False purely to save the PDF. Both calls are
    # kept as in the original script (the first is effectively a no-op when
    # run non-interactively, e.g. via `python bin2cell_script.py ...`).
    sc.pl.spatial(ddata, color=["bin_count", "labels_joint_source"], img_key=img_key_str, basis="spatial_cropped")
    sc.pl.spatial(ddata, color=["bin_count", "labels_joint_source"], img_key=img_key_str, basis="spatial_cropped",
                  show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_joint_bins.pdf', format="pdf")

    sc.pl.umap(ddata, color=["labels_joint_source", "clusters"], wspace=0.4, show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_joint_clusters_umap.pdf', format="pdf")

    sc.pl.spatial(ddata, img_key="hires", color=["labels_joint_source", "clusters"], show=False)
    plt.savefig(odir_pat + "/" + patient_id + '_labels_joint_clusters_spatial.pdf', format="pdf")

    print("Saving ddata trimmed object .h5ad after clustering")
    ddata.write_h5ad(odir_pat + "/" + patient_id + "_ddata_b2c_clust.h5ad")


if __name__ == "__main__":
    main()
