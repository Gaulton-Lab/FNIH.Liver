#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
import datetime
import bin2cell as b2c

#from bin2cell import read_visium

def main():
    sc.logging.print_header()

    # Read arguments from the command line
    if len(sys.argv) != 9:
        #print("Usage: bin2cell.py two_micron_path source_image_path spatial_dir outdir patient_id mpp")
        print("Usage: bin2cell.py patient_id library_id condition source_img_path mic2_path spatial_dir_path outdir_samp_path mpp")
        # python $Script $PATIENT_ID $LIBRARY_ID $CONDITION $SOURCE_IMG_PATH $MIC2_PATH $SPATIAL_DIR $OUTDIR_SAMP $MPP
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

    odir_pat = outdir + patient_id 
    print(odir_pat)
    #create directory for stardist input/output files
    os.makedirs(odir_pat, exist_ok=True)

    adata = b2c.read_visium(two_micro_p, 
                        source_image_path = img_path, 
                        spaceranger_image_path = spatial_dir
                       )
    adata.var_names_make_unique()
    adata

    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)
    adata

    ### Add sample metadata
    adata.obs['patient_id'] = patient_id
    adata.obs['library_id'] = library_id
    adata.obs['condition'] = condition

    print("Saving scaled image")
    mpp = float(mpp_inp)
    img_key_str = str(mpp) + "_mpp"
    b2c.scaled_he_image(adata, mpp=mpp, save_path=odir_pat+"/"+patient_id+"_he.tiff")

    b2c.destripe(adata)

    #define a mask to easily pull out this region of the object in the future
    mask = ((adata.obs['array_row'] >= 1450) & 
        (adata.obs['array_row'] <= 1550) & 
        (adata.obs['array_col'] >= 250) & 
        (adata.obs['array_col'] <= 450)
        )

    bdata = adata[mask]
    # Saving image with mask
    sc.pl.spatial(bdata, color=[None, "n_counts", "n_counts_adjusted"], img_key=img_key_str,
              basis="spatial_cropped", show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_mask_ncounts.pdf', format = "pdf")

   
    print("Running segmentation on HE")
    current_time = datetime.datetime.now()
    print(current_time)
    b2c.stardist(image_path=odir_pat+"/"+patient_id+"_he.tiff", 
             labels_npz_path=odir_pat+"/"+patient_id+"_he.npz", 
             stardist_model="2D_versatile_he", 
             prob_thresh=0.01
            )
    current_time = datetime.datetime.now()
    print(current_time)



    b2c.insert_labels(adata, 
            labels_npz_path=odir_pat+"/"+patient_id+"_he.npz", 
            basis="spatial", 
            spatial_key="spatial_cropped",
            mpp=mpp, 
            labels_key="labels_he"
            )
  
    # Saving image with he labels
    bdata = adata[mask]
    bdata = bdata[bdata.obs['labels_he']>0]
    bdata.obs['labels_he'] = bdata.obs['labels_he'].astype(str)

    labels_he_fn = odir_pat + "/" +  patient_id + '_labels_he.pdf'

    sc.pl.spatial(bdata, color=[None, "labels_he"], img_key=img_key_str, basis="spatial_cropped",
             show = False)
    plt.savefig(labels_he_fn, format = "pdf")
 
    print("Saving HE render")
    #the label viewer wants a crop of the processed image
    #get the corresponding coordinates spanning the subset object
    crop = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial_cropped", mpp=mpp)

    #if this errors about missing get_cmap(), downgrade matplotlib below 3.9.0
    rendered = b2c.view_stardist_labels(image_path=odir_pat+"/"+patient_id+"_he.tiff", 
                                    labels_npz_path=odir_pat+"/"+patient_id+"_he.npz", 
                                    crop=crop
                                   )
    plt.imshow(rendered)
    plt.savefig(odir_pat+"/"+ patient_id+'_rendered_he.pdf')

    b2c.expand_labels(adata, 
                  labels_key='labels_he', 
                  expanded_labels_key="labels_he_expanded"
                 )


    print("Saving expanded labels")
    bdata = adata[mask]

    bdata = bdata[bdata.obs['labels_he_expanded']>0]
    bdata.obs['labels_he_expanded'] = bdata.obs['labels_he_expanded'].astype(str)

    sc.pl.spatial(bdata, color=[None, "labels_he_expanded"], img_key=img_key_str,
              basis="spatial_cropped", show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_he_expanded.pdf', format = "pdf")

    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, 
               save_path=odir_pat+"/"+patient_id+"_gex.tiff")

    print("Saving adata with HE")
    # Saving adata intermediate before gene expression version of segmentation
    adata.write_h5ad(odir_pat+"/"+patient_id+"_adata_he_only.h5ad")
    adata

    print("Running segmentation on gene expression")
    current_time = datetime.datetime.now()
    print(current_time)
    b2c.stardist(image_path=odir_pat+"/"+patient_id+"_gex.tiff", 
             labels_npz_path=odir_pat+"/"+patient_id+"_gex.npz", 
             stardist_model="2D_versatile_fluo", 
             prob_thresh=0.05, 
             nms_thresh=0.5
            )
    current_time = datetime.datetime.now()
    print(current_time)    


    b2c.insert_labels(adata,
                  labels_npz_path=odir_pat+"/"+patient_id+"_gex.npz",
                  basis="array",
                  mpp=mpp,
                  labels_key="labels_gex"
                 )
    print("Saving labels GEX")
    bdata = adata[mask]

    #0 means unassigned
    bdata = bdata[bdata.obs['labels_gex']>0]
    bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(str)

    sc.pl.spatial(bdata, color=[None, "labels_gex"], img_key=img_key_str, basis="spatial_cropped",
             show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_gex.pdf', format = "pdf")

   
    print("Saving rendered GEX")
    #the label viewer wants a crop of the processed image
    #get the corresponding coordinates spanning the subset object
    crop = b2c.get_crop(bdata, basis="array", mpp=mpp)

    #if this errors about missing get_cmap(), downgrade matplotlib below 3.9.0
    rendered = b2c.view_stardist_labels(image_path=odir_pat+"/"+patient_id+"_gex.tiff", 
                                    labels_npz_path=odir_pat+"/"+patient_id+"_gex.npz", 
                                    crop=crop
                                   )
    plt.imshow(rendered)
    plt.savefig(odir_pat+"/"+ patient_id+'_rendered_gex.pdf')


    b2c.salvage_secondary_labels(adata,
                             primary_label="labels_he_expanded",
                             secondary_label="labels_gex",
                             labels_key="labels_joint"
                            )

    print("Saving joint labels")
    bdata = adata[mask]

    #0 means unassigned
    bdata = bdata[bdata.obs['labels_joint']>0]
    bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(str)

    sc.pl.spatial(bdata, color=[None, "labels_joint_source", "labels_joint"], img_key=img_key_str,
              basis="spatial_cropped", show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_joint.pdf', format = "pdf")

    print("Saving bdata mask gex .h5ad")
    # Saving bdata masked after gene expression version of segmentation
    bdata.write_h5ad(odir_pat+"/"+patient_id+"_bdata_gex.h5ad")

    cdata = b2c.bin_to_cell(adata, labels_key="labels_joint", 
                        spatial_keys=["spatial", "spatial_cropped"])

    # Saving cdata whole object after gene expression version of segmentation
    cdata.write_h5ad(odir_pat+"/"+patient_id+"_cdata_b2c.h5ad")

    print("Saving joint labels on whole tissue")
    sc.pl.spatial(cdata, color=["bin_count", "labels_joint_source"], img_key=img_key_str, basis="spatial_cropped",
             show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_joint_bins_all.pdf', format = "pdf")

    print("Saving cdata whole object .h5ad")
    # Saving cdata whole object after gene expression version of segmentation
    cdata.write_h5ad(odir_pat+"/"+patient_id+"_cdata_b2c.h5ad")

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


    sc.pl.umap(cdata, color=["labels_joint_source", "clusters"], wspace=0.4, show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_joint_clusters_umap_all.pdf', format = "pdf")

    sc.pl.spatial(cdata, img_key="hires", color=["labels_joint_source", "clusters"], show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_joint_clusters_spatial_all.pdf', format = "pdf")

    print("Saving cdata whole object .h5ad after clustering")
    cdata.write_h5ad(odir_pat+"/"+patient_id+"_cdata_b2c_clust.h5ad")
    
    print("Plotting trimmed image of bin2cell called cells")
    cell_mask = ((cdata.obs['array_row'] >= 1450) & 
             (cdata.obs['array_row'] <= 1550) & 
             (cdata.obs['array_col'] >= 250) & 
             (cdata.obs['array_col'] <= 450)
            )

    ddata = cdata[cell_mask]
    sc.pl.spatial(ddata, color=["bin_count", "labels_joint_source"], img_key=img_key_str, basis="spatial_cropped")
    sc.pl.spatial(ddata, color=["bin_count", "labels_joint_source"], img_key=img_key_str, basis="spatial_cropped",
             show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_joint_bins.pdf', format = "pdf")


    sc.pl.umap(ddata, color=["labels_joint_source", "clusters"], wspace=0.4, show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_joint_clusters_umap.pdf', format = "pdf")
    
    sc.pl.spatial(ddata, img_key="hires", color=["labels_joint_source", "clusters"], show = False)
    plt.savefig(odir_pat+"/"+ patient_id+'_labels_joint_clusters_spatial.pdf', format = "pdf")

    print("Saving ddata trimmed object .h5ad after clustering")
    ddata.write_h5ad(odir_pat+"/"+patient_id+"_ddata_b2c_clust.h5ad")




if __name__ == "__main__":
    main()
