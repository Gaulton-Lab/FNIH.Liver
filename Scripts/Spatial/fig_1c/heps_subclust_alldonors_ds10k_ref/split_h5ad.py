import anndata as ad
import os
import scanpy as sc
import datetime
import matplotlib.pyplot as plt
import sys

### Inputs:
## input path to spatial .h5ad file to split
## input path output directory to create subsets


def main():
    sc.logging.print_header()

    # Directory containing the .h5ad files

    if len(sys.argv) != 3:
        print("Usage: python split_h5ad.py spatial_fp dir_to_tangram_outputs")
        sys.exit(1)

    else:
        print("Continuing!")
    
    spatial_fp = sys.argv[1]
    tan_pat_odir = sys.argv[2]

    ad_sp_full = sc.read_h5ad(spatial_fp)

    os.makedirs(tan_pat_odir, exist_ok = True)
    subsets_odir = tan_pat_odir + "subsets/"

    os.makedirs(subsets_odir, exist_ok = True)

    # Assuming you already have an AnnData object: adata
    n_cells_per_subset = 10000

    # Get the total number of cells
    n_cells = ad_sp_full.shape[0]

    # Calculate the number of subsets
    n_subsets = (n_cells + n_cells_per_subset - 1) // n_cells_per_subset  # ceil division

    # Loop over each subset, slice, and save it
    for i in range(n_subsets):
        start = i * n_cells_per_subset
        end = min((i + 1) * n_cells_per_subset, n_cells)  # Make sure we don't go out of bounds
        print(start)
        print(end)
        subset_adata = ad_sp_full[start:end, :]
       
        padded_number = f"{i+1:03d}" 
        # Define the filename for each subset
        subset_filename = f"subset_0{padded_number}_of_{n_subsets}.h5ad"

        out_fn =  subsets_odir + subset_filename
        out_fn
        
        # Save the subset to a new file
        subset_adata.write(out_fn)

        print(f"Saved: {out_fn}")
        

  

    print("Done saving files")


if __name__ == "__main__":
    main()
