#!/usr/bin/env python3
"""
Split_h5ad.py

Split one large spatial `.h5ad` file into fixed-size chunks so each chunk
can be run through Tangram independently (in a SLURM array job).

Each chunk contains at most `n_cells_per_subset` cells/bins (in original
order, i.e. simple contiguous slicing -- not random), and is written to
`<tan_pat_odir>/subsets/subset_0<NNN>_of_<n_subsets>.h5ad`.

Usage
-----
    python Split_h5ad.py <spatial_fp> <dir_to_tangram_outputs>

Arguments
---------
spatial_fp             Path to the spatial `.h5ad` file to split.
dir_to_tangram_outputs  Output directory for this sample's Tangram run; a
                        `subsets/` subdirectory is created under it to hold
                        the split files.
"""

import os
import sys

import scanpy as sc

# Number of cells/bins per split file.
N_CELLS_PER_SUBSET = 10000


def main():
    sc.logging.print_header()

    if len(sys.argv) != 3:
        print("Usage: python Split_h5ad.py spatial_fp dir_to_tangram_outputs")
        sys.exit(1)
    else:
        print("Continuing!")

    spatial_fp = sys.argv[1]
    tan_pat_odir = sys.argv[2]

    ad_sp_full = sc.read_h5ad(spatial_fp)

    os.makedirs(tan_pat_odir, exist_ok=True)
    subsets_odir = tan_pat_odir + "subsets/"
    os.makedirs(subsets_odir, exist_ok=True)

    n_cells = ad_sp_full.shape[0]
    n_subsets = (n_cells + N_CELLS_PER_SUBSET - 1) // N_CELLS_PER_SUBSET  # ceil division

    for i in range(n_subsets):
        start = i * N_CELLS_PER_SUBSET
        end = min((i + 1) * N_CELLS_PER_SUBSET, n_cells)
        print(start)
        print(end)
        subset_adata = ad_sp_full[start:end, :]

        # NOTE: the leading "0" below is literal, in addition to the 3-digit
        # zero-padding from the f-string -- e.g. the first subset is named
        # "subset_0001_of_N.h5ad", not "subset_001_of_N.h5ad". Kept as-is
        # since downstream scripts (e.g. tans_tangram_train_sp_gpu_arr.sh)
        # just glob "*.h5ad" and don't depend on a specific digit count.
        padded_number = f"{i + 1:03d}"
        subset_filename = f"subset_0{padded_number}_of_{n_subsets}.h5ad"
        out_fn = subsets_odir + subset_filename

        subset_adata.write(out_fn)
        print(f"Saved: {out_fn}")

    print("Done saving files")


if __name__ == "__main__":
    main()
