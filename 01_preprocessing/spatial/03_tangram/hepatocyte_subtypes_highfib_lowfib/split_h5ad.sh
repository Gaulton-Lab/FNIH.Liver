#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=split_h5ad
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Heart/spatial/split_h5ad.sh.e
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Heart/spatial/split_h5ad.sh.o
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

# SLURM wrapper for Split_h5ad.py: splits a spatial .h5ad into fixed-size
# chunks (subsets/) ahead of a Tangram array run.
#
# Expected environment (set by the caller via `sbatch --export=...`):
#   SPA_FP   Path to the spatial .h5ad file to split.
#   OUT_DIR  Per-sample output directory; a subsets/ subdirectory is
#            created under it by Split_h5ad.py.

source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh
conda activate spatial_backup

Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/split_h5ad.py'

## Check paths before executing
echo "split_h5ad FP: $SPA_FP"
echo "split_h5ad OUT_DIR: $OUT_DIR"

## Usage: python Split_h5ad.py spatial_fp dir_to_tangram_outputs
python -u $Script $SPA_FP $OUT_DIR
