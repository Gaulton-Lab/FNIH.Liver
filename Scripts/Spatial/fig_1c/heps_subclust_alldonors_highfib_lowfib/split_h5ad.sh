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

source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh

conda activate spatial_backup
Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/split_h5ad.py'

## Check paths before executing
echo "split_h5ad FP: $SPA_FP"
echo "split_h5ad OUT_DIR: $OUT_DIR"

## usage: python tangram_train_sp.py spatialAnnData scRNAannData admap_output_name adge_output_name
python -u $Script $SPA_FP $OUT_DIR
