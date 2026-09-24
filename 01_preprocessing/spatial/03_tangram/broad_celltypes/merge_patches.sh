#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --cpus-per-task=16
#SBATCH --mem=700G
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=merge_patch
# -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/merge_patches.sh.e
# -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/merge_patches.sh.o
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

#usual in condo
#--cpus-per-task=16
#--mem=700G
#--time=01:00:00
#--account=csd772
#--partition=condo
#--qos=condo

source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh

conda activate spatial_backup
Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/merge_patches.py'
## Check paths before executing
echo "merge patches OUT_DIR: $OUT_DIR"

## usage: python tangram_train_sp.py spatialAnnData scRNAannData admap_output_name adge_output_name
python -u $Script $OUT_DIR
