#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=pp_arr
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

# SLURM array job: post-processes one Tangram output chunk (one array task
# per subset produced by Split_h5ad.py / tans_tangram_train_sp_gpu_arr.sh).
#
# Expected environment (set by the caller via `sbatch --export=...`):
#   OUT_DIR          Per-sample output directory (same one used throughout
#                     the pipeline); tangram_runs/ and subsets/ live under it.
#   CELLTYPE_MD_COL  Name of the .obs column holding cell-type labels in the
#                     reference scRNA object (e.g. "celltype" or "cellsubtype").
#
# --output / --error for this array job are passed explicitly by the caller
# on the sbatch command line (e.g. .../pp_arr_<PATIENT_ID>.o-%a), rather than
# being set here via #SBATCH, so each patient's run gets its own log files
# without editing this script.

source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh
conda activate spatial_backup

Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/Tangram_postPrediction_arr.py'

# Directory of per-subset Tangram outputs (admap/adge/Post_scRNA/Post_sp files).
TAN_DIR="${OUT_DIR}tangram_runs/"

# Directory of split spatial .h5ad subsets; used to identify which subset
# file this array task index corresponds to.
FP="${OUT_DIR}subsets/"

SPA_FP=$(ls "$FP"*.h5ad | sed -n "${SLURM_ARRAY_TASK_ID}p")
FILE_NAME=$(basename "$SPA_FP" .h5ad)

## Check paths before executing
echo "pp TAN_DIR: $TAN_DIR"
echo "pp FP: $FP"
echo "pp FILE_NAME: $FILE_NAME"
echo "pp CELLTYPE_MD_COL: $CELLTYPE_MD_COL"

## Usage: python Tangram_postPrediction_arr.py tan_dir file_name celltype_md_col
python -u $Script $TAN_DIR $FILE_NAME $CELLTYPE_MD_COL
