#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus=1
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --account=csd772
#SBATCH --partition=rtx6000
#SBATCH --qos=condo-gpu
#SBATCH --job-name=tan_arr
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/tans_script_nash_gpu_arr_heps_dep.e-%a
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/tans_script_nash_gpu_arr_heps_dep.o-%a
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

# NOTE: as in merge_patches.sh, the original "# -e" / "# -o" lines here were
# missing the "#SBATCH" prefix and so were not actually taking effect; fixed
# above.
#
# GPU/memory notes from prior runs, kept for reference when re-tuning:
#   - rtx6000: usually works.
#   - a100: 10G mem, 1 cpu-per-task is enough for ~20k single-cell x 10k
#     spatial cells (the settings used above).
#   - a100 (larger runs): have also needed to specify 200G mem, 1 cpu, qos
#     condo-gpu.

# SLURM array job: runs Tangram training + mapping for one spatial subset
# (one array task per file under OUT_DIR/subsets/, produced by Split_h5ad.py).
#
# Expected environment (set by the caller via `sbatch --export=...`):
#   OUT_DIR       Per-sample output directory; subsets/ and tangram_runs/
#                 live under it.
#   SC_FP         Path to the preprocessed reference scRNA .h5ad.
#   MARKERS       Path to the marker-gene list used for Tangram training.
#   CELL_DENSITY  Cell-density argument passed through to the Tangram script.

source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh
conda activate spatial_backup

## Take each split .h5ad file and run tangram
Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/tans_tangram_train_sp_gpu_arr.py'

FP="${OUT_DIR}subsets/"
SPA_FP=$(ls "$FP"*.h5ad | sed -n "${SLURM_ARRAY_TASK_ID}p")
base_name=$(basename "$SPA_FP" .h5ad)

# Output directory for Tangram runs.
FP_TAN_OUT_PATH="${OUT_DIR}tangram_runs/"
[[ -d $FP_TAN_OUT_PATH ]] || mkdir $FP_TAN_OUT_PATH

AD_MAP_NAME="${FP_TAN_OUT_PATH}admap_nash_trim_list_gpu_${base_name}"
GE_NAME="${FP_TAN_OUT_PATH}adge_nash_trim_list_gpu_${base_name}"
POST_RNA_NAME="${FP_TAN_OUT_PATH}Post_scRNA_${base_name}"
POST_SP_NAME="${FP_TAN_OUT_PATH}Post_sp_${base_name}"

## Check paths before executing
echo "tangram SPA_FP: $SPA_FP"
echo "tangram SC_FP: $SC_FP"
echo "tangram MARKERS: $MARKERS"
echo "tangram POST_RNA_NAME: $POST_RNA_NAME"
echo "tangram POST_SP_NAME: $POST_SP_NAME"
echo "tangram CELL DENSITY: $CELL_DENSITY"

## Usage: python tans_tangram_train_sp_gpu_arr.py spatial_fp sc_fp admap_name adge_name markers post_rna_name post_sp_name cell_density
python -u $Script $SPA_FP $SC_FP $AD_MAP_NAME $GE_NAME $MARKERS $POST_RNA_NAME $POST_SP_NAME $CELL_DENSITY
