#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus=1
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --account=csd772
#SBATCH --partition=rtx3090
#SBATCH --qos=condo-gpu
#SBATCH --job-name=tan_arr
# -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/tans_script_nash_gpu_arr_heps_dep.e-%a
# -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/tans_script_nash_gpu_arr_heps_dep.o-%a
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address
# a40
# for a100 can get away with 10gb mem 1cpu per task with 20k single cell and 10k spatial cells
# rtx6000 usually works

source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh

# for A100 have been specifying 200G memory, 1 cpu per task, 1 cpu, qos condo-gpu, 
conda activate spatial_backup
## take each split .h5ad file and run tangram
Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/tans_tangram_train_sp_gpu_arr.py'

FP="${OUT_DIR}subsets/"

SPA_FP=$(ls "$FP"*.h5ad | sed -n "${SLURM_ARRAY_TASK_ID}p")

base_name=$(basename "$SPA_FP" .h5ad)

# output for tangram runs
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

## usage: python tangram_train_sp.py spatialAnnData scRNAannData admap_output_name adge_output_name
python -u $Script $SPA_FP $SC_FP $AD_MAP_NAME $GE_NAME $MARKERS $POST_RNA_NAME $POST_SP_NAME $CELL_DENSITY
