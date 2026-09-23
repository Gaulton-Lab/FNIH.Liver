#!/bin/bash
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


#/tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/pp_arr_heps_dep_${PATIENT_ID}.sh.e-%a
source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh

conda activate spatial_backup

Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/Tangram_postPrediction_arr.py'

TAN_DIR="${OUT_DIR}tangram_runs/"
# location of files processed through tangram

# getting base name of split .h5ad to read in processed tangram h5s
FP="${OUT_DIR}subsets/"

SPA_FP=$(ls "$FP"*.h5ad | sed -n "${SLURM_ARRAY_TASK_ID}p")

FILE_NAME=$(basename "$SPA_FP" .h5ad)

## Check paths before executing
echo "pp TAN_DIR: $TAN_DIR"
echo "pp FP: $FP"
echo "pp FILE_NAME: $FILE_NAME"
echo "pp CELLTYPE_MD_COL: $CELLTYPE_MD_COL"

## usage: python tangram_train_sp.py spatialAnnData scRNAannData admap_output_name adge_output_name
python -u $Script $TAN_DIR $FILE_NAME $CELLTYPE_MD_COL
