#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=heps_tg
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/submit_jobs_dep_liver_tangram_heps.sh.e
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/submit_jobs_dep_liver_tangram_heps.sh.o
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

### File to subset

## Path to spatial .h5ad from bin2cell 
## ideally later will read the bin2cell with patient id too
## this spatial is not filtered, the raw output from b2c

# patients HL180809, HL160029 normal 3-4 batch array

TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/tangram/cleaned/heps_subclust_alldonors_ds10k_ref/250214_heps_donors_h5ad_tangram_alldonors_ds10k_ref.tsv'

SAMPLE_LINE=`cat $SAMPLE_MD_TABLE | sed -n $((TASK_ID + 1))p`

SPA_FP=$(echo "$SAMPLE_LINE" | cut -f1 ) 

## output path to put tangram outputs
TARG_DIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/sandbox/tangram/deep_seq_heps_alldonors_ds10k_ref/'

## sub-directory name
PATIENT_ID=$(echo "$SAMPLE_LINE" | cut -f2)

CELL_DENSITY=$(echo "$SAMPLE_LINE" | cut -f3)

## ran sc .h5ad run through preprocessing, downsampled to 20k cells
## in the sc .h5ad combined T and NK to be T_NK, removed mast and schwann, then downsampled to 20k cells
SC_FP=$(echo "$SAMPLE_LINE" | cut -f4)
# broad celltypes rank genes
MARKERS=$(echo "$SAMPLE_LINE" | cut -f5)

CELLTYPE_MD_COL="cellsubtype"

OUT_DIR=${TARG_DIR}${PATIENT_ID}"/"

JOB_OUT='/tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/' #logs for jobs assumes this exists already

echo "master_pipe SPA_FP: $SPA_FP"
echo "master_pipe TARG_DIR: $TARG_DIR"
echo "master_pipe PATIENT_ID: $PATIENT_ID"
echo "master_pipe CELL_DENSITY: $CELL_DENSITY"
echo "master_pipe SC_FP: $SC_FP"
echo "master_pipe MARKERS: $MARKERS"
echo "master_pipe CELLTYPE_MD_COL: $CELLTYPE_MD_COL"
echo "master_pipe OUT_DIR: $OUT_DIR"
echo "master_pipe JOB_OUT: $JOB_OUT"


[[ -d $TARG_DIR ]] || mkdir $TARG_DIR
[[ -d $OUT_DIR ]] || mkdir $OUT_DIR



# Get the number of files in outdir/subsets
NUM_FILES=$(ls "${OUT_DIR}subsets" | wc -l)
if [ "$NUM_FILES" -eq 0 ]; then
    echo "No files found in $OUT_DIR/subsets"
    exit 1
fi


#echo "num files main script: $NUM_FILES"
# Submit the first job array with --array=1-14 and capture the job ID --dependency=afterok:$JOBID_1
# 30 min in gpu, 15 for each dataset

## can't use backslashes unfortunately
echo "Submitting Tangram run"
JOBID_2=$(sbatch --export=OUT_DIR=${OUT_DIR},SC_FP=${SC_FP},MARKERS=${MARKERS},CELL_DENSITY=${CELL_DENSITY} --output="${JOB_OUT}tans_script_gpu_arr_${PATIENT_ID}.o-%a" --error="${JOB_OUT}tans_script_gpu_arr_${PATIENT_ID}.e-%a" --array=1-"$NUM_FILES" tans_tangram_train_sp_gpu_arr.sh | awk '{print $4}')
# manually 49-52 for broad
#JOBID_2=$(sbatch --export=OUT_DIR=${OUT_DIR},SC_FP=${SC_FP},MARKERS=${MARKERS},CELL_DENSITY=${CELL_DENSITY} --output="${JOB_OUT}tans_script_gpu_arr${PATIENT_ID}.o-%a" --error="${JOB_OUT}tans_script_gpu_arr${PATIENT_ID}.e-%a" --array=5,6,19,20 tans_tangram_train_sp_gpu_arr.sh | awk '{print $4}')
echo "JOBID_2: ${JOBID_2}"

# Submit the second job array with --array=1-14, with a dependency on the first job array # --dependency=afterok:$JOBID_2
# 13-19  min post prediction
echo "Submitting postprocessing run"
JOBID_3=$(sbatch --export=OUT_DIR=${OUT_DIR},CELLTYPE_MD_COL=${CELLTYPE_MD_COL},JOB_OUT=${JOB_OUT},PATIENT_ID=${PATIENT_ID} --error="${JOB_OUT}pp_arr_${PATIENT_ID}.e-%a" --output="${JOB_OUT}pp_arr_${PATIENT_ID}.o-%a" --dependency=afterok:$JOBID_2 --array=1-"$NUM_FILES" Tangram_postPrediction_arr.sh | awk '{print $4}')
#JOBID_3=$(sbatch --export=OUT_DIR=${OUT_DIR},CELLTYPE_MD_COL=${CELLTYPE_MD_COL},JOB_OUT=${JOB_OUT},PATIENT_ID=${PATIENT_ID} --output="${JOB_OUT}pp_arr_${PATIENT_ID}.e-%a" --output="${JOB_OUT}pp_arr_${PATIENT_ID}.o-%a" --array=1-"$NUM_FILES"  Tangram_postPrediction_arr.sh | awk '{print $4}')

# --dependency=afterok:$JOBID_2 
echo "JOBID_3: ${JOBID_3}"

# Step 4: Once the third job array completes, run the fourth script to merge all files
# 15-20 min to merge back together
JOBID_4=$(sbatch --export=OUT_DIR=${OUT_DIR}  --output="${JOB_OUT}merge_patches_${PATIENT_ID}.sh.o" --error="${JOB_OUT}merge_patches_${PATIENT_ID}.sh.e" --dependency=afterok:$JOBID_3 merge_patches.sh | awk '{print $4}')
#JOBID_4=$(sbatch --export=OUT_DIR=${OUT_DIR}  --output="${JOB_OUT}merge_patches_${PATIENT_ID}.sh.o" --error="${JOB_OUT}merge_patches_${PATIENT_ID}.sh.e" merge_patches.sh | awk '{print $4}')
echo "JOBID_4: ${JOBID_4}"

echo "All jobs submitted. Final merge will run after all tasks are completed."
# sacct -j 2787729 --format=State --noheader | head -1

