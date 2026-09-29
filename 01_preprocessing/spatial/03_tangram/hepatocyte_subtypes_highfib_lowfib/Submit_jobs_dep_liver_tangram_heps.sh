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

# Driver for the hepatocyte-subtype Tangram pipeline, highfib_lowfib variant
# (reference scRNA subset selected on high-fibrosis / low-fibrosis donors).
#
# This submits the 3 remaining pipeline stages for one sample (one array
# task = one line of the sample sheet). Unlike the broad-cell-type pipeline,
# the split step here is a separate one-off (see Submit_jobs_split_h5ad_arr.sh
# in the highfib_lowfib sibling folder) rather than launched from this file:
#   Stage 2: Tangram train/map array   (tans_tangram_train_sp_gpu_arr.sh)
#   Stage 3: Tangram postprocess array (Tangram_postPrediction_arr.sh, depends on stage 2)
#   Stage 4: merge patches             (merge_patches.sh, depends on stage 3)
#
# PIPELINE-STAGE TOGGLE: see the broad_celltypes/Submit_jobs_dep_*.sh drivers
# for the general pattern. All 3 stages are active below with their
# dependency chain intact.
#
# Cleanup note: removed a stage-2 dead alternate that hardcoded a one-off
# array subset (`--array=5,6,19,20`, for reprocessing specific patients) and
# a stage-3 dead alternate that only differed from the active line by
# omitting the useful "JOBID_3: ..." echo below it.
#
# Sample sheet columns used (tab-separated, 1-indexed via `cut -f`):
#   1: SPA_FP        path to the spatial .h5ad for this sample
#   2: PATIENT_ID    output subdirectory name for this sample
#   3: CELL_DENSITY  cell-density argument passed to the Tangram training script
#   4: SC_FP         path to the preprocessed reference scRNA .h5ad
#   5: MARKERS       path to the marker-gene list used for Tangram training

TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/tangram/cleaned/heps_subclust_alldonors_highfib_lowfib/250220_heps_h5ad_tangram_alldonors_highfib_lowfib.tsv'
SAMPLE_LINE=$(sed -n "$((TASK_ID + 1))p" "$SAMPLE_MD_TABLE")  # +1 to skip the header row

SPA_FP=$(echo "$SAMPLE_LINE" | cut -f1)
PATIENT_ID=$(echo "$SAMPLE_LINE" | cut -f2)
CELL_DENSITY=$(echo "$SAMPLE_LINE" | cut -f3)
SC_FP=$(echo "$SAMPLE_LINE" | cut -f4)
MARKERS=$(echo "$SAMPLE_LINE" | cut -f5)

## Output path to put Tangram outputs under.
TARG_DIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/sandbox/tangram/deep_seq_heps_alldonors_highfib_lowfib/'
OUT_DIR=${TARG_DIR}${PATIENT_ID}"/"

CELLTYPE_MD_COL="cellsubtype"
JOB_OUT='/tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/' # logs for jobs, assumes this dir already exists

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

# Number of split files already present in OUT_DIR/subsets (one array task
# per subset for stages 2-3 below).
NUM_FILES=$(ls "${OUT_DIR}subsets" | wc -l)
if [ "$NUM_FILES" -eq 0 ]; then
    echo "No files found in $OUT_DIR/subsets"
    exit 1
fi

# --- Stage 2: Tangram train/map (one array task per subset) ---
# ~30 min on GPU, ~15 min per dataset historically.
# NOTE: sbatch --export can't use backslash line continuations, so each
# call below is kept on one line.
echo "Submitting Tangram run"
JOBID_2=$(sbatch --export=OUT_DIR=${OUT_DIR},SC_FP=${SC_FP},MARKERS=${MARKERS},CELL_DENSITY=${CELL_DENSITY} --output="${JOB_OUT}tans_script_gpu_arr_${PATIENT_ID}.o-%a" --error="${JOB_OUT}tans_script_gpu_arr_${PATIENT_ID}.e-%a" --array=1-"$NUM_FILES" tans_tangram_train_sp_gpu_arr.sh | awk '{print $4}')
echo "JOBID_2: ${JOBID_2}"

# --- Stage 3: Tangram postprocess (one array task per subset; depends on stage 2) ---
# ~13-19 min historically.
echo "Submitting postprocessing run"
JOBID_3=$(sbatch --export=OUT_DIR=${OUT_DIR},CELLTYPE_MD_COL=${CELLTYPE_MD_COL},JOB_OUT=${JOB_OUT},PATIENT_ID=${PATIENT_ID} --error="${JOB_OUT}pp_arr_${PATIENT_ID}.e-%a" --output="${JOB_OUT}pp_arr_${PATIENT_ID}.o-%a" --dependency=afterok:$JOBID_2 --array=1-"$NUM_FILES" Tangram_postPrediction_arr.sh | awk '{print $4}')
echo "JOBID_3: ${JOBID_3}"

# --- Stage 4: merge patches back into whole-sample objects (depends on stage 3) ---
# ~15-20 min historically.
JOBID_4=$(sbatch --export=OUT_DIR=${OUT_DIR} --output="${JOB_OUT}merge_patches_${PATIENT_ID}.sh.o" --error="${JOB_OUT}merge_patches_${PATIENT_ID}.sh.e" --dependency=afterok:$JOBID_3 merge_patches.sh | awk '{print $4}')
echo "JOBID_4: ${JOBID_4}"

echo "All jobs submitted. Final merge will run after all tasks are completed."
# sacct -j 2787729 --format=State --noheader | head -1
