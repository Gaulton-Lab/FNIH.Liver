#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=masl_tg
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/sub_jobs_liver_masl_filt_ds.sh.e
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/sub_jobs_liver_masl_filt_ds.sh.o
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

# Driver for the broad-cell-type Tangram pipeline, MASL condition.
#
# This submits (a subset of) the 3 remaining pipeline stages for one sample
# (one array task = one line of the sample sheet), after Split_h5ad_b2c_filt_arr.sh
# has already produced OUT_DIR/subsets/:
#   Stage 2: Tangram train/map array   (tans_tangram_train_sp_gpu_arr.sh)
#   Stage 3: Tangram postprocess array (Tangram_postPrediction_arr.sh, depends on stage 2)
#   Stage 4: merge patches             (merge_patches.sh, depends on stage 3)
#
# PIPELINE-STAGE TOGGLE: previously only stage 4 (merge) was active here,
# assuming stages 2-3 had already been run separately. All 4 stages are now
# active below, with the full dependency chain (stage 3 waits on stage 2,
# stage 4 waits on stage 3) -- comment stages back out if you want to run
# them separately again (e.g. to rerun just the merge for a sample whose
# earlier stages already completed).
#
# Cleanup note: removed a second, hardcoded-array-range dead alternate for
# stage 2 (a one-off manual override for reprocessing a handful of files,
# e.g. "--array=1-4") and a stage-3 dead alternate that had a real bug
# (--output passed twice, once with an ".e-%a" suffix, instead of --error).
#
# Sample sheet columns used (tab-separated, 1-indexed via `cut -f`):
#   3: SPA_FP        (tsv column "spatial_liver_b2c_path_after_filtering") path to
#                    the filtered Bin2cell (b2c) spatial .h5ad (unused directly
#                    here, since subsets/ already exists by this point)
#   4: PATIENT_ID    (tsv column "patient_id_filt") output subdirectory name for this sample
#   5: CELL_DENSITY  (tsv column "cell_density") cell-density argument passed to the Tangram training script

TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/tangram/cleaned/241112_split_h5ad_liver_b2c.tsv'
SAMPLE_LINE=$(sed -n "$((TASK_ID + 1))p" "$SAMPLE_MD_TABLE")  # +1 to skip the header row

SPA_FP=$(echo "$SAMPLE_LINE" | cut -f3)  # filtered b2c spatial output (column 3, "..._after_filtering", per the sample sheet)
PATIENT_ID=$(echo "$SAMPLE_LINE" | cut -f4)
CELL_DENSITY=$(echo "$SAMPLE_LINE" | cut -f5)

## Output path to put Tangram outputs under.
TARG_DIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/sandbox/tangram/deep_seq/'
OUT_DIR=${TARG_DIR}${PATIENT_ID}"/"

## Reference scRNA object (preprocessed, downsampled to 20k cells; T and NK
## combined into T_NK, mast and Schwann removed) and matching broad-cell-type
## marker/rank-genes directory, both from the MASL preprocessing run.
SC_FP='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/outputs/sandbox/spatial/20241112_tangram_preprocessing_genelist50_liver_masl_rem_mast/Preprocessed_adsc.h5ad'
MARKERS='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/outputs/sandbox/spatial/20241112_tangram_preprocessing_genelist50_liver_masl_rem_mast/markers'

CELLTYPE_MD_COL="celltype"
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
