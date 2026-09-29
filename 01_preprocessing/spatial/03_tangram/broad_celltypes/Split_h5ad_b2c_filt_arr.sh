#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=sp_nf
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/sub_jobs_liver_b2c_split_filt.sh.e-%a
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/sub_jobs_liver_b2c_split_filt.sh.o-%a
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

# SLURM array entry point (broad cell types): for one sample (one array task
# = one line of the sample sheet below), submits split_h5ad.sh to split that
# sample's filtered Bin2cell spatial object into Tangram-sized chunks.
#
# This is the "step 1 of 4" launcher for the broad-cell-type Tangram
# pipeline; the remaining 3 stages (Tangram train, postprocess, merge) are
# submitted by the per-condition Submit_jobs_dep_*_liver_no_mast_filt_b2c.sh
# drivers once this step's outputs exist.
#
# Sample sheet columns used (tab-separated, 1-indexed via `cut -f`):
#   3: SPA_FP     (tsv column "spatial_liver_b2c_path_after_filtering") path to
#                the filtered Bin2cell (b2c) spatial .h5ad
#   4: PATIENT_ID (tsv column "patient_id_filt") output subdirectory name for this sample

TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/tangram/cleaned/241112_split_h5ad_liver_b2c.tsv'
SAMPLE_LINE=$(sed -n "$((TASK_ID + 1))p" "$SAMPLE_MD_TABLE")  # +1 to skip the header row

SPA_FP=$(echo "$SAMPLE_LINE" | cut -f3)   # filtered b2c spatial output
PATIENT_ID=$(echo "$SAMPLE_LINE" | cut -f4)

## Output path to put Tangram outputs under.
TARG_DIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/sandbox/tangram/deep_seq/'
OUT_DIR=${TARG_DIR}${PATIENT_ID}"/"

echo "master_pipe SPA_FP: $SPA_FP"
echo "master_pipe TARG_DIR: $TARG_DIR"
echo "master_pipe PATIENT_ID: $PATIENT_ID"
echo "master_pipe OUT_DIR: $OUT_DIR"

[[ -d $TARG_DIR ]] || mkdir $TARG_DIR
[[ -d $OUT_DIR ]] || mkdir $OUT_DIR

# Step 1: split this sample's spatial object into Tangram-sized chunks.
echo "Running the first job script to split files..."
# Historically ~23 min for this step.
JOBID_1=$(sbatch --export=SPA_FP=${SPA_FP},OUT_DIR=${OUT_DIR} split_h5ad.sh | awk '{print $4}')
