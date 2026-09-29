#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=spl_hep
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/sub_jobs_heps_split_filt.sh.e-%a
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/sub_jobs_heps_split_filt.sh.o-%a
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

# SLURM array entry point (hepatocyte subtypes, highfib_lowfib variant): for
# one sample (one array task = one line of the sample sheet below), submits
# split_h5ad.sh to split that sample's spatial object into Tangram-sized
# chunks.
#
# This is "step 1" for this pipeline variant; Submit_jobs_dep_liver_tangram_heps.sh
# (in this same folder) submits the remaining 3 stages (Tangram train,
# postprocess, merge) once this step's outputs exist.
#
# Sample sheet columns used (tab-separated, 1-indexed via `cut -f`):
#   1: SPA_FP     path to the (n_counts-filtered) spatial .h5ad for this sample
#   2: PATIENT_ID output subdirectory name for this sample

TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/tangram/cleaned/heps_subclust_alldonors_highfib_lowfib/250220_heps_h5ad_tangram_alldonors_highfib_lowfib.tsv'
SAMPLE_LINE=$(sed -n "$((TASK_ID + 1))p" "$SAMPLE_MD_TABLE")  # +1 to skip the header row

SPA_FP=$(echo "$SAMPLE_LINE" | cut -f1)     # filtered n_counts spatial output
PATIENT_ID=$(echo "$SAMPLE_LINE" | cut -f2) # filtered n_counts spatial output, patient_id

## Output path to put Tangram outputs under.
TARG_DIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/sandbox/tangram/deep_seq_heps_alldonors_highfib_lowfib/'
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
