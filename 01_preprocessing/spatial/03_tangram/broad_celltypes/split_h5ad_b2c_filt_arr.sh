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

### File to subset

## Path to spatial .h5ad from bin2cell 
## ideally later will read the bin2cell with patient id too
TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/tangram/cleaned/241112_split_h5ad_liver_b2c.tsv'

SAMPLE_LINE=`cat $SAMPLE_MD_TABLE | sed -n $((TASK_ID + 1))p`


SPA_FP=$(echo "$SAMPLE_LINE" | cut -f3 ) #using filtered b2c output

## output path to put tangram outputs
TARG_DIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/sandbox/tangram/deep_seq/'

## sub-directory name
PATIENT_ID=$(echo "$SAMPLE_LINE" | cut -f4) #using filtered b2c output, patient_id

OUT_DIR=${TARG_DIR}${PATIENT_ID}"/"


echo "master_pipe SPA_FP: $SPA_FP"
echo "master_pipe TARG_DIR: $TARG_DIR"
echo "master_pipe PATIENT_ID: $PATIENT_ID"
echo "master_pipe OUT_DIR: $OUT_DIR"


[[ -d $TARG_DIR ]] || mkdir $TARG_DIR
[[ -d $OUT_DIR ]] || mkdir $OUT_DIR


# Step 1: Run the first job script to split files
echo "Running the first job script to split files..."
# about 23 min for some reason
JOBID_1=$(sbatch --export=SPA_FP=${SPA_FP},OUT_DIR=${OUT_DIR} split_h5ad.sh | awk '{print $4}')  # Submit the first job script

