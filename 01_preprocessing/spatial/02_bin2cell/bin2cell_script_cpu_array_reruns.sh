#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=10:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=b2c_cpu
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/bin2cell_cpu_reruns.e-%a
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/bin2cell_cpu_reruns.o-%a
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu     # Optional, Send mail to this address

# SLURM array job: runs bin2cell_script.py for one library (one array task =
# one line of bin2cell_samp_md_reruns.txt).
#
# Sample sheet columns (tab-separated; see bin2cell_samp_md_reruns.txt),
# 1-indexed to match the file's actual column order:
#   1: library_id                      e.g. QY_2645_1_2_3
#   2: patient_id                      e.g. HL180809
#   3: condition
#   4: image_path                      brightfield image used as spaceranger input
#   5: spaceranger_micron2_path        square_002um binned-outputs dir
#   6: spaceranger_spatial_dir_path    spaceranger "spatial" outputs dir
#   7: mpp                             microns per pixel (0.5 default; 0.25 advised for liver)
#
# FIX (this pass): the sample sheet is tab-delimited, but this script
# previously parsed it with `cut -f<N> -d' '` (splitting on spaces). Since
# the file has no spaces to split on, that returned the entire line as a
# single field for every column, silently breaking every value below. Now
# parsed with a single tab-delimited `read`. This also surfaced that
# PATIENT_ID and LIBRARY_ID were previously assigned from the wrong columns
# (swapped relative to the sheet's real column order above) -- fixed to
# match the actual header. Also removed a broken `mkdir $OUTDIR_SAMP` call
# referencing a variable that was never assigned (its assignment line was
# commented out); bin2cell_script.py already creates its own per-patient
# output directory internally, so this wasn't needed.

# activate conda env
source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh

conda activate spatial

Script='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/bin2cell/bin2cell_script.py'
#Script='/tscc/lustre/ddn/scratch/tvashist/charlene/scripts/bin2cell_script.py'
## Usage: bin2cell_script.py patient_id library_id condition source_img_path mic2_path spatial_dir_path outdir_samp_path mpp

SAMP_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/bin2cell/bin2cell_samp_md_reruns.txt'
#SAMP_MD_TABLE='/tscc/lustre/ddn/scratch/tvashist/charlene/scripts/bin2cell_samp_md.txt'

TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_LINE=$(sed -n "$((TASK_ID + 1))p" "$SAMP_MD_TABLE")  # +1 to skip the header row

IFS=$'\t' read -r LIBRARY_ID PATIENT_ID CONDITION SOURCE_IMG_PATH MIC2_PATH SPATIAL_DIR MPP <<< "$SAMPLE_LINE"
#MPP='0.25'

OUTDIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/outputs/sandbox/spatial/bin2cell_array_cpu/'

[[ -d $OUTDIR ]] || mkdir $OUTDIR

echo $PATIENT_ID
echo $LIBRARY_ID
echo $CONDITION
echo $SOURCE_IMG_PATH
echo $MIC2_PATH
echo $SPATIAL_DIR
echo $OUTDIR
echo $MPP

python $Script $PATIENT_ID $LIBRARY_ID $CONDITION $SOURCE_IMG_PATH $MIC2_PATH $SPATIAL_DIR $OUTDIR $MPP
