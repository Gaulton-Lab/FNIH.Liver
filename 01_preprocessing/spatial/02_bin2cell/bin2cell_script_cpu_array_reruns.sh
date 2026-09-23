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


# activate cond env
source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh

conda activate spatial

Script='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/bin2cell/bin2cell_script.py'
#Script='/tscc/lustre/ddn/scratch/tvashist/charlene/scripts/bin2cell_script.py'
##bin2cell.py two_micron_path source_image_path spatial_dir outdir patient_id mpp

SAMP_MD_TABLE='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/scripts/bin2cell/bin2cell_samp_md_reruns.txt'
#SAMP_MD_TABLE='/tscc/lustre/ddn/scratch/tvashist/charlene/scripts/bin2cell_samp_md.txt'

TASK_ID=${SLURM_ARRAY_TASK_ID}

SAMPLE_LINE=`cat $SAMP_MD_TABLE | sed -n $((TASK_ID + 1))p`

# patient identifier for dataset ex. HL160029, will be used to name directory and output filess
#PATIENT_ID='HL180801'
PATIENT_ID=$(echo $SAMPLE_LINE | cut -f1 -d' ')


# library_id
#LIBRARY_ID='JL_101'
LIBRARY_ID=$(echo $SAMPLE_LINE | cut -f2 -d' ')

# condition
CONDITION=$(echo $SAMPLE_LINE | cut -f3 -d' ')

## brightfield image that is used as input to spaceranger
SOURCE_IMG_PATH=$(echo $SAMPLE_LINE | cut -f4 -d' ')
#SOURCE_IMG_PATH='/tscc/lustre/ddn/scratch/tvashist/charlene/input_images/MASL_HL180801_2nd_reg.tif'

# directory where two micron matrices are located (from spaceranger)
# usually formatted as path1/library_id/outs/binned_outputs/square_002uml
#MIC2_PATH='/tscc/projects/ps-epigen/10x_output/Space_Ranger/outputs4/loupe/JL_101/outs/binned_outputs/square_002um/'
MIC2_PATH=$(echo $SAMPLE_LINE | cut -f5 -d' ')

# directory where spatial files are located (from spaceranger)
# usually formatted as path1/library_id/outs/spatial
SPATIAL_DIR=$(echo $SAMPLE_LINE | cut -f6 -d' ')
#SPATIAL_DIR='/tscc/projects/ps-epigen/10x_output/Space_Ranger/outputs4/loupe/JL_101/outs/spatial'

OUTDIR='/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/outputs/sandbox/spatial/bin2cell_array_cpu/'

#OUTDIR_SAMP="${OUTDIR}${PATIENT_ID}"

# microns per pixel, default is 0.5 but advised to use 0.25 for liver
MPP=$(echo $SAMPLE_LINE | cut -f7 -d' ')
#MPP='0.25'

[[ -d $OUTDIR ]] || mkdir $OUTDIR
[[ -d $OUTDIR_SAMP ]] || mkdir $OUTDIR_SAMP

#echo $PATH
#echo $LD_LIBRARY_PATH

echo $PATIENT_ID
echo $LIBRARY_ID
echo $CONDITION
echo $SOURCE_IMG_PATH
echo $MIC2_PATH
echo $SPATIAL_DIR
echo $OUTDIR
echo $MPP

python $Script $PATIENT_ID $LIBRARY_ID $CONDITION $SOURCE_IMG_PATH $MIC2_PATH $SPATIAL_DIR $OUTDIR $MPP
#python $Script $MIC2_PATH $SOURCE_IMG_PATH $SPATIAL_DIR $OUTDIR $PATIENT_ID $MPP
	      
