#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --account=csd772
#SBATCH --cpus-per-task=16
#SBATCH --mem=700G
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --job-name=merge_patch
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/merge_patches.sh.e
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/merge_patches.sh.o
#SBATCH --mail-type BEGIN,END                      # Optional, Send mail when the job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu      # Optional, Send mail to this address

# NOTE: the original version of this script had "# -e ..." / "# -o ..."
# (missing the "#SBATCH" prefix), so those two lines were plain comments,
# not real SBATCH directives -- this job's logs would have fallen back to
# SLURM's default slurm-<jobid>.out naming instead. Fixed above so logs
# actually land at the intended path.
#
# --mem=700G is the value used historically in condo; 50G is reportedly
# enough if this run is only merging raw counts (i.e. no imputed/large
# per-gene layers to concatenate).

# SLURM job: merges all per-subset Tangram-processed patches for one sample
# back into whole-sample objects (the final step of the pipeline, after the
# split -> Tangram train -> postprocess array stages).
#
# Expected environment (set by the caller via `sbatch --export=...`):
#   OUT_DIR  Per-sample output directory (same one used throughout the
#            pipeline).

source /tscc/nfs/home/cmiciano/miniconda3/etc/profile.d/conda.sh
conda activate spatial_backup

Script='/tscc/projects/ps-epigen/users/cmiciano/useful/spatial/merge_patches.py'

## Check paths before executing
echo "merge patches OUT_DIR: $OUT_DIR"

## Usage: python merge_patches.py out_dir
python -u $Script $OUT_DIR
