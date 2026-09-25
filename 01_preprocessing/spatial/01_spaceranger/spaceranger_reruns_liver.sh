#!/bin/bash
#SBATCH -N 1                    #Number of nodes
#SBATCH -n 1                    #Total number of tasks
#SBATCH -c 18                   #Number of threads per process
#SBATCH -t 07:00:00              #Short for --time walltimelimit
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/spaceranger_liver_rerun_%A_%a.e
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/spaceranger_liver_rerun_%A_%a.o
#SBATCH -p condo                #Partition name
#SBATCH -q condo                #QOS name
#SBATCH -A csd772               #Allocation name
#SBATCH --mail-type ALL         #Optional, Send mail when job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu #Optional, Send mail to this address
#SBATCH --propagate=NONE
#SBATCH --job-name=sr_liv
#
# NOTE: --array is intentionally not set here (SLURM does not evaluate shell
# inside #SBATCH lines, so it can't be computed from the sample sheet at
# submission time). Pass it on the command line instead, sized to match the
# number of data rows (i.e. lines - 1 for the header) in the sample sheet:
#
#   sbatch --array=1-7 spaceranger_reruns_liver_array.sh spaceranger_samples.txt
#
# or, computed automatically:
#
#   N=$(($(wc -l < spaceranger_samples.txt) - 1))
#   sbatch --array=1-${N} spaceranger_reruns_liver_array.sh spaceranger_samples.txt

set -euo pipefail

export PATH=/opt/spaceranger-3.0.0:$PATH

# ---- constants shared by every sample (edit here, not per-row) ----
TRANSCRIPTOME=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A
PROBE_SET=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
OUTDIR=/tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/spaceranger

# ---- sample sheet (tab-separated; see spaceranger_samples.txt) ----
# Columns: id  sample  fastqs  cytaimage  slide  area  image  loupe_alignment
# "NA" in slide / area / loupe_alignment means that flag is omitted for this sample.
SAMPLE_SHEET="${1:?Usage: sbatch --array=1-N spaceranger_reruns_liver_array.sh <sample_sheet.txt>}"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit this script with --array=1-N (see NOTE above)." >&2
    exit 1
fi

# +1 to skip the header line; row 1 of data -> SLURM_ARRAY_TASK_ID=1 -> file line 2
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
ROW="$(sed -n "${LINE_NUM}p" "$SAMPLE_SHEET")"

if [ -z "$ROW" ]; then
    echo "ERROR: no row ${SLURM_ARRAY_TASK_ID} (file line ${LINE_NUM}) in ${SAMPLE_SHEET}." >&2
    exit 1
fi

IFS=$'\t' read -r ID SAMPLE FASTQS CYTAIMAGE SLIDE AREA IMAGE LOUPE_ALIGNMENT <<< "$ROW"

echo "=== Task ${SLURM_ARRAY_TASK_ID} (sample sheet line ${LINE_NUM}): id=${ID} ==="

cd "$OUTDIR"

# Build the argument list, only adding --slide/--area/--loupe-alignment when provided.
args=(
    count
    --id="$ID"
    --sample="$SAMPLE"
    --transcriptome="$TRANSCRIPTOME"
    --fastqs="$FASTQS"
    --probe-set="$PROBE_SET"
    --cytaimage="$CYTAIMAGE"
    --create-bam=true
    --image="$IMAGE"
)
[ "$SLIDE" != "NA" ] && args+=(--slide="$SLIDE")
[ "$AREA" != "NA" ] && args+=(--area="$AREA")
[ "$LOUPE_ALIGNMENT" != "NA" ] && args+=(--loupe-alignment="$LOUPE_ALIGNMENT")

echo "Running: spaceranger ${args[*]}"
spaceranger "${args[@]}"
