#!/usr/bin/env bash
# Emit one sbatch script per (cell type x condition) loop-calling group, plus a
# submit_all.sh driver. Headers copied from the project's existing TSCC job scripts
# (-p condo -q condo -A csd788). Each job runs the group's Snakefile_master, which
# runs the whole scHiCluster call_loop pipeline end-to-end.
#
# WRITE_BASE = where .sbatch files are written now (mediator mount).
# PATH_BASE  = TSCC compute-node path embedded inside the scripts. Set equal if run on TSCC.
set -euo pipefail

WRITE_BASE=./impute/10K
PATH_BASE=./impute/10K
CELLTYPES="B Cholangiocyte Endothelial HSC Mast Myeloid NK Schwann T"
CPU=16
MEM=64G
TIME=2-00:00:00
PART=condo
QOS=condo
ACCT=csd788

SUBMIT=$WRITE_BASE/submit_all_percondition_loops.sh
printf '#!/usr/bin/env bash\n# Submit every per-condition loop job. Run from a TSCC login node (has sbatch).\nset -euo pipefail\n\n' > "$SUBMIT"

njobs=0
for CT in $CELLTYPES; do
  CMDS=$WRITE_BASE/$CT/loop/snakemake_cmds.txt
  [ -s "$CMDS" ] || continue
  # snakemake_cmds.txt already excludes groups below MIN_DONORS; its paths are PATH_BASE-rooted.
  while read -r _sm _d GDIR_PATH _rest; do
    GROUP=$(basename "$GDIR_PATH")
    GDIR_FS=$WRITE_BASE/$CT/loop/$GROUP
    SB_FS=$GDIR_FS/run_loop.sbatch
    SB_PATH=$PATH_BASE/$CT/loop/$GROUP/run_loop.sbatch
    cat > "$SB_FS" <<EOF
#!/bin/bash
#SBATCH -J loop_${GROUP}
#SBATCH -p ${PART}
#SBATCH -q ${QOS}
#SBATCH -A ${ACCT}
#SBATCH -N 1
#SBATCH -c ${CPU}
#SBATCH --mem=${MEM}
#SBATCH -t ${TIME}
#SBATCH -o ${GDIR_PATH}/slurm-%j.out
#SBATCH -e ${GDIR_PATH}/slurm-%j.err
set -eo pipefail   # NOT -u: conda's compiler (de)activation scripts use unbound vars
source /tscc/nfs/home/y2xie/anaconda3/etc/profile.d/conda.sh
conda activate schicluster
snakemake -d ${GDIR_PATH} -s ${GDIR_PATH}/Snakefile_master -j ${CPU} --rerun-incomplete
EOF
    echo "sbatch ${SB_PATH}" >> "$SUBMIT"
    njobs=$((njobs+1))
  done < "$CMDS"
done
chmod +x "$SUBMIT"
echo "Wrote $njobs sbatch scripts under $WRITE_BASE ; paths inside them reference $PATH_BASE"
echo "Submit driver: $SUBMIT"
