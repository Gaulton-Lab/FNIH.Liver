#!/usr/bin/env bash
# Generate per-condition (Control/MASH/MASL/MetALD) scHiCluster loop-calling
# scaffolds for every cell type, mirroring the existing Hepatocytes run.
#
# Reproduces:  impute/10K/Hepatocytes/loop/Hepatocytes_<cond>/  ->  <CellType>/loop/<CellType>_<cond>/
# Each group is one pseudobulk = all donors of that cell type with that clinical condition.
#
# Path duality: this script may run on the mediator node (project tree mounted at
# /mnt/tscc2/...) but the generated files must reference the TSCC compute-node path
# (/tscc/projects/ps-renlab2/...). WRITE_BASE = where files are written/read now;
# PATH_BASE = the path string embedded INSIDE the generated files (valid at run time).
# If you run this directly on a TSCC node, set WRITE_BASE=PATH_BASE.
#
# After running this, launch loop calling on a compute node (schicluster env) via the
# generated impute/10K/<CellType>/loop/run_all.sh or .sbatch scripts.
set -euo pipefail

WRITE_BASE=./impute/10K
PATH_BASE=./impute/10K

COOL_FS=$WRITE_BASE/cool                                   # filesystem (existence checks, donor map)
COOL_PATH=$PATH_BASE/cool                                  # embedded in cell_table.csv
CHROM=~/genome_ref/hg38.main.chrom.sizes   # embedded
CHROM_FS=~/genome_ref/hg38.main.chrom.sizes      # check now
BLACKLIST=$PATH_BASE/loop/hg38.main.rowsum300.10kb.bed     # embedded
BLACKLIST_FS=$WRITE_BASE/loop/hg38.main.rowsum300.10kb.bed # check now
RES=10000
CPU=16
MIN_DONORS=${MIN_DONORS:-5}          # groups with fewer donors are scaffolded but left OUT of the run driver
CELLTYPES="B Cholangiocyte Endothelial HSC Mast Myeloid NK Schwann T"   # Hepatocytes already done

# ---- sanity (against the now-accessible filesystem) ----
[ -f "$CHROM_FS" ]     || { echo "ERROR: chrom.sizes not found: $CHROM_FS" >&2; exit 1; }
[ -f "$BLACKLIST_FS" ] || { echo "ERROR: blacklist not found: $BLACKLIST_FS" >&2; exit 1; }
[ -d "$COOL_FS/Hepatocytes" ] || { echo "ERROR: cool dir not found: $COOL_FS/Hepatocytes" >&2; exit 1; }

# ---- donor -> condition map, derived from the existing Hepatocytes condition-split cools ----
DC=$(mktemp)
ls "$COOL_FS/Hepatocytes/" \
  | sed -nE 's/^(HL[0-9]+)_Hepatocytes_([A-Za-z-]+)_10kb\.cool$/\1\t\2/p' \
  | sort -u > "$DC"
echo "donor->condition map: $(wc -l < "$DC") donors"

write_snakefile_master () {  # $1 = groupdir (filesystem path to write into)
  cat > "$1/Snakefile_master" <<EOF
cpu = $CPU
resolution = $RES
chrom_size_path = "$CHROM"
black_list_path = "$BLACKLIST"

from schicluster.loop import call_loop

rule call_loop:
    input:
        'cell_table.csv'
    output:
        'Success'
    run:
        call_loop(
            cell_table_path='cell_table.csv',
            output_dir='./',
            chrom_size_path=chrom_size_path,
            shuffle=True,
            chunk_size=200,
            dist=5050000,
            cap=5,
            pad=5,
            gap=2,
            resolution=resolution,
            min_cutoff=1e-06,
            keep_cell_matrix=False,
            cpu_per_job=cpu,
            log_e=True,
            raw_resolution_str=None,
            downsample_shuffle=500,
            black_list_path=black_list_path,
            fdr_pad=7,
            fdr_min_dist=5,
            fdr_max_dist=500,
            fdr_thres=0.1,
            dist_thres=20000,
            size_thres=1)
EOF
}

for CT in $CELLTYPES; do
  LOOPDIR_FS=$WRITE_BASE/$CT/loop
  LOOPDIR_PATH=$PATH_BASE/$CT/loop
  mkdir -p "$LOOPDIR_FS"
  CMDS=$LOOPDIR_FS/snakemake_cmds.txt        ; : > "$CMDS"
  RUNALL=$LOOPDIR_FS/run_all.sh              ; : > "$RUNALL"
  SKIP=$LOOPDIR_FS/SKIPPED_small_groups.txt  ; : > "$SKIP"
  printf '#!/usr/bin/env bash\nset -eo pipefail\nsource /tscc/nfs/home/y2xie/anaconda3/etc/profile.d/conda.sh\nconda activate schicluster\n' > "$RUNALL"

  for COND in Control MASH MASL MetALD; do
    GROUP=${CT}_${COND}
    GDIR_FS=$LOOPDIR_FS/$GROUP
    GDIR_PATH=$LOOPDIR_PATH/$GROUP
    mkdir -p "$GDIR_FS"
    CT_TABLE=$GDIR_FS/cell_table.csv ; : > "$CT_TABLE"
    n=0
    while IFS=$'\t' read -r donor cond; do
      [ "$cond" = "$COND" ] || continue
      [ -f "$COOL_FS/${donor}_${CT}_10kb.cool" ] || continue
      printf '%s,%s,%s\n' "${donor}_${GROUP}" "$COOL_PATH/${donor}_${CT}_10kb.cool" "$GROUP" >> "$CT_TABLE"
      n=$((n+1))
    done < "$DC"
    write_snakefile_master "$GDIR_FS"

    cmd="snakemake -d $GDIR_PATH -s $GDIR_PATH/Snakefile_master -j $CPU"
    if [ "$n" -ge "$MIN_DONORS" ]; then
      echo "$cmd" >> "$CMDS"
      echo "$cmd" >> "$RUNALL"
    else
      echo "$GROUP  donors=$n  (< $MIN_DONORS, excluded from run driver)" >> "$SKIP"
    fi
    printf '  %-22s donors=%s\n' "$GROUP" "$n"
  done
  chmod +x "$RUNALL"
  echo "[$CT] driver: $CMDS  ($(wc -l < "$CMDS") groups to run; skipped: $(wc -l < "$SKIP"))"
done

rm -f "$DC"
echo "DONE. Files written under $WRITE_BASE ; paths inside them reference $PATH_BASE"
