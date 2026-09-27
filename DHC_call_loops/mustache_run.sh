#!/usr/bin/env bash
# Run mustache on all cool/mcool files (celltype + disease) at 10kb.
source /home/y2xie/miniconda3/etc/profile.d/conda.sh
conda activate peakachu
set -eo pipefail

ROOT_DIS=./04.matrices/raw/disease
ODIR=./04.matrices/raw/disease/loop/mustache_new
LDIR=$ODIR/logs
CT_COOLS=$ODIR/cools_10kb   # extracted from mcool via cooler cp
mkdir -p "$ODIR" "$LDIR"

RES=10000
PT=0.1     # p-value threshold (mustache default 0.2; 0.1 is stricter and matches prior runs)

# 10 celltype 10kb cools (extracted from mcools)
declare -a CT_FILES
for cl in "$CT_COOLS"/*_10000.cool; do
  base=$(basename "$cl" _10000.cool)
  CT_FILES+=("${base}|${cl}")
done

# 41 disease 10kb cools
declare -a DIS_FILES
for cl in "$ROOT_DIS"/*_10000.cool; do
  base=$(basename "$cl" _10000.cool)
  DIS_FILES+=("${base}_disease|${cl}")
done

ALL=("${CT_FILES[@]}" "${DIS_FILES[@]}")
echo "[$(date '+%F %T')] mustache pipeline -- ${#ALL[@]} files" | tee "$LDIR/run.log"

for entry in "${ALL[@]}"; do
  IFS='|' read -r label path <<< "$entry"
  out=$ODIR/${label}.mustache.tsv
  log=$LDIR/${label}.log
  if [[ -s "$out" ]]; then
    echo "[$(date '+%F %T')] $label  already done" | tee -a "$LDIR/run.log"
    continue
  fi
  echo "[$(date '+%F %T')] $label  running" | tee -a "$LDIR/run.log"
  set +e
  mustache -f "$path" -r $RES -pt $PT -p 4 -o "$out" > "$log" 2>&1
  rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then
    echo "[$(date '+%F %T')] $label  FAILED rc=$rc -- see $log" | tee -a "$LDIR/run.log"
  else
    n=$(($(wc -l < "$out" 2>/dev/null || echo 1) - 1))
    echo "[$(date '+%F %T')] $label  done -> $n loops" | tee -a "$LDIR/run.log"
  fi
done

echo "[$(date '+%F %T')] ALL DONE" | tee -a "$LDIR/run.log"
