#!/usr/bin/env bash
# Run HiCCUPS on all .hic files (celltype + Hep disease) at 10kb.
# Uses juicer_tools 1.22.01 in --cpu mode.

set -eo pipefail
JAR=~/packages/juicer/juicer_tools_1.22.01.jar
ROOT_CT=./04.matrices/raw/celltype
ROOT_DIS=./04.matrices/raw/disease
ODIR=./04.matrices/raw/disease/loop/hiccups_new
LDIR=$ODIR/logs
mkdir -p "$ODIR" "$LDIR"

# Standard HiCCUPS at 10kb (matches Rao 2014 paper recommendations for low-depth)
RES=10000
HICCUPS_OPTS="--cpu --threads 4 --ignore-sparsity -r $RES -f 0.1 -p 4 -i 7 -t 0.02,1.5,1.75,2 -d 20000"

# All 15 .hic files
FILES=(
  "$ROOT_CT"/*_10000.hic
  "$ROOT_DIS"/*_10000.hic
)

echo "[$(date '+%F %T')] HiCCUPS pipeline -- ${#FILES[@]} files" | tee "$LDIR/run.log"

for hic in "${FILES[@]}"; do
  base=$(basename "$hic" _10000.hic)
  out=$ODIR/$base
  log=$LDIR/${base}.log
  # output is a directory; if merged_loops.bedpe exists, skip
  if [[ -s "$out/merged_loops.bedpe" ]]; then
    echo "[$(date '+%F %T')] $base  already done, skipping" | tee -a "$LDIR/run.log"
    continue
  fi
  rm -rf "$out"
  echo "[$(date '+%F %T')] $base  running" | tee -a "$LDIR/run.log"
  set +e
  java -Xmx32g -jar "$JAR" hiccups $HICCUPS_OPTS "$hic" "$out" > "$log" 2>&1
  rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then
    echo "[$(date '+%F %T')] $base  FAILED rc=$rc" | tee -a "$LDIR/run.log"
  else
    n=$(wc -l < "$out/merged_loops.bedpe" 2>/dev/null || echo 0)
    echo "[$(date '+%F %T')] $base  done -> $n loops" | tee -a "$LDIR/run.log"
  fi
done

echo "[$(date '+%F %T')] ALL DONE" | tee -a "$LDIR/run.log"
