#!/usr/bin/env bash

source /home/y2xie/miniconda3/etc/profile.d/conda.sh
conda activate peakachu
set -eo pipefail

ROOT=/mnt/tscc2/y2xie/projects/77.LC/92.FNIH_DHC_IGM_241121/04.matrices/raw/celltype
MDIR=$ROOT/loop/models
ODIR=$ROOT/loop
LDIR=$ROOT/loop/logs
mkdir -p "$ODIR/scores" "$ODIR/loops" "$LDIR"

declare -A MODEL=(
  [B]=5
  [Cholangiocyte]=5
  [Endothelial]=50
  [Hepatocytes]=100
  [HSC]=50
  [Mast]=5
  [Myeloid]=10
  [NK]=5
  [Schwann]=5
  [T]=5
)

for ct in B Cholangiocyte Endothelial Hepatocytes HSC Mast Myeloid NK Schwann T; do
  m=${MODEL[$ct]}
  mcool=$ROOT/${ct}.mcool::/resolutions/10000
  model=$MDIR/high-confidence.${m}million.10kb.w6.pkl
  scored=$ODIR/scores/${ct}.scores.bedpe
  pooled=$ODIR/loops/${ct}.loops.bedpe
  log=$LDIR/${ct}.log

  echo "[$(date '+%F %T')] $ct  (model: ${m}M)" | tee -a "$LDIR/run.log"
  if [[ -s "$scored" ]]; then
    echo "  -> scores already exist, skipping score_genome" | tee -a "$LDIR/run.log"
  else
    set +e
    peakachu score_genome \
      -r 10000 \
      -p "$mcool" \
      --clr-weight-name weight \
      -m "$model" \
      -O "$scored" \
      > "$log" 2>&1
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
      echo "[$(date '+%F %T')] $ct  SCORE FAILED (rc=$rc) — see $log" | tee -a "$LDIR/run.log"
      rm -f "$scored"
      continue
    fi
  fi

  set +e
  peakachu pool \
    -r 10000 \
    -i "$scored" \
    -o "$pooled" \
    -t 0.9 \
    >> "$log" 2>&1
  set -e

  n=$(wc -l < "$pooled" 2>/dev/null || echo 0)
  echo "[$(date '+%F %T')] $ct  done -> $n loops" | tee -a "$LDIR/run.log"
done

echo "[$(date '+%F %T')] ALL DONE" | tee -a "$LDIR/run.log"
