#!/usr/bin/env bash
# 2-stage Hepatocytes loop pipeline: Peakachu t=0.7 -> FitHiC2 FDR<0.01 filter.
# Adapted from /projects/ps-renlab2/y2xie/projects/BICAN/ref/hba_snm3c/peakachu/scripts/ (100kb)
# to 10kb resolution for FNIH samples.

source /home/y2xie/miniconda3/etc/profile.d/conda.sh
conda activate peakachu
set -eo pipefail

RES=10000
MIN_DIST=60000        # 60 kb lower bound (matches Peakachu --lower 6 bins at 10kb)
MAX_DIST=3000000      # 3 Mb upper bound (matches Peakachu --upper 300 bins at 10kb)
SCRIPTS=peakachu
ROOT=./04.matrices/raw/disease
CT_ROOT=./04.matrices/raw/celltype
F2_ROOT=$ROOT/loop/fithic2
INP=$F2_ROOT/inputs
OUT=$F2_ROOT/output
FLT=$F2_ROOT/filter
LDIR=$F2_ROOT/logs
mkdir -p "$INP" "$OUT" "$FLT" "$LDIR"

# sample | cool/mcool URI | peakachu_t0.7_bedpe   (using ORIGINAL scoring, not rescore_cislong)
SAMPLES=(
  "Hepatocytes|$CT_ROOT/Hepatocytes.mcool::/resolutions/$RES|$CT_ROOT/loop/loops_t0.7/Hepatocytes.loops.bedpe"
  "Hepatocytes_Control|$ROOT/Hepatocytes_Control_10000.cool|$ROOT/loop/loops_t0.7/Hepatocytes_Control.loops.bedpe"
  "Hepatocytes_MASH|$ROOT/Hepatocytes_MASH_10000.cool|$ROOT/loop/loops_t0.7/Hepatocytes_MASH.loops.bedpe"
  "Hepatocytes_MASL|$ROOT/Hepatocytes_MASL_10000.cool|$ROOT/loop/loops_t0.7/Hepatocytes_MASL.loops.bedpe"
  "Hepatocytes_MetALD|$ROOT/Hepatocytes_MetALD_10000.cool|$ROOT/loop/loops_t0.7/Hepatocytes_MetALD.loops.bedpe"
  "Hepatocytes_allMASH|$ROOT/Hepatocytes_allMASH_10000.cool|$ROOT/loop/loops_t0.7/Hepatocytes_allMASH.loops.bedpe"
  "Endothelial|$CT_ROOT/Endothelial.mcool::/resolutions/$RES|$CT_ROOT/loop/loops_t0.7/Endothelial.loops.bedpe"
  "Endothelial_Control|$ROOT/Endothelial_Control_10000.cool|$ROOT/loop/loops_t0.7/Endothelial_Control.loops.bedpe"
  "Endothelial_MASH|$ROOT/Endothelial_MASH_10000.cool|$ROOT/loop/loops_t0.7/Endothelial_MASH.loops.bedpe"
  "Endothelial_MASL|$ROOT/Endothelial_MASL_10000.cool|$ROOT/loop/loops_t0.7/Endothelial_MASL.loops.bedpe"
  "Endothelial_MetALD|$ROOT/Endothelial_MetALD_10000.cool|$ROOT/loop/loops_t0.7/Endothelial_MetALD.loops.bedpe"
  "HSC|$CT_ROOT/HSC.mcool::/resolutions/$RES|$CT_ROOT/loop/loops_t0.7/HSC.loops.bedpe"
  "HSC_Control|$ROOT/HSC_Control_10000.cool|$ROOT/loop/loops_t0.7/HSC_Control.loops.bedpe"
  "HSC_MASH|$ROOT/HSC_MASH_10000.cool|$ROOT/loop/loops_t0.7/HSC_MASH.loops.bedpe"
  "HSC_MASL|$ROOT/HSC_MASL_10000.cool|$ROOT/loop/loops_t0.7/HSC_MASL.loops.bedpe"
  "HSC_MetALD|$ROOT/HSC_MetALD_10000.cool|$ROOT/loop/loops_t0.7/HSC_MetALD.loops.bedpe"
  "Myeloid|$CT_ROOT/Myeloid.mcool::/resolutions/$RES|$CT_ROOT/loop/loops_t0.7/Myeloid.loops.bedpe"
  "Myeloid_Control|$ROOT/Myeloid_Control_10000.cool|$ROOT/loop/loops_t0.7/Myeloid_Control.loops.bedpe"
  "Myeloid_MASH|$ROOT/Myeloid_MASH_10000.cool|$ROOT/loop/loops_t0.7/Myeloid_MASH.loops.bedpe"
  "Myeloid_MASL|$ROOT/Myeloid_MASL_10000.cool|$ROOT/loop/loops_t0.7/Myeloid_MASL.loops.bedpe"
  "Myeloid_MetALD|$ROOT/Myeloid_MetALD_10000.cool|$ROOT/loop/loops_t0.7/Myeloid_MetALD.loops.bedpe"
)

echo "[$(date '+%F %T')] Hepatocytes FitHiC2 pipeline — ${#SAMPLES[@]} samples at ${RES}bp, max_dist=${MAX_DIST}" \
  | tee "$LDIR/run.log"

for entry in "${SAMPLES[@]}"; do
  IFS='|' read -r label uri peakachu_bedpe <<< "$entry"
  SDIR=$OUT/$label
  mkdir -p "$SDIR"
  echo "[$(date '+%F %T')] === $label ===" | tee -a "$LDIR/run.log"
  echo "  uri:           $uri"             | tee -a "$LDIR/run.log"
  echo "  peakachu(t0.7): $peakachu_bedpe ($(wc -l < $peakachu_bedpe 2>/dev/null) loops)" | tee -a "$LDIR/run.log"

  contact=$INP/${label}.contact.gz
  frag=$INP/${label}.fragments.gz
  bias=$INP/${label}.bias.gz

  # ---- 1) contact file: chrom1 mid1 chrom2 mid2 count (exclude diagonal, intra only) ----
  if [[ -s "$contact" ]]; then
    echo "  contact: already exists, skipping"             | tee -a "$LDIR/run.log"
  else
    echo "[$(date '+%F %T')] $label  building contact.gz"  | tee -a "$LDIR/run.log"
    cooler dump --join "$uri" \
      | awk -v OFS='\t' '$1==$4 && $2!=$5 {
            m1 = int((($2+$3)/2));
            m2 = int((($5+$6)/2));
            printf "%s\t%d\t%s\t%d\t%d\n", $1, m1, $4, m2, $7+0
          }' \
      | gzip - > "$contact"
  fi

  # ---- 2) fragments file: chrom 0 mid total_contacts 1 ----
  if [[ -s "$frag" ]]; then
    echo "  fragments: already exists, skipping"  | tee -a "$LDIR/run.log"
  else
    echo "[$(date '+%F %T')] $label  building fragments.gz" | tee -a "$LDIR/run.log"
    tmp=$INP/${label}.fragments.tmp
    python "$SCRIPTS/scripts/summarize_intra_contacts.py" --cool "$uri" --out "$tmp" \
      >> "$LDIR/${label}.log" 2>&1
    awk -v OFS='\t' '{
        m = int((($2+$3)/2));
        printf "%s\t%d\t%d\t%d\t%d\n", $1, 0, m, $4, 1
      }' "$tmp" | gzip - > "$frag"
    rm -f "$tmp"
  fi

  # ---- 3) bias file: chrom mid bias (1/weight, mean-normalized; NaN->-1) ----
  if [[ -s "$bias" ]]; then
    echo "  bias: already exists, skipping"  | tee -a "$LDIR/run.log"
  else
    echo "[$(date '+%F %T')] $label  building bias.gz"  | tee -a "$LDIR/run.log"
    tmp=$INP/${label}.bias.tmp
    python "$SCRIPTS/scripts/mcool_to_fithic2_bias.py" --cool "$uri" --out "$tmp" \
      >> "$LDIR/${label}.log" 2>&1
    awk -v OFS='\t' '{
        m = int((($2+$3)/2));
        printf "%s\t%d\t%g\n", $1, m, $4
      }' "$tmp" | gzip - > "$bias"
    rm -f "$tmp"
  fi

  # ---- 4) run fithic2 ----
  if [[ -s "$SDIR/FitHiC.spline_pass1.res${RES}.significances.txt.gz" ]]; then
    echo "  fithic: already ran, skipping"  | tee -a "$LDIR/run.log"
  else
    echo "[$(date '+%F %T')] $label  fithic2 running..."  | tee -a "$LDIR/run.log"
    fithic \
      -i "$contact" \
      -f "$frag" \
      -t "$bias" \
      -o "$SDIR" \
      -r $RES \
      -L $MIN_DIST \
      -U $MAX_DIST \
      >> "$LDIR/${label}.log" 2>&1
  fi

  # ---- 5) filter Peakachu loops by FitHiC2 FDR<0.01 ----
  sig_bedpe=$SDIR/FitHiC.spline_pass1.res${RES}.significances001.bedpe
  filt=$FLT/${label}.loops.bedpe
  if [[ -s "$sig_bedpe" ]]; then
    n_fithic=$(wc -l < "$sig_bedpe")
    echo "[$(date '+%F %T')] $label  fithic FDR<0.01 sig pairs: $n_fithic"  | tee -a "$LDIR/run.log"
    cut -f1-6 "$sig_bedpe" \
      | bedtools pairtopair -a - -b "$peakachu_bedpe" \
      | cut -f7- - \
      | sort -k1,1 -k2,2n -u > "$filt"
    n_kept=$(wc -l < "$filt")
    n_peak=$(wc -l < "$peakachu_bedpe")
    echo "[$(date '+%F %T')] $label  peakachu(t0.7)=$n_peak  -> after fithic2 filter: $n_kept ($(awk "BEGIN{printf \"%.1f\", $n_kept*100/$n_peak}")%)" \
      | tee -a "$LDIR/run.log"
  else
    echo "[$(date '+%F %T')] $label  no FitHiC significances bedpe found -- skipping filter step"  | tee -a "$LDIR/run.log"
  fi
done

echo "[$(date '+%F %T')] ALL DONE"  | tee -a "$LDIR/run.log"
