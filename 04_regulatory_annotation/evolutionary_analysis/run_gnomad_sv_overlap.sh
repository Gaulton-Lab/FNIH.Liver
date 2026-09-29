#!/usr/bin/env bash
# Overlap of cREs with gnomAD-SV v4 CNVs, stratified by allele frequency.
set -euo pipefail

WORK=./05.R/cnv_overlap
CRE_DIR=./05.R/atac_peaks
BG_BED=$CRE_DIR/All_sig_cor_DEfibrosis_atac_peaks.bed   # or your full ATAC peak universe
GENOME=$WORK/hg38.chrom.sizes                            # fetch from UCSC
mkdir -p "$WORK" && cd "$WORK"

# ---- 1. gnomAD-SV v4 sites VCF (hg38) ----
# Public release: https://gnomad.broadinstitute.org/downloads#v4-structural-variants
GNOMAD_VCF=gnomad.v4.1.sv.sites.vcf.gz
[ -f "$GNOMAD_VCF" ] || \
  wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/$GNOMAD_VCF{,.tbi}

# ---- 2. CNVs only, PASS, with AF ----
# SVTYPE in {DEL, DUP, CNV}; bcftools writes BED with CHROM/POS0/END + AF + SVTYPE
bcftools view -f PASS -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="CNV"' "$GNOMAD_VCF" \
  | bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\t%INFO/AF\n' \
  | awk 'BEGIN{OFS="\t"} $5!="."' \
  | sort -k1,1 -k2,2n > gnomad_sv_cnv.bed

# ---- 3. AF bins ----
awk 'BEGIN{OFS="\t"} {
  af=$5+0
  if      (af <  0.001) bin="singleton"
  else if (af <  0.01 ) bin="rare"
  else if (af <  0.05 ) bin="low"
  else                  bin="common"
  print $1,$2,$3,$4,$5,bin > "gnomad_sv_cnv."bin".bed"
}' gnomad_sv_cnv.bed

# ---- 4. Intersect each cRE set against each AF bin ----
echo -e "celltype\tdirection\tn_cre\taf_bin\tn_overlap\tfrac_overlap" > overlap_summary.tsv
for cre in "$CRE_DIR"/*_DEfibrosis_atac_peaks.bed; do
  base=$(basename "$cre" _DEfibrosis_atac_peaks.bed)
  # base like "Hepatocytes_upATAC" or "Hepatocytes"
  ct=${base%_*ATAC}; dir=${base##*_}
  [ "$ct" = "$base" ] && dir="all"
  n=$(wc -l < "$cre")
  for bin in singleton rare low common; do
    k=$(bedtools intersect -u -a "$cre" -b gnomad_sv_cnv.$bin.bed | wc -l)
    awk -v ct="$ct" -v d="$dir" -v n="$n" -v b="$bin" -v k="$k" \
      'BEGIN{OFS="\t"; printf "%s\t%s\t%d\t%s\t%d\t%.4f\n", ct,d,n,b,k,(n? k/n:0)}' \
      >> overlap_summary.tsv
  done
done

# ---- 5. Enrichment vs. shuffled background ----
# For each cRE set, shuffle within the ATAC peak universe (preserves chrom + size dist)
# and recompute overlap N times to get a null. Report empirical p and fold-enrichment.
N_SHUF=1000
for cre in "$CRE_DIR"/*_upATAC_DEfibrosis_atac_peaks.bed "$CRE_DIR"/*_downATAC_DEfibrosis_atac_peaks.bed; do
  base=$(basename "$cre" .bed)
  for bin in rare low common; do
    obs=$(bedtools intersect -u -a "$cre" -b gnomad_sv_cnv.$bin.bed | wc -l)
    for i in $(seq 1 $N_SHUF); do
      bedtools shuffle -i "$cre" -g "$GENOME" -incl "$BG_BED" -noOverlapping \
        | bedtools intersect -u -a - -b gnomad_sv_cnv.$bin.bed | wc -l
    done > null.$base.$bin.txt
    awk -v o="$obs" 'BEGIN{n=0;ge=0;s=0} {s+=$1; n++; if($1>=o) ge++}
      END{printf "%s\t%s\tobs=%d\tmean_null=%.2f\tFE=%.2f\tp=%.4g\n",
           "'$base'","'$bin'",o,s/n,(s/n? o/(s/n):0),(ge+1)/(n+1)}' null.$base.$bin.txt \
      >> enrichment.tsv
  done
done
