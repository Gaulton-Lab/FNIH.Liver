#!/usr/bin/env bash
# Extract PASS CNVs (DEL/DUP/CNV) from gnomAD-SV v4.1 VCF -> BED with AF + SVTYPE + SVLEN
set -euo pipefail
cd ./05.R/cnv_overlap

export PATH=/home/y2xie/miniconda3/envs/seurat/bin:$PATH

VCF=gnomad.v4.1.sv.sites.vcf.gz
OUT=gnomad_sv_cnv.bed

bcftools view -f PASS \
    -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="CNV"' \
    "$VCF" \
  | bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AF\n' \
  | awk 'BEGIN{OFS="\t"} $6!="." && $6!="" && $3!="."' \
  | sort -k1,1 -k2,2n > "$OUT"

wc -l "$OUT"

# Stratify by AF
awk 'BEGIN{OFS="\t"} {
  af=$6+0
  if      (af <  0.001) bin="singleton"
  else if (af <  0.01 ) bin="rare"
  else if (af <  0.05 ) bin="low"
  else                  bin="common"
  print $1,$2,$3,$4,$5,$6,bin > "gnomad_sv_cnv."bin".bed"
}' "$OUT"

for bin in singleton rare low common; do
  echo -n "$bin: "; wc -l gnomad_sv_cnv.$bin.bed
done
