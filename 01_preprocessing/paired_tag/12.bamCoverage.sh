#!/bin/bash

### genearte bw with RPGC norm. Add more options later

usage() { 
cat <<EOF
Usage: 12.bamCoverage.sh [-h] [-i input.bam] [-o out_dir] [-g hg38 / mm10 [-m dna/rna]
Description: bam to bigwig RPGC norm
Options:    
    -h           Print help and exit
    -i           input bam 
    -o           output path
    -g           genome used (hg38 / mm10)
    -m           DNA or RNA (dna / rna)
EOF
    exit 1
}

smooth="F"
mode="dna"
while getopts ":i:o:g:m:h" flag
do
    case "${flag}" in
        i) f=${OPTARG};;
        o) out_path=${OPTARG};;
		g) genome=${OPTARG};;
		m) mode=${OPTARG};;
		h) usage;;
    esac
done

if [ -z "${f}" ] || [ -z "${out_path}" ] || [ -z "${genome}" ]; then
    usage
fi

if [[ "${genome}" == "mm10" ]]
then 
    bl=~/genome/mm10/mm10_CnR_blacklist.bed
    size=2652783500
elif [[ "${genome}" == "hg38" ]]
then
    bl=~/genome/hg38/hg38-blacklist.merge.bed
    size=2913022398
fi

fname=`basename $f .bam`
if [[ $mode != "rna" ]]
then
	bamCoverage -b ${f} -o ${out_path}/${fname}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize ${size} --extendReads 300 -p 2 -bl ${bl} --minFragmentLength 20
else
	bamCoverage -b ${f} -o ${out_path}/${fname}.bw --normalizeUsing RPKM
fi

