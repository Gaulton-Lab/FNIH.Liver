#!/bin/bash
while getopts t:i:s:o:v:g: flag
do
    case "${flag}" in
        t) TMPdir=${OPTARG};;
        i) cellrangerOuts=${OPTARG};;
        s) sample=${OPTARG};;
        o) outDir=${OPTARG};;
        v) genotypeVCF=${OPTARG};;
        g) goodBarcodes=${OPTARG};;
    esac
done

# Set the temporary directory path specific to your script
#tmp_dir=  #"/nfs/lab/rlmelton/demuxlet_tmp"

# Create the temporary directory
mkdir -p "$TMPdir"

# Set the TMPDIR environment variable to the temporary directory path
export TMPDIR="$TMPdir"

# create demux output dir
#mkdir $outDir

echo "[[demuxlet script]]: making mod bams" `date`
# Make remove chr bams
samtools view -H ${cellrangerOuts}/atac_possorted_bam.bam | sed -e 's/SN:chr1/SN:1/' | sed -e 's/SN:chr2/SN:2/' | sed -e 's/SN:chr3/SN:3/' | sed -e 's/SN:chr4/SN:4/' | sed -e 's/SN:chr5/SN:5/' | sed -e 's/SN:chr6/SN:6/' | sed -e 's/SN:chr7/SN:7/' | sed -e 's/SN:chr8/SN:8/' | sed -e 's/SN:chr9/SN:9/' | sed -e 's/SN:chr10/SN:10/' | sed -e 's/SN:chr11/SN:11/' | sed -e 's/SN:chr12/SN:12/' | sed -e 's/SN:chr13/SN:13/' | sed -e 's/SN:chr14/SN:14/' | sed -e 's/SN:chr15/SN:15/' | sed -e 's/SN:chr16/SN:16/' | sed -e 's/SN:chr17/SN:17/' | sed -e 's/SN:chr18/SN:18/' | sed -e 's/SN:chr19/SN:19/' | sed -e 's/SN:chr20/SN:20/' | sed -e 's/SN:chr21/SN:21/' | sed -e 's/SN:chr22/SN:22/' | sed -e 's/SN:chrX/SN:X/' | sed -e 's/SN:chrY/SN:Y/' | sed -e 's/SN:chrM/SN:MT/' | samtools reheader - ${cellrangerOuts}/atac_possorted_bam.bam > ${outDir}/${sample}_atac_possorted_bam_remove_chr.bam

samtools view -H ${cellrangerOuts}/gex_possorted_bam.bam | sed -e 's/SN:chr1/SN:1/' | sed -e 's/SN:chr2/SN:2/' | sed -e 's/SN:chr3/SN:3/' | sed -e 's/SN:chr4/SN:4/' | sed -e 's/SN:chr5/SN:5/' | sed -e 's/SN:chr6/SN:6/' | sed -e 's/SN:chr7/SN:7/' | sed -e 's/SN:chr8/SN:8/' | sed -e 's/SN:chr9/SN:9/' | sed -e 's/SN:chr10/SN:10/' | sed -e 's/SN:chr11/SN:11/' | sed -e 's/SN:chr12/SN:12/' | sed -e 's/SN:chr13/SN:13/' | sed -e 's/SN:chr14/SN:14/' | sed -e 's/SN:chr15/SN:15/' | sed -e 's/SN:chr16/SN:16/' | sed -e 's/SN:chr17/SN:17/' | sed -e 's/SN:chr18/SN:18/' | sed -e 's/SN:chr19/SN:19/' | sed -e 's/SN:chr20/SN:20/' | sed -e 's/SN:chr21/SN:21/' | sed -e 's/SN:chr22/SN:22/' | sed -e 's/SN:chrX/SN:X/' | sed -e 's/SN:chrY/SN:Y/' | sed -e 's/SN:chrM/SN:MT/' | samtools reheader - ${cellrangerOuts}/gex_possorted_bam.bam > ${outDir}/${sample}_gex_possorted_bam_remove_chr.bam

ATACbam=${outDir}/${sample}_atac_possorted_bam_remove_chr.bam
RNAbam=${outDir}/${sample}_gex_possorted_bam_remove_chr.bam

echo "[[demuxlet script]]: made mod bams" `date`

# merge bams
echo "[[demuxlet script]]: merging bams" `date`
mergeBam="${outDir}/${sample}_merge.bam"

samtools merge $mergeBam $ATACbam $RNAbam
echo "[[demuxlet script]]: done merging bams" `date`

demuxOut="${outDir}/${sample}"

echo "[[demuxlet script]]: starting demuxlet" `date`

# Replace <demuxlet command> with the actual command you want to run
popscle demuxlet --sam $mergeBam --vcf $genotypeVCF --field GT --out $demuxOut --group-list $goodBarcodes

# Clean up and remove the temporary directory
rm -rf "$TMPdir"
echo "[[demuxlet script]]: done with demuxlet" `date`
echo "[[demuxlet script]]: removing un-merged remove_chr bams" `date`

rm $ATACbam
rm $RNAbam

echo "[[demuxlet script]]: done with full demuxlet script" `date`