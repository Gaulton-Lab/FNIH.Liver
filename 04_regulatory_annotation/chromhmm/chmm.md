
### split fragment file into bed:
extract10X_from_fragment.py --cluster fragment.cluster_info.txt --indir fragment_10X/ --outprfx bed/FNIH_Liver --prefix TRUE

### create meta:
cd bed; for f in *bed; do fname=`basename $f .bed`; echo ${fname} | awk -v file=${f} -F'[_.]' '{print $3"\t"$4"\t"file}' - >> cellmarkfiletable; done

### chromhmm binarize
current="./"
indir=${current}/05.R/chromhmm/bed
meta=${current}/05.R/chromhmm/cellmarkfiletable
outdir=${current}/05.R/chromhmm/binary
java -mx12800M -jar ~/packages/ChromHMM/ChromHMM.jar BinarizeBed -b 200 -gzip -center ~/packages/ChromHMM/CHROMSIZES/hg38.txt ${indir} ${meta} ${outdir}

### chromhmm learnmodel
for r in {5..8}
do         
	mkdir -p ${current}/output/bin1k_${r}        
	nohup java -mx9600M -jar ~/packages/ChromHMM/ChromHMM.jar LearnModel -b 200 -gzip -r 1000 -p 8 -color 0,0,255 -l ~/packages/ChromHMM/CHROMSIZES/hg38.txt ${outdir} ${current}/05.R/chromhmm/output/bin1k_${r}/ ${r} hg38 & 
done

### make segment:
java -mx9600M -jar ~/packages/ChromHMM/ChromHMM.jar MakeSegmentation -b 200 -l ~/packages/ChromHMM/CHROMSIZES/hg38.txt output/bin200_5/model_5.txt qry_binary qry_output

### annotate peaks
for f in ../output/bin200_5/*_5_segments.bed.gz; do fname=`basename $f _5_segments.bed.gz`; bedtools intersect -wao -f 0.499 -a ../../../reference/241022_Peaks/consensus/LiverUnionPeaks.bed -b $f |  awk '{print $1"\t"$2"\t"$3"\t"$14}' - > Liver_UnionPeak_${fname}_annotate.bed; done
for f in disease_output/*bed; do fname=`basename $f ~/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/regions/tss_1500bp.bed -b ${f} | awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$11}' - > disease_tss/${fname}_annotate.bed; done
