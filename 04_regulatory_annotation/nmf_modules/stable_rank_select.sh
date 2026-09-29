#!/bin/bash

path="./reference/241022_Peaks/NMF"
path2script=~/scripts/git/snATACutils/bin
npz=FNIH_Liver_pool.ATAC_celltype.LiverUnionPeaks.npz
fname=`basename $npz .npz`

cd ${path}/res 

for i in cellSparse entropy; do echo $i >> sta.lite.header; done
cut -f 1 ${fname}.r10n10.box.sta | sed "1 s/contributes/rank/g" > sta.box.lite.header

for r in `seq 3 10`;
do 
	echo $r;
	n=10
	prefix=${fname}.r${r}n${n}
	# calculate cell sparseness and entropy using the statH file
	Rscript ${path2script}/nmfATAC.statBox.R -i ${fname}.r${r}n10.statH -o ${fname}.r${r}n10 >> ${fname}.r${r}n10.sta.txt
	cat ${prefix}.sta.txt | sed -e "s/^/${fname}\t${r}\t/g" | paste - sta.lite.header > ${prefix}.sta.tmp
	sed '1d' ${prefix}.box.sta | cut -f 2 | sed -e "1i ${r}" > ${prefix}.box.contributes
	sed '1d' ${prefix}.box.sta | cut -f 3 | sed -e "1i ${r}" > ${prefix}.box.sparseness
	sed '1d' ${prefix}.box.sta | cut -f 4 | sed -e "1i ${r}" > ${prefix}.box.entropy
	sed '1d' ${prefix}.statH | sed -e "s/^/${fname}\t${r}\t/g" > ${prefix}.statH.tmp
done;


cat *.sta.tmp | sed -e "1i samples\tranks\tval\tstat" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$4,$3}' > ${fname}.sta.txt
rm -rf *.sta.tmp
paste sta.box.lite.header *.box.contributes > ${fname}.contributes.sta.txt
paste sta.box.lite.header *.box.sparseness > ${fname}.sparseness.sta.txt
paste sta.box.lite.header *.box.entropy > ${fname}.entropy.sta.txt
rm -rf *.box.contributes *.box.sparseness *.box.entropy
cat *.statH.tmp | sed -e "1i samples\tranks\txgi\tindex\tclass0\tclass1\tcontributes\tsparseness\tentropy" > ${fname}.statH.sta.txt
rm -rf *.statH.tmp
Rscript ${path2script}/nmfATAC.plotBox.R -i ${fname}.statH.sta.txt -o ${fname}.statH.sta




