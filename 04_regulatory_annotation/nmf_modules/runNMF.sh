#!/bin/bash
#SBATCH -J runNMF
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 1
#SBATCH -A csd788
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 12:00:00

module purge
module load cpu slurm gcc

path="./reference/241022_Peaks/NMF"
path2script=~/scripts/git/snATACutils/bin
npz=FNIH_Liver_pool.ATAC_celltype.LiverUnionPeaks.npz
fname=`basename $npz .npz`
cd ${path}
mkdir ${path}/res

~/anaconda3/envs/seurat/bin/python ${path2script}/nmfATAC.cluster2peak.lite.py -i ${npz} -x ${fname}.xgi -y ${fname}.ygi -r ${r} -n 10 -o res/${fname}.r${r}n10 &> res/${fname}.r${r}n10.log
~/anaconda3/envs/seurat/bin/Rscript ${path2script}/nmfATAC.plotH.R -i res/${fname}.r${r}n10.H.mx -o res/${fname}.r${r}n10
~/anaconda3/envs/seurat/bin/Rscript ${path2script}/nmfATAC.plotW.R -i res/${fname}.r${r}n10.W.mx -o res/${fname}.r${r}n10
~/anaconda3/envs/seurat/bin/python ${path2script}/nmfATAC.stat.py -m ${npz} -x ${fname}.xgi -y ${fname}.ygi --basis res/${fname}.r${r}n10.W.mx --coef res/${fname}.r${r}n10.H.mx -c 0.2 -o res/${fname}.r${r}n10
~/anaconda3/envs/seurat/bin/Rscript ${path2script}/nmfATAC.statBox.R -i res/${fname}.r${r}n10.statH -o res/${fname}.r${r}n10 >> res/${fname}.r${r}n10.sta.txt

