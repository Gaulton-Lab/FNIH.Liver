#!/bin/bash

path=./reference/241022_Peaks/NMF

file=${path}/FNIH_Liver_pool.ATAC_celltype.LiverUnionPeaks.cpm.tsv
fname=`basename $file .cpm.tsv`
path2script=~/scripts/git/snATACutils/bin

function loadavg {
   while [ cat /proc/loadavg | awk '{print int($1)}' -gt 50 ]; do sleep 120; date; done;
}

cut -f 1 ${file} | sed '1d' > ${path}/${fname}.xgi
head -n 1 ${file} | tr '\t' '\n' > ${path}/${fname}.ygi
sed '1d' ${file} | cut -f 2- > ${path}/${fname}_label.tmp
cd ${path}/
echo "
import numpy as np
from scipy import sparse
from scipy.sparse import save_npz, load_npz

data = np.loadtxt('"${fname}_label.tmp"')
data_sp = sparse.csr_matrix(data)
npz_file = '"${fname}.npz"'
save_npz(npz_file, data_sp)
" | python


