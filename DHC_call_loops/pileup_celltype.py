#!/usr/bin/env python
"""Cross-sample loop pileup: 4 cell-type aggregates × 4 cell-type cools."""
import os, sys, time, json
import numpy as np
import pandas as pd
import cooler
import cooltools
from coolpuppy.coolpup import pileup

ROOT = './04.matrices/raw/celltype'
LOOP_DIR = f'{ROOT}/loop/loops_t0.9'
OUT_DIR  = f'{ROOT}/loop/pileup_celltype'
os.makedirs(OUT_DIR, exist_ok=True)

CELLS = ['Endothelial', 'HSC', 'Hepatocytes', 'Myeloid']
RES = 10_000
FLANK = 100_000       # 100kb -> 21×21 pixel pileups at 10kb

t0 = time.time()
print('Loading coolers and computing cis-expected per sample…', flush=True)
clrs = {}
expected = {}
for c in CELLS:
    uri = f'{ROOT}/{c}.mcool::/resolutions/{RES}'
    clrs[c] = cooler.Cooler(uri)
    print(f'  {c}: cis-expected …', flush=True)
    expected[c] = cooltools.expected_cis(
        clrs[c],
        ignore_diags=2,
        chunksize=1_000_000,
        nproc=4,
    )
    expected[c].to_csv(f'{OUT_DIR}/expected.{c}.tsv', sep='\t', index=False)
print(f'  expected done in {time.time()-t0:.0f}s', flush=True)

print('Loading loops (t=0.9)…', flush=True)
loops = {}
for c in CELLS:
    fp = f'{LOOP_DIR}/{c}.loops.bedpe'
    df = pd.read_csv(fp, sep='\t', header=None,
                     names=['chrom1','start1','end1','chrom2','start2','end2','prob','raw'])
    loops[c] = df[['chrom1','start1','end1','chrom2','start2','end2']]
    print(f'  {c}: {len(df)} loops', flush=True)

def central_enrichment(arr):
    """Mean of 3×3 center / mean of four 3×3 corner patches."""
    n = arr.shape[0]
    c = n // 2
    center = arr[c-1:c+2, c-1:c+2]
    tl = arr[:3, :3]
    tr = arr[:3, -3:]
    bl = arr[-3:, :3]
    br = arr[-3:, -3:]
    corners = np.concatenate([tl.flatten(), tr.flatten(), bl.flatten(), br.flatten()])
    return float(np.nanmean(center) / np.nanmean(corners))

mat = pd.DataFrame(index=CELLS, columns=CELLS, dtype=float)
pile_store = {}
for src in CELLS:
    for tgt in CELLS:
        t1 = time.time()
        result = pileup(
            clr=clrs[tgt],
            features=loops[src],
            features_format='bedpe',
            expected_df=expected[tgt],
            flank=FLANK,
            nproc=4,
        )
        arr = result.iloc[0]['data']
        n_used = int(result.iloc[0]['n'])
        ce = central_enrichment(arr)
        mat.loc[src, tgt] = ce
        pile_store[f'{src}__{tgt}'] = arr
        print(f'  src={src:<13} tgt={tgt:<13} n={n_used:<6} central_enrich={ce:.3f}  ({time.time()-t1:.0f}s)', flush=True)

mat.to_csv(f'{OUT_DIR}/central_enrichment_4x4.tsv', sep='\t')
np.savez_compressed(f'{OUT_DIR}/pileups_arrays.npz', **pile_store)
print('\nN×N central enrichment (rows=loop source, cols=contact-map sample):')
print(mat.to_string(float_format='%.3f'))
print(f'\nWrote:\n  {OUT_DIR}/central_enrichment_4x4.tsv\n  {OUT_DIR}/pileups_arrays.npz')
print(f'Total time: {(time.time()-t0)/60:.1f} min')
