#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='filter fragments by barcodes')
parser.add_argument('--cluster', type=str, dest="cluster", help='four columns meta: sample.fragments.gz, library, barcode, cluster')
parser.add_argument('--indir', type=str, dest="indir", help='directory path of input bedpe')
parser.add_argument('--outprfx', type=str, dest="outprfx", help='output directory')
parser.add_argument('--prefix', type=str, default="TRUE", dest="addprefix", help='Add prefix to barcode or not [TRUE/FALSE]; default is TRUE')

args = parser.parse_args()

import numpy as np
import pandas as pd
import gzip
import os
from time import perf_counter as pc

def run():
    start_time = pc()
    clusterf = args.cluster
    dirPrefix = str(args.indir)
    outPrefix = args.outprfx
    addPfx = args.addprefix
    print("filtering fragments by clusters...")
    generate_bedpe(clusterf, dirPrefix, outPrefix, addPfx)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_bedpe(clusterf, dirPrefix, outPrefix, addPfx):
    clust_dat = pd.read_csv(clusterf, sep="\t") ### four columns meta: sample.bam, library, barcode, cluster
    if addPfx == "TRUE":
        clust_dat["uniq_barcode"] = clust_dat[['library', 'barcode']].apply(lambda x: ':'.join(x), axis=1)
    else:
        clust_dat["uniq_barcode"] = clust_dat['barcode']

    clust_dat_dict = pd.Series(clust_dat["cluster"].values, index=clust_dat["uniq_barcode"]).to_dict()
    outf_dict = dict()
    sample_id = pd.unique(clust_dat["sample"])
    for cls_id in pd.unique(clust_dat["cluster"]):
        outbfname = "".join((outPrefix, "_", str(cls_id), ".bed"))
        outbf = open(outbfname, "a")
        outf_dict[cls_id] = outbf
    if len(sample_id) == 1:
        sample = sample_id[0]
        library_id = pd.unique(clust_dat.loc[clust_dat['sample'] == sample]['library'])
        library = library_id[0]
        generate_bedpe_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, addPfx)
    else:
        for sample in sample_id:
            library_id = pd.unique(clust_dat.loc[clust_dat['sample'] == sample]['library'])
            library = library_id[0]
            generate_bedpe_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, addPfx)
    [outf_dict[ofile].close() for ofile in outf_dict]

def generate_bedpe_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, addPfx):
    print("extract reads from", sample)
    iiput = "".join((dirPrefix, "/", sample))
    with oppf(iiput, mode='rb') as infile:
        for line in infile:
            try:
                dline = line.decode('utf8')
            except AttributeError:
                dline = line
            if dline.startswith('#'):
                continue
            else:
                ### chr1    9977    10132   GATCAGTAGGAACAAG-1      1
                fields = dline.strip("\n").split('\t')
                barcode = fields[3]
                if addPfx == "TRUE":
                    new_barcode = ":".join((library, barcode))
                else:
                    new_barcode = barcode
                    
                if new_barcode in clust_dat_dict:
                    fields[3] = new_barcode
                    wline = '\t'.join((fields))
                    cls = clust_dat_dict[new_barcode]
                    opair = outf_dict[cls]
                    opair.write(wline + "\n")
                else:
                    continue

    infile.close()

def oppf(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

if __name__ == "__main__":
    run()


