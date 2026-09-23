#!/bin/bash
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks
#SBATCH -c 18 #Number of threads per process
#SBATCH -t 07:00:00 #Short for --time walltimelimit
#SBATCH -e /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/spaceranger_liver_rerun.sh.e
#SBATCH -o /tscc/nfs/home/cmiciano/job_outs/Liver/RNA/spatial/spaceranger_liver_rerun.sh.o
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd772 #Allocation name
#SBATCH --mail-type ALL #Optional, Send mail when job ends
#SBATCH --mail-user cnmiciano@health.ucsd.edu #Optional, Send mail to this address
#SBATCH --propagate=NONE
#SBATCH --job-name=sr_liv


# -t 48:00:00 #Short for --time walltimelimit
# -p gold #Partition name
# -q hcg-csd772 #QOS name

export PATH=/opt/spaceranger-3.0.0:$PATH

# sbatch: error: NodeNames=tscc-11-32 CPUs=56 match no Sockets, Sockets*CoresPerSocket or Sockets*CoresPerSocket*ThreadsPerCore. Resetting CPUs.
#cd /tscc/projects/ps-epigen/10x_output/Space_Ranger/outputs6/reruns
cd /tscc/projects/ps-epigen/users/cmiciano/Liver/RNA/spatial/outputs/spaceranger


# check that qy_2645 has extra sequences in november folder? added 4th run in november folder
# outputs6 QY_2645_1_2_3 # use new alignment, using 'QY_2645_HL_180809c' instead of 'QY_2645_HL_180809'
#/tscc/projects/ps-epigen/10x_output/Space_Ranger/spaceranger-3.0.0/spaceranger count --id=QY_2645_1_2_3_4 \
#--sample=QY_2645,QY_2645_2,QY_2645_3 \
#--transcriptome=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A \ 
#--fastqs=/tscc/projects/ps-epigen/10x_output/Space_Ranger/fastqs/HFHKTAFX7/HFHKTAFX7,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage.ucsd.edu/240829_LH00444_0187_A22KK7TLT3,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241002_LH00444_0211_B22FJCLLT4,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241122_LH00444_0243_A22J3NKLT4 \
#--probe-set=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
#--cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2645_HL_180809/CAVG10717_2024-08-15_19-59-56_2024-08-15_19-36-45_H1-GXJ7QTG_A1_hl19002.tif \
#--create-bam=true \
#--slide=H1-GXJ7QTG \
#--area=A1 \
#--image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2645_HL_180809/NORMAL_180809_2nd\ reg.tif \
#--loupe-alignment=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2645_HL_180809c/H1-GXJ7QTG-A1-fiducials-image-registration.json 


#/tscc/projects/ps-epigen/10x_output/Space_Ranger/spaceranger-3.0.0/spaceranger count --id=QY_2665_1_2_3 \
#--sample=QY_2665,QY_2665_2 \
#--transcriptome=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A  \
#--fastqs=/tscc/projects/ps-epigen/10x_output/Space_Ranger/fastqs/HFHJNAFX7/HFHJNAFX7,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241002_LH00444_0211_B22FJCLLT4,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241122_LH00444_0243_A22J3NKLT4 \
#--probe-set=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
#--cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2665_HL1610029/CAVG10717_2024-09-04_21-48-29_2024-09-04_21-25-24_H1-9BBQ7R9_D1_hl170058.tif \
#--create-bam=true --slide=H1-9BBQ7R9 --area=D1  \
#--image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2665_HL1610029/NORMAL_HL160029.tif

#/tscc/projects/ps-epigen/10x_output/Space_Ranger/spaceranger-3.0.0/spaceranger count --id=QY_2666_1_2_3 \
#--sample=QY_2666,QY_2666_2 \
#--transcriptome=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A \
#--fastqs=/tscc/projects/ps-epigen/10x_output/Space_Ranger/fastqs/HFHJNAFX7/HFHJNAFX7,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241002_LH00444_0211_B22FJCLLT4,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241122_LH00444_0243_A22J3NKLT4 \
#--probe-set=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
#--cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2666_HL170058/CAVG10717_2024-09-04_21-48-29_2024-09-04_21-25-24_H1-9BBQ7R9_A1_hl160029.tif \
#--create-bam=true --slide=H1-9BBQ7R9 --area=A1 \
#--image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2666_HL170058/MASH_HL170058.tif

# for new version use
#--cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2666_HL170058c/CAVG10717_2024-09-04_21-48-29_2024-09-04_21-25-24_H1-9BBQ7R9_A1_hl160029.tif \
#--loupe-alignment=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2666_HL170058c/H1-9BBQ7R9-A1-fiducials-image-registration.json
#--image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2666_HL170058c/MASH_HL170058\(RE\).flippedtif.modified.tif

#/tscc/projects/ps-epigen/10x_output/Space_Ranger/spaceranger-3.0.0/spaceranger count --id=QY_2673_1_2_3 \
#--sample=QY_2673,QY_2673_2 \
#--transcriptome=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A \
#--fastqs=/tscc/projects/ps-epigen/10x_output/Space_Ranger/fastqs/HF7MJAFX7/HF7MJAFX7,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241002_LH00444_0211_B22FJCLLT4,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241122_LH00444_0243_A22J3NKLT4 \
#--probe-set=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
#--cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2673_HL230324_c/CAVG10717_2024-09-12_21-09-04_2024-09-12_20-40-05_H1-PT3Z29V_A1_sample-hl230325.tif \
#--create-bam=true --slide=H1-PT3Z29V --area=A1 \
#--loupe-alignment=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2673_HL230324_c/H1-PT3Z29V-A1-fiducials-image-registration.json \
#--image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2673_HL230324_c/MASL_HL230324.flipped.Hyun.tif



/tscc/projects/ps-epigen/10x_output/Space_Ranger/spaceranger-3.0.0/spaceranger count --id=QY_2674_1_2_3 \
 --sample=QY_2674,QY_2674_2 \
 --transcriptome=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A \
 --fastqs=/tscc/projects/ps-epigen/10x_output/Space_Ranger/fastqs/HF7MJAFX7/HF7MJAFX7,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241002_LH00444_0211_B22FJCLLT4,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241122_LH00444_0243_A22J3NKLT4 \
 --probe-set=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
 --cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2674_HL230325/CAVG10717_2024-09-12_21-09-04_2024-09-12_20-40-05_H1-PT3Z29V_D1_sample-hl230324.tif \
 --create-bam=true --image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2674_HL230325/MASH_HL230325\ \(1\)_flipped.tif \
 --loupe-alignment=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2674_HL230325/H1-PT3Z29V-D1-fiducials-image-registration.json 

 #/tscc/projects/ps-epigen/10x_output/Space_Ranger/spaceranger-3.0.0/spaceranger count --id=JL_104_2_3 \
 #--sample=JL_104,JL_104_2 \
 #--transcriptome=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A \
 #--fastqs=/tscc/projects/ps-epigen/10x_output/Space_Ranger/fastqs/HF7MJAFX7/HF7MJAFX7,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241002_LH00444_0211_B22FJCLLT4,/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage2.ucsd.edu/241122_LH00444_0243_A22J3NKLT4 \
 #--probe-set=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
 #--cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/JL_104_HL221006/CAVG10717_2024-09-18_21-01-30_2024-09-18_20-30-47_H1-9XPFWYD_D1_SandwichClosed.tif \
 #--create-bam=true \
 #--image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/JL_104_HL221006/MASH_HL221006_flipped.tif \
 #--loupe-alignment=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/JL_104_HL221006/H1-9XPFWYD-D1-fiducials-image-registration.json 


/tscc/projects/ps-epigen/10x_output/Space_Ranger/spaceranger-3.0.0/spaceranger count --id=QY_2480_Liver_MASH --sample=QY_2480_Liver_MASH --transcriptome=/tscc/projects/ps-epigen/10x_output/Space_Ranger/refdata-gex-GRCh38-2020-A --fastqs=/tscc/projects/ps-epigen/10x_output/IGM-Storage/igm-storage.ucsd.edu/240605_LH00444_0124_A22K3T3LT3 --probe-set=/tscc/projects/ps-epigen/10x_output/Space_Ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv --cytaimage=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2480_HL_2201019/CAVG10717_2024-05-01_20-41-17_2024-05-01_20-16-08_H1-X2FTDX6_D1_220825.tif --create-bam=true --slide=H1-X2FTDX6 --area=D1 --image=/tscc/projects/ps-epigen/10x_output/Space_Ranger/inputs/QY_2480_HL_2201019/HL221019 (2).image J.tif

