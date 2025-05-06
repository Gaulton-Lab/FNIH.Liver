

# EFHD1 chrombpnet
screen -S chrombp
srun --pty --nodes=1 --ntasks-per-node=2 --cpus-per-task=4 --gpus 2 --mem=128G --account=csd854 -t 24:00:00 -p rtx6000 -q condo-gpu --wait 0 /bin/bash

conda activate chrombpnet.old
module load cpu
module load bcftools

bcftools view /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/tensorQTL/04_eQTLs/genotypes/Liver.QTL.r209.maf05.matchednames.vcf.gz -i 'CHROM=="chr2" && POS==232645051' -Oz -o /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/rs7604422.vcf.gz
cd /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/
tabix rs7604422.vcf.gz
bcftools consensus -f /tscc/projects/ps-gaultonlab/ref/ref_genome_GRCh38.p14/hg38.p14.fa /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/rs7604422.vcf.gz -o /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/rs7604422.fa

module load gpu
module load cuda/11.2
module load cudnn/8.1

chrombpnet contribs_bw -m /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/chrombpnet_nobias.h5 -r /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/250303_WE_MVP_hg38_credible_set_1k.narrowPeak -g /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/rs7604422.fa -c /tscc/projects/ps-gaultonlab/ref/hg38.chrom.sizes -op /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/contribs_bw_rs7604422_genome -pc {counts,profile}
chrombpnet contribs_bw -m /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/chrombpnet_nobias.h5 -r /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/250303_WE_MVP_hg38_credible_set_1k.narrowPeak -g /tscc/projects/ps-gaultonlab/ref/ref_genome_GRCh38.p14/hg38.p14.fa -c /tscc/projects/ps-gaultonlab/ref/hg38.chrom.sizes -op /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/contribs_bw_ref_genome -pc {counts,profile}
chrombpnet pred_bw -bm /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/bias_model_scaled.h5 -cm /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/chrombpnet.h5 -cmb /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/chrombpnet_nobias.h5 -r /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/250303_WE_MVP_hg38_credible_set_1k.narrowPeak -g /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/rs7604422.fa -c /tscc/projects/ps-gaultonlab/ref/hg38.chrom.sizes -op /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/preds_bw_rs7604422_genome 
chrombpnet pred_bw -bm /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/bias_model_scaled.h5 -cm /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/chrombpnet.h5 -cmb /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/Results.Hep.Disease/Control/models/chrombpnet_nobias.h5 -r /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/250303_WE_MVP_hg38_credible_set_1k.narrowPeak -g /tscc/projects/ps-gaultonlab/ref/ref_genome_GRCh38.p14/hg38.p14.fa -c /tscc/projects/ps-gaultonlab/ref/hg38.chrom.sizes -op /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/13_chromBPNet_Peaks/variant_prediction/chrombpnet_plotting/preds_bw_ref_genome 