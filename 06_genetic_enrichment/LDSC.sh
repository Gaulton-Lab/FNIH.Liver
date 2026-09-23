# liftover peaks to hg19

for cell in $(cat /nfs/lab/tscc/welison/FNIH.Liver/LDSC/celltypes.txt); do 
		#cat /nfs/lab/tscc/welison/FNIH.Liver/peak-call-pipeline/celltype.peaks/${cell}.0.9dectile.noheader.bed | awk -F'\t' 'BEGIN {OFS = FS} {print $1, $2, $3}' > /nfs/lab/tscc/welison/FNIH.Liver/LDSC/debug/peaks/${cell}_peaks.narrowPeak.awk.86map.pipeline.bed
		liftOver /nfs/lab/tscc/welison/FNIH.Liver/peak-call-pipeline/celltype.peaks/${cell}.0.8dectile.noheader.bed /nfs/lab/liftOver/hg38ToHg19.over.chain.gz /nfs/lab/tscc/welison/FNIH.Liver/LDSC/debug/peaks/$cell.86pool.0.8dectile.hg19.pipeline.bed /nfs/lab/tscc/welison/FNIH.Liver/LDSC/debug/peaks/$cell.86pool.0.8dectile.hg19.pipeline.unmapped.bed
done
# Don't forget the background
liftOver /nfs/lab/tscc/welison/FNIH.Liver/peak-call-pipeline/celltype.peaks/all.merged.peaks.bed /nfs/lab/liftOver/hg38ToHg19.over.chain.gz /nfs/lab/tscc/welison/FNIH.Liver/LDSC/debug/peaks/Merge.86pool.hg19.pipeline.bed /nfs/lab/tscc/welison/FNIH.Liver/LDSC/debug/peaks/Merge.86pool.hg19.pipeline.unmapped.bed

# Make annotations
for annot in $(cat /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/celltypes.txt) Merge; do 
  for i in {1..22}; do
    srun --job-name=$annot.$i \
         --output=/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/slurm.out/$annot.$i.%j.out \
         --error=/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/slurm.out/$annot.$i.%j.err \
         --nodes=1 \
         --ntasks-per-node=1 \
         --cpus-per-task=10 \
         --mem=32G \
         --account=csd854 \
         -t 24:00:00 \
         -p condo \.2
         -q condo \
         python /tscc/nfs/home/welison/git.packages/ldsc/make_annot.py \
		    --bed-file /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/peaks/$annot.86pool.hg19.pipeline.bed \
		    --bimfile /tscc/projects/ps-gaultonlab/ref/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim \
		    --annot-file /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/annot/$annot.ophelia.86pool.hg19.pipeline.${i}.annot.gz &
  done
done
# Wait for all background jobs to complete
wait

# Run the LD Regression
for annot in $(cat /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/celltypes.txt) Merge; do 
  for i in {1..22}; do
    srun --job-name=$annot.$i \
         --output=/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/slurm.out/$annot.$i.%j.out \
         --error=/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/slurm.out/$annot.$i.%j.err \
         --nodes=1 \
         --ntasks-per-node=1 \
         --cpus-per-task=10 \
         --mem=32G \
         --account=csd854 \
         -t 24:00:00 \
         -p condo \
         -q condo \
         python /tscc/nfs/home/welison/git.packages/ldsc/ldsc.py \
         --print-snps /tscc/projects/ps-gaultonlab/ref/LDSC/1000G_EUR_Phase3_baseline_snps/hm.${i}.snp \
         --ld-wind-cm 1.0 \
         --out /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/annot/$annot.ophelia.86pool.hg19.pipeline.${i} \
         --bfile /tscc/projects/ps-gaultonlab/ref/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} \
         --thin-annot \
         --annot /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/annot/$annot.ophelia.86pool.hg19.pipeline.${i}.annot.gz \
         --l2 &
  done
  wait
done

# Partition the heritability
for annot in $(cat /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/celltypes.txt); do echo -e "$annot\t/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/annot/$annot.ophelia.86pool.hg19.pipeline.,/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/annot/Merge.ophelia.86pool.hg19.pipeline."; done > /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/opelia.86pool.hg19.pipeline.cts.txt

srun --job-name=86bam.all \
         --output=/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/slurm.out/86bam.all.%j.out \                                         
         --error=/tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/slurm.out/86bam.all.%j.err \                                          
         --nodes=1 \                                                                                                                            
         --ntasks-per-node=1 \                                                                                                                  
         --cpus-per-task=10 \                                                                                                                   
         --mem=32G \                                                                                                                            
         --account=csd854 \                                                                                                                     
         -t 24:00:00 \                                                                                                                          
         -p condo \                                                                                                                             
         -q condo \                                                                                                                             
         python /tscc/nfs/home/welison/git.packages/ldsc/ldsc.py \
				--h2-cts /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/trait/NAFLD.EUR.MVP.2021.sumstats.gz \
				--ref-ld-chr-cts /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/opelia.86pool.hg19.pipeline.cts.txt \
				--ref-ld-chr /tscc/projects/ps-gaultonlab/ref/LDSC/1000G_EUR_Phase3_baseline/baseline. \
				--out /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/NEW.ALT.GWAS.h2cts.ophelia.86pool.hg19.pipeline \
				--w-ld-chr /tscc/projects/ps-gaultonlab/ref/LDSC/weights_hm3_no_hla/weights. &     
              
wait   

srun --job-name=Cirrhosis.coheritability \
         --output=Cirrhosis.coheritability.%j.out \
         --error=Cirrhosis.coheritability.%j.err \
         --nodes=1 \
         --ntasks-per-node=1 \
         --cpus-per-task=10 \
         --mem=32G \
         --account=csd854 \
         -t 24:00:00 \
         -p condo \
         -q condo \
         python /tscc/nfs/home/welison/git.packages/ldsc/ldsc.py \
--rg /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/trait/NAFLD.EUR.MVP.2021.sumstats.gz,CHIRHEP_NAS_meta_out_filtered.na_filt.column_filt.ldsc.sumstats.gz \
--ref-ld-chr /tscc/projects/ps-gaultonlab/ref/LDSC/eur_w_ld_chr/ \
--w-ld-chr /tscc/projects/ps-gaultonlab/ref/LDSC/eur_w_ld_chr/ \
--out cALT.Cirrhosis.coheritability

srun --job-name=Cirrhosis.coheritability \
         --output=Cirrhosis.coheritability.%j.out \
         --error=Cirrhosis.coheritability.%j.err \
         --nodes=1 \
         --ntasks-per-node=1 \
         --cpus-per-task=10 \
         --mem=32G \
         --account=csd854 \
         -t 24:00:00 \
         -p condo \
         -q condo \
         python /tscc/nfs/home/welison/git.packages/ldsc/ldsc.py \
--rg /tscc/projects/ps-gaultonlab/welison/FNIH.Liver/LDSC/debug/trait/NAFLD.EUR.MVP.2021.sumstats.gz,finngen_R11_NAFLD_meta_out_filtered.na_filt.column_filt.ldsc.sumstats.gz \
--ref-ld-chr /tscc/projects/ps-gaultonlab/ref/LDSC/eur_w_ld_chr/ \
--w-ld-chr /tscc/projects/ps-gaultonlab/ref/LDSC/eur_w_ld_chr/ \
--out cALT.NAFLD.coheritability
