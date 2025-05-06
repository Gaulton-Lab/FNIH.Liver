# Genotypes - TopMed v3

### Array Data

Compress array data with 7zip to transfer to machine with GenomeStudios

```bash
cp /nfs/lab/projects/nash_nafld_liver/genotypes/array/pool.1
/nfs/lab/welison/packages/7zz a 207062880027.7z 207062880027
/nfs/lab/welison/packages/7zz a 207550350003.7z 207550350003
cd /nfs/lab/projects/nash_nafld_liver/genotypes/array/pool.2
/nfs/lab/welison/packages/7zz a 207700400017.7z 207700400017
cd /nfs/lab/projects/nash_nafld_liver/genotypes/array/pool.3.4
/nfs/lab/welison/packages/7zz a 207710410015.7z 207710410015
/nfs/lab/welison/packages/7zz a 207700400013.7z 207700400013
/nfs/lab/welison/packages/7zz a 207710410009.7z 207710410009
```

Download via WinSCP

Open GenomeStudios

1. Make a new project
2. Select "Use sample sheet to load sample intensities"
3. Provide:
    1. Sample sheet-Provided by IGM in the folder you downloaded
    2. Data directory-The folder in the same directory as the sample sheet (sometimes there are multiple and you will have to run this for each folder)
    3. Manifest file-The folder containing the CNV manifest file from illumina  (The sample sheet has the array used, download from illumina) **THE ONE ON THE ILLUMINA SITE IS hg19, THE LAB COMPUTER HAS hg38, OR CONTACT ILLUMINA. ~~ALSO NEED hg38 CLUSTER FILE~~**
4. Browse for your cluster file and Check the boxes 
    1. “import cluster positions from cluster file” → cluster file can be downloaded form illumina
    2. Calculate sample and SNP statistics
5. Finish —> runs program
6. Navigate to Paired Sample Table and click “Export displayed data to a file”
7. Find the Report Wizard button and click “Custom Report” and navigate to the PLINK Input Report Plugin
******SELECT FORWARD STRAND******
8. Transfer back using WinSCP (did some renaming to make sure files match the array)

### Prep for TopMed

1. Merge plink files
    
    ```bash
    cd /nfs/lab/projects/nash_nafld_liver/genotypes/raw_data
    # Repeat for each array
    echo /nfs/lab/projects/nash_nafld_liver/genotypes/array/pool.3.4/PLINK_207710410015_hg38/pool.3.4.ped /nfs/lab/projects/nash_nafld_liver/genotypes/array/pool.3.4/PLINK_207710410015_hg38/pool.3.4.map >> merge_list.txt
    plink1.9 --merge-list merge_list.txt --make-bed --out FNIH.Liver.Raw.Genotypes
    ```
    
2. Prep for McCarthy tools
    
    ```bash
    plink1.9 --bfile FNIH.Liver.Raw.Genotypes --geno 0.05 --maf 0.01 --hwe 1e-5 -keep-allele-order --make-bed --out FNIH.Liver.Raw.Genotypes.filt
    plink1.9 --bfile FNIH.Liver.Raw.Genotypes.filt --recode vcf --keep-allele-order --out FNIH.Liver.Raw.Genotypes.filt
    plink1.9 --freq --bfile FNIH.Liver.Raw.Genotypes.filt --keep-allele-order --out FNIH.Liver.Raw.Genotypes.filt
    ```
    
3. Run McCarthy Tools
    
    ```bash
    perl /nfs/lab/welison/packages/HRC-1000G-check-bim.pl -b FNIH.Liver.Raw.Genotypes.filt.bim -f FNIH.Liver.Raw.Genotypes.filt.frq -r /nfs/lab/cmcgrail/lung_230615/all_lung_genotypes_merged/PASS.Variantsbravo-dbsnp-all.tab.gz -h
    sed s/plink/plink1.9/g Run-plink.sh > New_Run-plink.sh
    bash New_Run-plink.sh
    ```
    
4. Convert to bgzipped vcf
    
    ```bash
    for chr in {1..23}; do plink1.9 --bfile FNIH.Liver.Raw.Genotypes.filt-updated-chr$chr --recode vcf bgz --keep-allele-order --output-chr chrM --out FNIH.Liver.Raw.Genotypes.filt-updated-chr$chr; done
    ```
    
5. Download via WinSCP
6. Create TOPMed account and login
7. Create a new job:

- Proper Array Build (Currently we’re using Hg38)
- rsq filter- off
- Phasing: Eagle v2.4
- Population: Skip
- Mode: Quality control and imputation

TOPMed will send you an email when the files are ready. To get them on the server do the following:  
- Follow the link in your email and find the curl command under imputation results and link for downloading all files



### Merge, filter, and format

```bash
# Unzip and check files
# Replace qvMG1IL@Gq0gj with the password in the topmed email
for file in *.zip; do unzip -P qvMG1IL@Gq0gj $file; done
md5sum *.zip
# Compare to results.zip
```

```bash
# Filter by R^2 and minor allele frequency (MAF)
for chr in {1..22}; do
bcftools view chr${chr}.dose.vcf.gz -i 'R2>0.9 & MAF>0.01' -Oz -o chr${chr}.dose.r209.maf01.vcf.gz
done
#chr=X
#bcftools view chr${chr}.dose.vcf.gz -i 'R2>0.9 & MAF>0.01' -Oz -o chr${chr}.dose.r209.maf01.vcf.gz
```

```bash
# Merge chromosomes
bcftools concat chr*.dose.r209.maf01.vcf.gz -Oz -o merged.dose.r209.maf01.vcf.gz
# Change labels from rsID to position
bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' merged.dose.r209.maf01.vcf.gz -Oz -o FNIH.Liver.89merged.dose.r209.maf01.vcf.gz
# Rename samples
bcftools reheader FNIH.Liver.89merged.dose.r209.maf01.vcf.gz -s rename.merged.dose.r209.maf01.vcf.gz.txt -Oz -o FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz 
tabix FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz
```

### For Demuxlet (change based on modality, regions, and samples)

From here we need to

1. Select regions and samples using bedtools
2. Format (sort) to match demuxlet expected input

```bash
# This is the command with the region filter (-R, a bed file), 
# but for scHi-C we are using all regions currently
# For PairedTag we will subset by regions
# bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -R /nfs/lab/projects/nash_nafld_liver/paired-tag/H3K27ac/assets/paired-tag_Pool2/LV003_LV005_peaks.broadPeak -S samples.txt -Oz -o Liver_Paired-Tag_Pool2_16samples_H3K27ac_Peaks.vcf.gz
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S HiC.Pool1.samples.txt -Ov -o Liver_HiC_Pool1_All_snps.vcf
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S HiC.Pool2.samples.txt -Ov -o Liver_HiC_Pool2_All_snps.vcf
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S HiC.Pool3.samples.txt -Ov -o Liver_HiC_Pool3_All_snps.vcf
# HiC Pools 1-3 are the same donors as Paired-Tag Pools 1-3, for 4 HiC is missing one donor.
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S HiC.Pool4.samples.txt -Ov -o Liver_HiC_Pool4_All_snps.vcf
```

```bash
# This is the command with the region filter (-R, a bed file), 
# but for scHi-C we are using all regions currently
# For PairedTag we will subset by regions
# bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -R /nfs/lab/projects/nash_nafld_liver/paired-tag/H3K27ac/assets/paired-tag_Pool2/LV003_LV005_peaks.broadPeak -S samples.txt -Oz -o Liver_Paired-Tag_Pool2_16samples_H3K27ac_Peaks.vcf.gz
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S PairedTag.Pool1.samples.txt -Ov -o Liver_PairedTag_Pool1_All_snps.vcf
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S PairedTag.Pool2.samples.txt -Ov -o Liver_PairedTag_Pool2_All_snps.vcf
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S PairedTag.Pool3.samples.txt -Ov -o Liver_PairedTag_Pool3_All_snps.vcf
# HiC Pools 1-3 are the same donors as Paired-Tag Pools 1-3, for 4 HiC is missing one donor.
bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -S PairedTag.Pool4.samples.txt -Ov -o Liver_PairedTag_Pool4_All_snps.vcf
```

2.

```bash
# For HiC
#vcfs=("Liver_HiC_Pool1_All_snps" "Liver_HiC_Pool2_All_snps" "Liver_HiC_Pool3_All_snps" "Liver_HiC_Pool4_All_snps")
vcfs=("Liver_HiC_Pool1_All_snps" "Liver_HiC_Pool3_All_snps" "Liver_HiC_Pool4_All_snps")

for vcf in ${vcfs[@]};
do grep '^#' $vcf.vcf > $vcf.ordered.vcf && grep -v '^#' $vcf.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> $vcf.ordered.vcf; 
cat $vcf.ordered.vcf | awk '{print $1}' | uniq;
awk '{ if(match($0,/(##contig=<ID=)([1-9].*|MT.*|X.*|Y.*)/,m)) print m[1]"chr"m[2]; else print $0}' $vcf.ordered.vcf > $vcf.ordered.contig.vcf;
sed '/^##contig/d' $vcf.ordered.vcf > $vcf.ordered.nocontig.vcf;
done
```

```bash
# Change for PairedTag
vcfs=(" Liver_PairedTag_Pool1_All_snps" " Liver_PairedTag_Pool2_All_snps" " Liver_PairedTag_Pool3_All_snps" " Liver_PairedTag_Pool4_All_snps")

for vcf in ${vcfs[@]};
do grep '^#' $vcf.vcf > $vcf.ordered.vcf && grep -v '^#' $vcf.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> $vcf.ordered.vcf; 
cat $vcf.ordered.vcf | awk '{print $1}' | uniq;
sed '/^##contig/d' $vcf.ordered.vcf > $vcf.ordered.nocontig.vcf;
#awk '{ if(match($0,/(##contig=<ID=)([1-9].*|MT.*|X.*|Y.*)/,m)) print m[1]"chr"m[2]; else print $0}' $vcf.ordered.vcf > $vcf.ordered.contig.vcf;
#rm $vcf.ordered.vcf;
done
```

```bash
#All sample vcf

bcftools view FNIH.Liver.89merged.dose.r209.maf01.names.fixed.vcf.gz -Ov -o Liver_HiC_All_samples_All_snps.vcf

vcf=Liver_HiC_All_samples_All_snps
grep '^#' $vcf.vcf > $vcf.ordered.vcf && grep -v '^#' $vcf.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> $vcf.ordered.vcf
cat $vcf.ordered.vcf | awk '{print $1}' | uniq
sed '/^##contig/d' $vcf.ordered.vcf > $vcf.ordered.nocontig.vcf
```