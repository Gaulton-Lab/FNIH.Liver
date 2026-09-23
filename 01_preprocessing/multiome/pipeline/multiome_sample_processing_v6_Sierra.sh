#!/usr/bin/env bash:

### This is the overall wrapper script for the 10X Multiome data processing pipeline.
### This script was originally designed for the multiomic_islet project and thus has a few references paths
### hardcoded. When applying this script to other data make sure to change these (indicated in line).


### Read in dynamic inputs from flags using getopts
while getopts s:d:r:o:t:v:m:n:e:g: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
        d) sample_dir=${OPTARG};;
        r) reticulate_dir=${OPTARG};;
        o) overall_output_dir=${OPTARG};;
        t) TMPdir=${OPTARG};;
        v) genotypeVCF=${OPTARG};;
        m) min_RNA=${OPTARG};;
        n) min_ATAC=${OPTARG};;
        e) exp_donors=${OPTARG};;
        g) marker_gene_file=${OPTARG};;
        ## add in flag for number of samples in pool so we can check that against vcf
		## possibly add in marker gene list as arg
    esac
done

### HARDCODED PATHS AND VARIABLES -- CHANGE AS NECESSARY
# locations of the component scripts
filtering_Rscript_fp="/nfs/lab/rlmelton/scripts/multiome_pipeline_v4/github/1sample_metrics_filtering_v4.R"
ATAClfm_Pyscript_fp="/nfs/lab/rlmelton/scripts/multiome_pipeline_v4/github/2ATAC_processing_v4.py"
demuxlet_script_fp="/nfs/lab/rlmelton/scripts/multiome_pipeline_v4/github/Popscle_demuxletAllinOneScript_final_v1.sh"
soupX_Rscript_fp="/nfs/lab/rlmelton/scripts/multiome_pipeline_v4/github/3sample_postfilter_SoupX_ATACwindows_demuxlet_v5_fNIHmod.R"

# locations of the marker gene lists
marker_file=$marker_gene_file #"/nfs/lab/projects/multiomic_islet/references/islet_markers.txt"
short_marker_file=$marker_gene_file #"/nfs/lab/projects/multiomic_islet/references/islet_markers_shortlist.txt"

# ALSO NOTE THAT THE PYTHON INSTALL TO USE TO RUN THE PYTHON SCRIPT IS CURRENTLY HARDCODED BELOW

### Set up necessary file structure
# make the output directory for this sample
cd $overall_output_dir
outs_dir="$(pwd)/$sample"
if [ ! -d "$outs_dir" ]; then mkdir $outs_dir; fi

# make the plotting subfolder
cd $outs_dir
plots_dir="$(pwd)/plots"
if [ ! -d "$plots_dir" ]; then mkdir $plots_dir; fi

#create a log file
log_file="$outs_dir/log_file.txt"
if [ ! -d "$log_file" ]; then touch $log_file; fi

echo "Starting script "`date` >> $log_file

### Run the scripts in order, with necessary inputs
echo "Individual Sample Processing for Sample "$sample >> $log_file
echo "" >> $log_file

## checks vcfs against expected number of samples in pooled samples
echo "First check is the number of samples in vcf matches the number of expected donors!" >> "$log_file"
num_in_vcf=$(awk '{print NF}' "$genotypeVCF" | sort -nu | tail -n 1)

# Subtract 9 from num_in_vcf: there are 9 columns always in a vcf, then 1 column per donor
num_in_vcf=$((num_in_vcf - 9))

echo "$num_in_vcf" >> "$log_file"
echo "$exp_donors" >> "$log_file"

if [ "$num_in_vcf" -eq "$exp_donors" ]; then
    echo "The VCF has the expected number of donors" >> "$log_file"
else
    echo "The VCF is missing donors" >> "$log_file"
    exit 1
fi

# Run the metric filtering script with static thresholds
#echo "Starting Rscript #1 "`date` >> $log_file
echo "Starting R script #1: "`date` >> $log_file
Rscript $filtering_Rscript_fp $sample_dir $outs_dir $min_RNA $min_ATAC $reticulate_dir $marker_file >> $log_file
echo "" >> $log_file

# Make ATAC long format windows matrix script
cd $sample_dir
echo "Starting Pyscript #1: "`date` >> $log_file
### IMPORTANT NOTE: I've hardcoded in the python install on Ophelia we need to use bc my PATH defaults to
### a different install for some reason. If your system is setup well you can probably run the commented line instead
python $ATAClfm_Pyscript_fp --bam "atac_possorted_bam.bam" --keep $outs_dir"/filtered_barcodes.txt" --outdir $outs_dir"/" --fragments "atac_fragments.tsv.gz"  >> $log_file &
 #/usr/bin/python3 $ATAClfm_Pyscript_fp --bam "atac_possorted_bam.bam" --keep $outs_dir"/filtered_barcodes.txt" --outdir $outs_dir --fragments "atac_fragments.tsv.gz"  >> $log_file
echo "" >> $log_file

# good barcodes file
goodBarcodes="${outs_dir}/filtered_barcodes.txt"
# create demuxlet output directory
demux_out="${outs_dir}/demux_out"
if [ ! -d "$demux_out" ]; then mkdir $demux_out; fi

# Run demuxlet script at the same time as the atac windows script
echo "Starting Demuxlet script: "`date` >> $log_file
bash $demuxlet_script_fp -t $TMPdir -i $sample_dir -s $sample -o $demux_out -v $genotypeVCF -g $goodBarcodes  >> $log_file &
echo "" >> $log_file

# have wrapper wait unitl both background scripts are done before starting r script #2
wait 
echo "Both Pyscript #1 and demuxlet script are done "`date` >> $log_file

# Adjust rwx of the 3 output files of ^ (also gzip the lfm matrix file)
chmod -R ugo+rwx "$outs_dir/atac_fragments.filtered_barcode.tsv.gz"
chmod -R ugo+rwx "$outs_dir/atac_possorted_reads.filtered_barcode.bed.gz"
chmod -R ugo+rwx "$outs_dir/atac_possorted_bam.filt.rmdup.bam"
gzip "$outs_dir/atac.long_fmt.filtered_barcode.mtx"
chmod -R ugo+rwx "$outs_dir/atac.long_fmt.filtered_barcode.mtx.gz"

#remove bams created in the atac processing 
rm "$outs_dir/atac_possorted_bam.compiled.filt.bam"
rm "$outs_dir/atac_possorted_bam.filt.md.bam"

# SoupX and Windows script
cd $outs_dir
echo "Starting Rscript #2: "`date` >> $log_file
Rscript $soupX_Rscript_fp $sample_dir $outs_dir $reticulate_dir $marker_file $short_marker_file $short_marker_file $demux_out >> $log_file
echo "Done with all scripts: "`date` >> $log_file

echo "Ending: "`date` >> $log_file
