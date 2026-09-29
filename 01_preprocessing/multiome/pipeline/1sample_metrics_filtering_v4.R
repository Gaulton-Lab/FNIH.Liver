# In this script we read in the raw RNA and ATAC counts and use input metrics
# to filter BCs for cells. We then cluster and produce the normal multiome plots.
# Finally we output an rds file and list of barcodes for use in future steps.

# This is v2 of the pipeline and the following changes have been implemented:
# - The marker genes used in the dotplot are read in from a file


### Set up steps
#Read in necessary inputs
args = commandArgs(trailingOnly=TRUE)
sample_dir <- args[1]
output_dir <- args[2] #sample specific subdir of overall outputs dir
rna_min_features <- as.numeric(args[3])
atac_min_fragments <- as.numeric(args[4])
reticulate_path <- args[5]
marker_file <- args[6]

#Set reticulate python
#These are essential reticulate functions that allow us to use python packages in R
#you'll need a conda env with 'leidenalg' and 'pandas' installed to do this
#then route reticulate to the python installed in that conda env with the below functions
Sys.setenv(RETICULATE_PYTHON=sprintf("%s/bin/python",reticulate_path))
library(reticulate)
reticulate::use_python(sprintf("%s/bin/python",reticulate_path))
reticulate::use_condaenv(reticulate_path)

#Load necessary libraries
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
warnLevel <- getOption('warn')
options(warn = -1)


### Read in raw data
print(paste("[[R script 1]]:Reading in raw data:",Sys.time()))
wd <- sample_dir
inputdata.10x <- Read10X_h5(file.path(wd, 'raw_feature_bc_matrix.h5'))
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$`Peaks`
adata <- CreateSeuratObject(counts=rna_counts)
adata[['percent.mt']] <- PercentageFeatureSet(adata, pattern = '^MT-')
print(paste("[[R script 1]]:",length(colnames(adata[["RNA"]])),"total BCs"))


### Filter by RNA nFeatures
adata_sub <- subset(
  x = adata,
  subset = nFeature_RNA >= rna_min_features
)
print(paste("[[R script 1]]:",length(colnames(adata_sub[["RNA"]])),"total BCs after RNA threshold"))


### Add QC data, plot UMIs/BC, then filter by ATAC number fragments
#cut down qc and make adata subset
print(paste("[[R script 1]]: Reading in QC data:",Sys.time()))
qc <- read.table(file.path(wd, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
qc <- as.data.frame(qc)
rownames(qc) <- qc$gex_barcode
qc <- qc[Cells(adata_sub), 1:length(colnames(qc))]

#filter by ATAC number fragments
adata_sub <- AddMetaData(adata_sub, qc)
adata_sub_atac <- subset(
  x = adata_sub,
  subset = atac_fragments >= atac_min_fragments
)
print(paste("[[R script 1]]:",length(colnames(adata_sub_atac[["RNA"]])),"total BCs after ATAC threshold"))

#remove multiplets using QC data
adata_sub_multiplet <- subset(
  x = adata_sub_atac,
  subset = excluded_reason != 1 
)
print(paste("[[R script 1]]:",length(colnames(adata_sub_multiplet[["RNA"]])),"total BCs after multiplet removal"))
print(paste("[[R script 1]]: Metric filtering done:",Sys.time()))

#make the raw RNA and ATAC UMIs/BC plots
par(mfrow=c(1,2))
bc_summary <- read.table(file.path(wd, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
bc_summary["new_gex_umis_count"] = bc_summary["gex_exonic_umis"] + bc_summary["gex_intronic_umis"]
bc_summary <- bc_summary[order(bc_summary$new_gex_umis_count, decreasing=TRUE),]

pdf_path0 <- file.path(output_dir, 'plots/plot_d.pdf')
pdf(paste(pdf_path0))
plot(bc_summary$new_gex_umis_count, log='xy', type='l', main="RNA UMIs",
     xlab="BC Rank", ylab="Number RNA UMIs", col='dark green', lwd=2)

pdf_path1 <- file.path(output_dir, 'plots/plot_e.pdf')
pdf(paste(pdf_path1))
bc_summariy <- bc_summary[order(bc_summary$atac_fragments, decreasing=TRUE),]
plot(bc_summary$atac_fragments, log='xy', type='l', main="ATAC UMIs",
     xlab="BC Rank", ylab="Number ATAC UMIs", col='orange', lwd=2)
garbage <- dev.off()
invisible(gc())

### Add in ATAC data for the barcodes that passed both filters
print(paste("[[R script 1]]: Adding in ATAC data:",Sys.time()))
atac_counts <- atac_counts[,colnames(adata_sub_multiplet)]
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(':', '-'))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'hg38'

frag.file <- file.path(wd, 'atac_fragments.tsv.gz')
suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_counts, sep=c(':', '-'), genome='hg38', fragments=frag.file, min.cells=-1, min.features=-1, annotation=annotations))
adata_sub_multiplet[['ATAC']] <- chrom_assay
gc()

# Adding raw RNA assay to the RDS file
DefaultAssay(adata) <- 'RNA'
raw_rna <- GetAssayData(object = adata, slot = 'counts')
raw_rna_assay <- CreateAssayObject(counts = raw_rna)
adata[['RNA_raw']] <- raw_rna_assay

# Save intermediate RDS 
saveRDS(raw_rna_assay, file = file.path(output_dir,sprintf("raw_rna_assay_intermediate.rds")))

### Compute intermediate bulk metrics
print(paste("[[R script 1]]: Computing intermediate bulk metrics:",Sys.time()))
adata <- adata_sub_multiplet
print(paste("[[R script 1]]: ATAC median high-quality fragments per cell (atac_fragments): ",median(adata[[]][,'atac_fragments'])))
DefaultAssay(adata) <- 'ATAC'
adata <- TSSEnrichment(adata)
print(paste("[[R script 1]]: ATAC median TSSe:",median(adata[[]][,'TSS.enrichment'])))
print(paste("[[R script 1]]: RNA median genes per cell:",median(adata[[]][,'nFeature_RNA'])))


### Cluster the filtered object
print(paste("[[R script 1]]: Analysing and clustering filtered object:",Sys.time()))
#RNA analysis
DefaultAssay(adata) <- 'RNA'
suppressMessages(adata <- SCTransform(adata, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_'))

#ATAC analysis
#We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(adata) <- 'ATAC'
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff='q0')
adata <- RunSVD(adata)
adata <- RunUMAP(adata, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')

#Multimodal analysis
adata <- FindMultiModalNeighbors(adata, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
adata <- RunUMAP(adata, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
adata <- FindClusters(adata, graph.name='wsnn', algorithm=4, resolution = .5, verbose=FALSE, method="igraph")

#Add in log metadata values
adata$log_nCount_ATAC = log(adata$nCount_ATAC)
adata$log_nCount_SCT = log(adata$nCount_SCT)
adata$log_nFeature_ATAC = log(adata$nFeature_ATAC)
adata$log_nFeature_SCT = log(adata$nFeature_SCT)


### Plotting results (as pdfs)
print(paste("[[R script 1]]: Outputting Plots:",Sys.time()))

#3 UMAPs
p1 <- DimPlot(adata, reduction='umap.rna', group.by='seurat_clusters', label=TRUE, label.size=8, repel=TRUE) + ggtitle('RNA')
p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('RNA only')
p2 <- DimPlot(adata, reduction='umap.atac', group.by='seurat_clusters', label=TRUE, label.size=8, repel=TRUE) + ggtitle('ATAC')
p2 <- p2 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('ATAC only')
p3 <- DimPlot(adata, reduction='umap.wnn', group.by='seurat_clusters', label=TRUE, label.size=8, repel=TRUE) + ggtitle('WNN')
p3 <- p3 + xlab('UMAP 1')+ ylab('UMAP 2') + ggtitle('Combined')
figure <- ggarrange(p1 & NoLegend(), p2 & NoLegend(), p3 & NoLegend(), ncol = 3, nrow = 1)
figure <- annotate_figure(figure, top = text_grob("Intermediate Filtered UMAPS - R Script #1", color = "red", face = "bold", size = 25))
pdf_path2 <- file.path(output_dir, '/plots/plot_a.pdf')
ggsave(pdf_path2,figure,width=10, height=3, units="in", scale=3)

#Violin plots of metrics
p1 <- VlnPlot(adata, features='nCount_SCT', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata$nCount_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p2 <- VlnPlot(adata, features='nFeature_SCT', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata$nFeature_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p3 <- VlnPlot(adata, features='nCount_ATAC', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata$nCount_ATAC), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p4 <- VlnPlot(adata, features='nFeature_ATAC', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata$nFeature_ATAC), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p5 <- VlnPlot(adata, features='nCount_RNA', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata$nCount_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p6 <- VlnPlot(adata, features='nFeature_RNA', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata$nFeature_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
figure <- ggarrange(p1, p3, p5, p2, p4, p6, ncol = 3, nrow = 2,
                    common.legend = TRUE,legend="none")
figure
out_figure <- annotate_figure(figure, top = text_grob("Intermediate Filtered Metrics Violin Plots - R Script #1", color = "red", face = "bold", size = 25))
pdf_path3 <- file.path(output_dir, '/plots/plot_b.pdf')
ggsave(pdf_path3, out_figure, width=7, height=3, units="in", scale=3)

#Marker gene plot
marker.genes.long <- scan(marker_file,what="",sep="\n") #read in marker genes list from a file
unique_markers <- unique(marker.genes.long)
options(repr.plot.width=12, repr.plot.height=6)
p1 <- DotPlot(adata, assay='SCT', features=unique_markers, cluster.idents=TRUE)
p1 <- p1 + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab('') + ylab('')
p1
figure <- annotate_figure(p1, top = text_grob("DotPlot of Marker Genes - R Script #1", color = "red", face = "bold", size = 25))
pdf_path4 <- file.path(output_dir, '/plots/plot_c.pdf')
ggsave(pdf_path4,figure,width=3, height=2, units="in", scale=3)


### Save final Seurat object to an .rds file and output the list of keep BCs
print(paste("[[R script 1]]: Saving RDS and barcodes list:",Sys.time()))
saveRDS(adata, file = file.path(output_dir,'intermediate_filtered.rds'))
filtered_bcs <- colnames(adata[["RNA"]])
write(filtered_bcs, file=(file.path(output_dir,'filtered_barcodes.txt')),sep='\n')

print(paste("[[R script 1]]: Done:",Sys.time()))
