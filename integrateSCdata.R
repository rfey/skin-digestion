# rmf 5.31.2023, last modified 12.19.2023

# script assumes the following directories exist:
# plots_analysis/   (for storing plots related to the analysis)
# integratedAnalysis/   (for storing data objects created during analysis)

#################
### FUNCTIONS ###
#################

readRSECfiles <- function(fileList_RSEC){
  myFiles <- readLines(con = fileList_RSEC)
  
  samples <- list()
  counts_container <- list()
  i <- 1
  for (f in myFiles){
    sample <- str_match(f, ".*/(.*?)_RSEC_MolsPerCell.csv")[[2]]
    print(sample)
    samples[[i]] <- sample    
    
    print("Reading RSEC file...")
    counts <- read.table(f, sep = ",", comment.char = "#", header = T, row.names = 1, check.names = FALSE)
    print("Transposing RSEC file...")
    counts_transposed <- as.data.frame(t(counts))
    counts_container[[i]] <- counts_transposed
    i <- i+1
  }
  names(counts_container) <- samples
  return(counts_container)
}

readSampleTagFiles <- function(fileList_SMK){
  myTags <- readLines(con = fileList_SMK)
  
  samples <- list()
  tags_container <- list()
  i <- 1
  for (f in myTags){
    sample <- str_match(f, ".*/(.*?)_Sample_Tag_Calls.csv")[[2]]
    print(sample)
    samples[[i]] <- sample
    
    tags <- read.table(f, sep = ",", comment.char = "#", header = T, row.names = 1)
    tags_container[[i]] <- tags
    i <- i+1
  }
  names(tags_container) <- samples
  return(tags_container)
}

plotViolinQC <- function(seurat_container, outbase){
  for (name in names(seurat_container)){
    x <- seurat_container[[name]]
    plist <- VlnPlot(x, features = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo"), 
                     combine = FALSE, pt.size = 0)
    for(i in 1:length(plist)) {
      print(i)
      plist[[i]] <- plist[[i]] + NoLegend() + 
        theme(axis.title = element_blank(), axis.title.x = element_blank(), 
              axis.text.x=element_blank()) +
        theme(plot.subtitle = element_text(size=8))
    }
    png(file = paste("plots_analysis/violinPlot_",outbase,"_",name,"_fresh.png", sep = ""), bg = "white")
    p <- plot_grid(plotlist = plist, ncol = 2) + plot_annotation(title = name)
    print(p)
    dev.off()

    pdf(file = paste("plots_analysis/violinPlot_",outbase,"_",name,"_fresh.pdf", sep = ""))
    p <- plot_grid(plotlist = plist, ncol = 3) + plot_annotation(title = name)
    print(p)
    dev.off()
  }
}

############
### MAIN ###
############

library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(clustree)
library(dittoSeq)
library(scProportionTest)
library(harmony)

# read input files
fileList_RSEC <- "fileList.RSEC_fresh_test.txt"  # text file with name of each RSEC file on a separate line
fileList_SMK <- "fileList.sampleTag_fresh_test.txt"  # text file with name of each sampleTag file on a separate line
metadataFile <- "metadata_fresh.txt"
patientMetadataFile <- "metadata_fresh_patientLevel.txt"

RSEC_container <- readRSECfiles(fileList_RSEC) # this takes a long time
names(RSEC_container)
names(RSEC_container) <- gsub("-", "_", names(RSEC_container))
names(RSEC_container)

#saveRDS(RSEC_container, file = "integratedAnalysis/data_RSECdataframes_fresh.RDS")
RSEC_container <- readRDS(file = "integratedAnalysis/data_RSECdataframes_fresh.RDS")

sampleTag_container <- readSampleTagFiles(fileList_SMK)
names(sampleTag_container)
names(sampleTag_container) <- gsub("-", "_", names(sampleTag_container))
names(sampleTag_container)

metadata <- read.table(metadataFile, comment.char = "")
colnames(metadata) <- gsub("#", "", metadata[1,])
metadata <- metadata[2:nrow(metadata),]
metadata$dataset <- gsub("-", "_", metadata$dataset)

patient_metadata <- read.table(patientMetadataFile, comment.char = "")
colnames(patient_metadata) <- patient_metadata[1,]
patient_metadata <- patient_metadata[2:nrow(patient_metadata),]
patient_metadata$Sample_Tag <- paste("SampleTag0", patient_metadata$SMK, "_hs", sep = "")
patient_metadata$Sample_Tag[patient_metadata$Sample_Tag == "SampleTag0Undetermined_hs"] <- "Undetermined"
patient_metadata$Sample_Tag[patient_metadata$Sample_Tag == "SampleTag0Multiplet_hs"] <- "Multiplet"

######################
### pre-processing ###
######################

# split out abseqs from experiments with abseq
abseq_container <- list()
rnaseq_container <- list()
abseq_names <- list()
i <- 1
for (RSEC in RSEC_container){
  # if statement explanation: grep for "|", returns a logical vector the same length as the number of row names
  # get unique values of the logical vector: if no abseqs, there will be one unique value (FALSE), 
  # if there are abseqs there will be two unique values (TRUE for abseq, FALSE for gene names)
  # check the number of unique values with length()
  # if there is more than one value, parse out abseqs
  if (length(unique(grepl("|", row.names(RSEC_container[[i]]), fixed = TRUE))) > 1) {
    name <- names(RSEC_container)[[i]]
    print(name)
    abseqs <- RSEC_container[[i]][grep("|", row.names(RSEC_container[[i]]), fixed = TRUE),]
    rna <- RSEC_container[[i]][which(grepl("|", row.names(RSEC_container[[i]]), fixed = TRUE) == FALSE),]
    j <- length(abseq_container) + 1
    print(j)
    abseq_container[[j]] <- abseqs
    abseq_names[[j]] <- name
  } else {
    rna <- RSEC_container[[i]]
  }
  rnaseq_container[[i]] <- rna
  i <- i + 1
}
names(rnaseq_container) <- names(RSEC_container)
names(abseq_container) <- abseq_names
rm(RSEC_container)

# create seurat objects with metadata
seurat_container <- list()
i <- 1
for (rna_df in rnaseq_container){
  name <- names(rnaseq_container)[[i]]

  meta <- metadata[(metadata$dataset == name),] # get metadata for this dataset
  meta_expanded <- meta[rep(1,ncol(rna_df)),] # expand df to number of cells for this dataset
  row.names(meta_expanded) <- colnames(rna_df) # set metadata row names as cell names from counts data
  
  print("Creating Seurat object...")
  obj <- CreateSeuratObject(counts = rna_df, project = name, meta.data = meta_expanded)
  
  # add abseq info if available-- not all samples have abseq info
  if (!(is.null(abseq_container[[name]]))) {
    print("Adding abseq information...")
    obj[["ADT"]] <- CreateAssayObject(counts = abseq_container[[name]]) 
  }
  
  ### add sample tags ###
  print("Adding sample tag information...")
  sample_tag_df <- sampleTag_container[[name]] # there are two cols, they are redundant so take only the first
  sample_tag_df$cellnames <- row.names(sample_tag_df)
  
  # get patient-level metadata for just this dataset
  patient_meta <- patient_metadata[patient_metadata$sample_group == name,]
  
  # merge patient-level metadata by sample tag
  sample_meta <- merge(sample_tag_df, patient_meta, by = "Sample_Tag")
  row.names(sample_meta) <- sample_meta$cellnames
  print(unique(sample_meta$Sample_Tag))

  obj <- AddMetaData(obj, metadata = sample_meta)
  
  # remove cells called Multiplet, and Undetermined if there are other actual SMKs
  # if Undetermined is the only one, keep it bc it means SMKs were not added in the experimental workflow
  # but the processing was run with the SMK reference file
  samples <- unique(obj@meta.data$Sample_Tag)
  if (length(samples) != 1){
    print(samples)
    samples <- samples[!(samples == "Undetermined" | samples == "Multiplet")]
    obj <- subset(x = obj, subset = Sample_Tag %in% samples)
  } else {
    print(samples)
    samples <- samples[!(samples == "Multiplet")]
    obj <- subset(x = obj, subset = Sample_Tag %in% samples)
  }
  
  # add %mitochondrial genes to metadata
  obj$percent_mito <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$percent_ribo <- PercentageFeatureSet(obj, pattern = "RP[L|S]")

  seurat_container[[i]] <- obj
  i <- i + 1
}
names(seurat_container) <- names(rnaseq_container)

#saveRDS(seurat_container, file = "integratedAnalysis/data_listOfSeuratObjects.RDS")
seurat_container <- readRDS("integratedAnalysis/data_listOfSeuratObjects.RDS")

#######################
### quality control ###
#######################

seurat_container <- readRDS("integratedAnalysis/data_listOfSeuratObjects.RDS")

# merge into one seurat object for QC plotting
obj <- merge(seurat_container[[1]], seurat_container[2:length(seurat_container)])

# plot the number of cells per dataset before filtering
bars <- c()
datasets <- c()
i <- 1
for (name in names(seurat_container)){
  datasets[[i]] <- name
  bars[[i]] <- nrow(seurat_container[[name]]@meta.data)
  i <- i + 1
}
df <- data.frame(unlist(datasets), unlist(bars))
colnames(df) <- c("datasets", "bars")
df$datasets <- gsub("-large-input","",df$datasets,fixed=TRUE)
ggplot(df, aes(x=datasets, y=bars)) + 
  geom_bar(stat = "identity") +
  ggtitle("Number of Cells per Dataset") +
  geom_text(aes(label=bars), vjust=1.5, color="cyan", size=3.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # rotate axis labels
ggsave(filename = "plots_analysis/barPlotCellCountPerDataset_fresh.pdf")
ggsave(filename = "plots_analysis/barPlotCellCountPerDataset_fresh.png", bg = "white")

# plot number of cells per sample tag for each dataset
for (name in names(seurat_container)) {
  bars <- c()
  tags <- c()
  x <- seurat_container[[name]]
  i <- 1
  for (tag in unique(x$Sample_Tag)){
    y <- subset(x, subset = Sample_Tag == tag)
    bars[[i]] <- ncol(y)
    tags[[i]] <- tag
    i <- i +1
  }
  df <- data.frame(unlist(tags), unlist(bars))
  colnames(df) <- c("tags","bars")
  
  ggplot(df, aes(x=tags, y=bars)) + 
    geom_bar(stat = "identity") +
    ggtitle(paste("Number of Cells per SMK",name,sep=": ")) +
    geom_text(aes(label=bars), vjust=1.5, color="cyan", size=3.5) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # rotate axis labels
  ggsave(filename = paste("plots_analysis/barPlotCellCountPerSMK_",name,"_fresh.pdf",sep=""))
  ggsave(filename = paste("plots_analysis/barPlotCellCountPerSMK_",name,"_fresh.png",sep=""), bg = "white")
}

# plot number of reads per cell in each dataset
plotViolinQC(seurat_container, "prefilter")

VlnPlot(object = obj, features = "nFeature_RNA", pt.size = 0) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = 200) +
  geom_text(aes(0.8,200,label = '200', vjust = -0.5)) +
  geom_hline(yintercept = 3000) +
  geom_text(aes(0.8,3000,label = '3000', vjust = -0.5)) +
  geom_hline(yintercept = 4000) +
  geom_text(aes(0.8,4000,label = '4000', vjust = -0.5)) +
  geom_hline(yintercept = 6000) +
  geom_text(aes(0.8,6000,label = '6000', vjust = -0.5))
ggsave("plots_analysis/violinPlot_prefilter_allDatasets_nFeature_fresh.png", height = 5, width = 8, bg = "white")
ggsave("plots_analysis/violinPlot_prefilter_allDatasets_nFeature_fresh.pdf")

VlnPlot(object = obj, features = "percent_mito", pt.size = 0) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = 20) +
  geom_text(aes(3.4,20,label = '20', vjust = -0.5)) +
  geom_hline(yintercept = 15) +
  geom_text(aes(3.4,15,label = '15', vjust = -0.5))
ggsave("plots_analysis/violinPlot_prefilter_allDatasets_percentMito_fresh.png", height = 5, width = 8, bg = "white")
ggsave("plots_analysis/violinPlot_prefilter_allDatasets_percentMito_fresh.pdf")

# by digest method-- CD45 pos
obj_pos <- merge(seurat_container[["DE1_CD45pos"]], seurat_container[["DE2_CD45pos"]])

unique(Idents(obj_pos))
Idents(obj_pos) <- "digest_method"
unique(Idents(obj_pos))
# explicitly set order of idents
obj_pos@active.ident <- factor(x = obj_pos@active.ident, levels = c("Sequential","Simultaneous"))
unique(Idents(obj_pos))
feats = c("nCount_RNA","nFeature_RNA","percent_mito")
titles = c("Number of\nMolecules","Number of\nUnique Genes","Percent\nMitochondrial Genes")
plist <- list()
for (i in seq(1,length(feats))){
  p <- VlnPlot(obj_pos, features = feats[[i]],
          pt.size = 0, cols = c("Sequential"="coral","Simultaneous"="turquoise")) + 
    NoLegend() + ggtitle(titles[[i]]) + xlab("") +
    theme(axis.text.x = element_blank())
  plist[[i]] <- p
}
legend <- get_legend(VlnPlot(obj_pos, features = feats[[i]],
                             cols = c("Sequential"="coral","Simultaneous"="turquoise")))
blank_plot <- ggplot() + theme_void()
plist[[4]] <- blank_plot
plist[[5]] <- legend
plist[[6]] <- blank_plot
plot_grid(plotlist = plist, nrow = 2, rel_heights = c(6,1))
ggsave("plots_analysis/violin_QC_digestMethod_prefilter_pos.pdf",
       units = "in", width = 8, height = 6)

# number of cells per digest method
obj_pos@meta.data %>% group_by(digest_method) %>% summarise(n())

# by digest method-- CD45 neg
obj_neg <- seurat_container[["DE2_CD45neg"]]

unique(Idents(obj_neg))
Idents(obj_neg) <- "digest_method"
unique(Idents(obj_neg))
# explicitly set order of idents
obj_neg@active.ident <- factor(x = obj_neg@active.ident, levels = c("Sequential","Simultaneous"))
unique(Idents(obj_neg))
feats = c("nCount_RNA","nFeature_RNA","percent_mito")
titles = c("Number of\nMolecules","Number of\nUnique Genes","Percent\nMitochondrial Genes")
plist <- list()
for (i in seq(1,length(feats))){
  p <- VlnPlot(obj_neg, features = feats[[i]],
               pt.size = 0, cols = c("Sequential"="coral","Simultaneous"="turquoise")) + 
    ggtitle(titles[[i]]) + xlab("") + NoLegend() +
    theme(axis.text.x = element_blank())
  plist[[i]] <- p
}
legend <- get_legend(VlnPlot(obj_neg, features = feats[[i]],
                             cols = c("Sequential"="coral","Simultaneous"="turquoise")))
blank_plot <- ggplot() + theme_void()
plist[[4]] <- blank_plot
plist[[5]] <- legend
plist[[6]] <- blank_plot
plot_grid(plotlist = plist, nrow = 2, rel_heights = c(6,1))
ggsave("plots_analysis/violin_QC_digestMethod_prefilter_neg.pdf", 
       units= "in", width = 8, height = 6)

# number of cells per digest method
obj_neg@meta.data %>% group_by(digest_method) %>% summarise(n())

### quality control: filter cells and reads ###

filtered_seurat <- list()
names(seurat_container) # DOUBLE CHECK YOUR ORDER
filtered_seurat[["DE1_CD45pos"]] <- subset(seurat_container[["DE1_CD45pos"]], subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mito < 25)
filtered_seurat[["DE2_CD45neg"]] <- subset(seurat_container[["DE2_CD45neg"]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 25)
filtered_seurat[["DE2_CD45pos"]] <- subset(seurat_container[["DE2_CD45pos"]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_mito < 25)
names(filtered_seurat) # check names

#saveRDS(filtered_seurat, file = "data_listOfSeuratObjects_filtered.RDS")

### post filtering plotting ###
# plot the number of cells per dataset after filtering
bars <- c()
datasets <- c()
i <- 1
for (name in names(filtered_seurat)){
  datasets[[i]] <- name
  bars[[i]] <- nrow(filtered_seurat[[name]]@meta.data)
  i <- i + 1
}
df <- data.frame(unlist(datasets), unlist(bars))
colnames(df) <- c("datasets", "bars")
df$datasets <- gsub("-large-input","",df$datasets,fixed=TRUE)
ggplot(df, aes(x=datasets, y=bars)) + 
  geom_bar(stat = "identity") +
  ggtitle("Number of Cells per Dataset, Filtered") +
  geom_text(aes(label=bars), vjust=1.5, color="cyan", size=3.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # rotate axis labels
ggsave(filename = "plots_analysis/barPlotCellCountPerDataset_fresh_filtered.pdf")
ggsave(filename = "plots_analysis/barPlotCellCountPerDataset_fresh_filtered.png", bg = "white")

# plot number of cells per sample tag for each dataset
for (name in names(filtered_seurat)) {
  bars <- c()
  tags <- c()
  x <- filtered_seurat[[name]]
  i <- 1
  for (tag in unique(x$Sample_Tag)){
    y <- subset(x, subset = Sample_Tag == tag)
    bars[[i]] <- ncol(y)
    tags[[i]] <- tag
    i <- i +1
  }
  df <- data.frame(unlist(tags), unlist(bars))
  colnames(df) <- c("tags","bars")
  ggplot(df, aes(x=tags, y=bars)) + 
    geom_bar(stat = "identity") +
    ggtitle(paste("Number of Cells per SMK ", name, ", filtered", sep = "")) +
    geom_text(aes(label=bars), vjust=1.5, color="cyan", size=3.5) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # rotate axis labels
  ggsave(filename = paste("plots_analysis/barPlotCellCountPerSMK_",name,"_fresh_filtered.pdf",sep=""))
  ggsave(filename = paste("plots_analysis/barPlotCellCountPerSMK_",name,"_fresh_filtered.png",sep=""), bg = "white")
}

# plot number of reads per cell in each dataset after filtering
obj_filtered <- merge(filtered_seurat[[1]], filtered_seurat[2:length(filtered_seurat)])

plotViolinQC(filtered_seurat, "postfilter")

VlnPlot(object = obj_filtered, features = "nFeature_RNA", pt.size = 0) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("plots_analysis/violinPlot_postfilter_allDatasets_nFeature_fresh.png", bg = "white")
ggsave("plots_analysis/violinPlot_postfilter_allDatasets_nFeature_fresh.pdf")

VlnPlot(object = obj_filtered, features = "percent_mito", pt.size = 0) +
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
ggsave("plots_analysis/violinPlot_postfilter_allDatasets_percentMito_fresh.png", bg = "white")
ggsave("plots_analysis/violinPlot_postfilter_allDatasets_percentMito_fresh.pdf")

##########################################
### CD45 pos: integration and analysis ###
##########################################

filtered_seurat <- readRDS(file = "integratedAnalysis/data_listOfSeuratObjects_filtered.RDS")

# get just CD45 positive datasets
pos <- list()
pos[["DE1_CD45pos"]] <- filtered_seurat[["DE1_CD45pos"]]
pos[["DE2_CD45pos"]] <- filtered_seurat[["DE2_CD45pos"]]

### RNA integration ###
# normalize and find variable features for each seurat object in the list
norm_list <- lapply(X = pos, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# standard integration
features <- SelectIntegrationFeatures(object.list = norm_list)
anchors <- FindIntegrationAnchors(object.list = norm_list, anchor.features = features)
pos <- IntegrateData(anchorset = anchors)

#saveRDS(pos, file = "integratedAnalysis/data_pos_integrated.RDS")
pos <- readRDS(file = "integratedAnalysis/data_pos_integrated.RDS")

### standard integrated analysis for CD45 positive ###
DefaultAssay(pos) <- "integrated"
pos <- ScaleData(pos) 
pos <- RunPCA(pos)
ElbowPlot(pos, ndims = 20)
pos <- RunUMAP(pos, dims = 1:10, reduction = "pca")
pos <- FindNeighbors(pos, reduction = "pca", dims = 1:10)
res_range <- c(seq(0.1,1.0,0.1), seq(1.2,2.0,0.2))
res_range
pos <- FindClusters(pos, resolution = res_range)

# use clustree to find best resolution
head(pos[[]])
clustree(pos, prefix = "integrated_snn_res.", show_axis = TRUE) +
  theme(legend.position="none")
ggsave("plots_analysis/clustree_prebatch_fresh_pos.png", bg = "white", width = 7, height = 12, units = "in")
ggsave("plots_analysis/clustree_prebatch_fresh_pos.pdf", width = 7, height = 12, units = "in")

# set optimal resolution
Idents(pos) <- pos$integrated_snn_res.0.2

#saveRDS(pos, file = "integratedAnalysis/data_pos_analyzed.RDS")
pos <- readRDS(file = "integratedAnalysis/data_pos_analyzed.RDS")

dittoDimPlot(pos, var = Idents(pos), reduction.use = "umap") +
  ggtitle("CD45+ Pre-batch correction")
ggsave("plots_analysis/dimPlot_seurat_clusters_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_seurat_clusters_fresh_pos.pdf")

dittoDimPlot(pos, var = Idents(pos), split.by = "dataset") +
  ggtitle("CD45+ Pre-batch correction: dataset")
ggsave("plots_analysis/dimPlot_seurat_clusters_splitDataset_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_seurat_clusters_splitDataset_fresh_pos.pdf")

dittoDimPlot(pos, var = "dataset", order = "increasing") +
  ggtitle("CD45+ Pre-batch correction: dataset, increasing")
ggsave("plots_analysis/dimPlot_dataset_incr_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_dataset_incr_fresh_pos.pdf")

dittoDimPlot(pos, var = "dataset", order = "decreasing") +
  ggtitle("CD45+ Pre-batch correction: dataset, decreasing")
ggsave("plots_analysis/dimPlot_dataset_decr_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_dataset_decr_fresh_pos.pdf")

# visualize other metadata to see if batch correction is needed
dittoDimPlot(pos, var = "digest_method") +
  ggtitle("CD45+ Pre-batch correction: digestion method")
ggsave("plots_analysis/dimPlot_digest_method_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_digest_method_fresh_pos.pdf")

dittoDimPlot(pos, var = "digest_method", split.by = "digest_method") +
  theme(legend.position = "none") +
  ggtitle("CD45+ Pre-batch correction: digestion method")
ggsave("plots_analysis/dimPlot_splitDigestMethod_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_splitDigestMethod_fresh_pos.pdf")

dittoDimPlot(pos, var = "patientID") +
  ggtitle("CD45+ Pre-batch correction: patient ID")
ggsave("plots_analysis/dimPlot_patientID_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_patientID_fresh_pos.pdf")

dittoDimPlot(pos, var = "patientID", split.by = "patientID") +
  ggtitle("CD45+ Pre-batch correction: patient ID") +
  theme(legend.position = "none")
ggsave("plots_analysis/dimPlot_splitPatientID_fresh_pos.png", bg = "white")
ggsave("plots_analysis/dimPlot_splitPatientID_fresh_pos.pdf")

# normalize ADT assay
DefaultAssay(pos) <- "ADT"
pos <- NormalizeData(pos, normalization.method = "CLR", margin = 2)

#saveRDS(pos, file = "integratedAnalysis/data_pos_analyzed_ADT.RDS")
pos <- readRDS(file = "integratedAnalysis/data_pos_analyzed_ADT.RDS")

# check out a few markers for visualization
row.names(pos@assays$ADT)

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD3:SK7-CD3E-AHS0033-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD3 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD3E") + ggtitle("CD3 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD19:SJ25C1-CD19-AHS0030-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD19") + ggtitle("CD19 RNA")
p1 | p2

####################################
### CD45 pos: cluster annotation ###
####################################

pos <- readRDS(file = "integratedAnalysis/data_pos_analyzed_ADT.RDS")

### use AbSeq data to help annotate clusters ###
row.names(pos@assays$ADT)

# T cells
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD3:SK7-CD3E-AHS0033-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD3 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD3E") + ggtitle("CD3 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD8:SK1-CD8A-AHS0228-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD8 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD8A") + ggtitle("CD8 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD4:RPA-T4-CD4-AHS0227-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD4 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD4") + ggtitle("CD4 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD279:MIH4-PDCD1-AHS0190-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD279 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "PDCD1") + ggtitle("CD279 RNA")
p1 | p2

# B cells
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD19:SJ25C1-CD19-AHS0030-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD19") + ggtitle("CD19 RNA")
p1 | p2

FeaturePlot(pos, features = "MS4A1")

# myeloid cells
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD15-FUT4-AHS0196-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD15 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "FUT4") + ggtitle("CD15 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD14:M5E2-CD14-AHS0173-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD14 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD14") + ggtitle("CD14 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD16:3G8-FCGR3A-AHS0053-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD16 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "FCGR3A") + ggtitle("CD16 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD11b:M1-70-ITGAM-AHS0005-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD11b protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "ITGAM") + ggtitle("CD11b RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"   # DCs
p1 <- FeaturePlot(pos, "CD11c:B-LY6-ITGAX-AHS0056-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD11c protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "ITGAX") + ggtitle("CD11c RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD206-MRC1-AHS0072-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD206 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "MRC1") + ggtitle("CD206 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD80-CD80-AHS0046-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD80 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD80") + ggtitle("CD80 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD86:2331-CD86-AHS0057-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD86 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD86") + ggtitle("CD86 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "HLA-DR-CD74-AHS0035-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD74 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD74") + ggtitle("CD74 RNA")
p1 | p2

FeaturePlot(pos, features = "LYZ")
FeaturePlot(pos, features = "HLA-DRA")
FeaturePlot(pos, features = "HLA-DPA1")

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD56:B159-NCAM1-AHS0257-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD56 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "NCAM1") + ggtitle("CD56 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "Tim3-HAVCR2-AHS0016-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("Tim3 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "HAVCR2") + ggtitle("Tim3 RNA")
p1 | p2

# CD1a Langerhans cells
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD1a-CD1A-AHS0067-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD1a protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD1A") + ggtitle("CD1a RNA")
p1 | p2

FeaturePlot(pos, features = "CD207")

# DCs
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD11c:B-LY6-ITGAX-AHS0056-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD11c protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "ITGAX") + ggtitle("CD11c RNA")
p1 | p2

# mast cells and basophils
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD117:104D2-KIT-AHS0165-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD117 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "KIT") + ggtitle("CD117 RNA")
p1 | p2

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD69-CD69-AHS0010-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD69 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CD69") + ggtitle("CD69 RNA")
p1 | p2

FeaturePlot(pos, features = "CPA3")
FeaturePlot(pos, features = "FCGR2A")
FeaturePlot(pos, features = "FCER1A")
FeaturePlot(pos, features = "IL13")

# mast cells
FeaturePlot(pos, features = "TPSAB1")
FeaturePlot(pos, features = "CMA1")
FeaturePlot(pos, features = "MITF")

# basophils
FeaturePlot(pos, features = "ENPP3") # CD203c
FeaturePlot(pos, features = "CD63")
FeaturePlot(pos, features = "IL3RA") # CD123
FeaturePlot(pos, features = "CCR3")
FeaturePlot(pos, features = "IL4")
FeaturePlot(pos, features = "CEBPA")

# fibroblasts
FeaturePlot(pos, features = "CD34")
FeaturePlot(pos, features = "COL1A1")
FeaturePlot(pos, features = "LOXL1")
FeaturePlot(pos, features = "LUM")
FeaturePlot(pos, features = "FBLN1")
FeaturePlot(pos, features = "FBLN2")

# immune cells (CD45+)
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD45RA:HI100-PTPRC-AHS0009-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD45RA protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "PTPRC") + ggtitle("CD45 RNA")
p1 | p2
ggsave("plots_analysis/featurePlot_CD45RA_RNA_ADT.png")
ggsave("plots_analysis/featurePlot_CD45RA_RNA_ADT.pdf")

DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD45RO-PTPRC-AHS0036-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD45RO protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "PTPRC") + ggtitle("CD45 RNA")
p1 | p2
ggsave("plots_analysis/featurePlot_CD45RO_RNA_ADT.png")
ggsave("plots_analysis/featurePlot_CD45RO_RNA_ADT.pdf")

# CD193 basophils, eosinophils, Th1 & Th2 T cells
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD193-CCR3-AHS0159-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD193 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "CCR3") + ggtitle("CD193 RNA")
p1 | p2

# CD138 plasma cells
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD138-SDC1-AHS0121-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD138 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "SDC1") + ggtitle("CD138 RNA")
p1 | p2

FeaturePlot(pos, features = "CD27")
FeaturePlot(pos, features = "SLAMF7")

# CD161 Th17 cells (should be mixed with T cells at this clustering resolution)
DefaultAssay(pos) <- "ADT"
p1 <- FeaturePlot(pos, "CD161:HP-3G10-KLRB1-AHS0205-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD161 protein")
DefaultAssay(pos) <- "RNA"
p2 <- FeaturePlot(pos, "KLRB1") + ggtitle("CD161 RNA")
p1 | p2

FeaturePlot(pos, features = "RORC")

# endothelial cells
FeaturePlot(pos, features = "CD34")
FeaturePlot(pos, features = "PECAM1")
FeaturePlot(pos, features = "CDH5")
FeaturePlot(pos, features = "KDR")

# vascular smooth muscle cells
FeaturePlot(pos, features = "ACTA2")
FeaturePlot(pos, features = "CNN1")
FeaturePlot(pos, features = "TAGLN")
FeaturePlot(pos, features = "MYH11")

### find DE markers between mast cell clusters ###
mastDEmarkers <- FindMarkers(pos, ident.1 = "4", ident.2 = "10")
# Positive LFC = feature more highly expressed in the FIRST group.
head(mastDEmarkers, n = 10)

# get significant only
sigMastDEmarkers <- mastDEmarkers[mastDEmarkers$p_val_adj <= 0.05,]

# ascending order
sigMastDEmarkers %>%
  arrange(avg_log2FC) %>%
  slice(1:10)

# descending order
sigMastDEmarkers %>%
  arrange(-avg_log2FC) %>%
  slice(1:10)

### run find markers ###
DefaultAssay(pos) <- "RNA"
markers <- FindAllMarkers(pos, assay = "RNA")
#saveRDS(markers, file = "integratedAnalysis/data_pos_markers.RDS")
markers <- readRDS("integratedAnalysis/data_pos_markers.RDS")

# view top genes by cluster number
head(markers[markers$cluster == "10",], n = 10)

unique(pos$seurat_clusters)
pos$seurat_clusters <- Idents(pos)
unique(pos$seurat_clusters)

pos <- RenameIdents(pos, 
                    `0` = "T cells",
                    `1` = "Myeloid cells",
                    `2` = "Myeloid cells",
                    `3` = "Fibroblasts",
                    `4` = "Mast cells",
                    `5` = "Myeloid cells",
                    `6` = "Endothelial cells",
                    `7` = "Plasma cells",
                    `8` = "B cells",
                    `9` = "Vascular smooth muscle cells",
                    `10` = "Antigen-presenting mast cells")

pos$annot_clusters <- Idents(pos)
unique(pos$annot_clusters)

#saveRDS(pos, file = "integratedAnalysis/data_pos_annotated_sepMastClusters.RDS")

### visualize ###
DefaultAssay(pos) <- "RNA"
cluster_markers <- c("PTPRC","CD3E","CD8A",
                     "CD14","ITGAX","MRC1",
                     "CD86","CD74","LYZ","HLA-DRA",
                     "COL1A1","FBLN1","LUM",
                     "KIT","CD69","CMA1","TPSAB1",
                     "PECAM1","CDH5","KDR","CD34",
                     "SDC1","CD27","SLAMF7",
                     "CD19","MS4A1",
                     "ACTA2","CNN1","TAGLN","MYH11")

legend <- get_legend(DotPlot(pos, features = cluster_markers) +
                       theme(legend.title=element_text(size=10)))
p <- DotPlot(pos, features = cluster_markers) +
  RotatedAxis() +
  ggtitle("CD45+ Marker Gene Expression") +
  theme(legend.position = "",
        axis.title.y = element_blank())
plot_grid(p, legend, nrow = 1, ncol = 2, rel_widths = c(6,1))
ggsave("plots_analysis/dotplot_pos_RNA_sepMastClusters.png", bg = "white", units = "in", width = 10, height = 5)
ggsave("plots_analysis/dotplot_pos_RNA_sepMastClusters.pdf", units = "in", width = 10, height = 5)

DefaultAssay(pos) <- "ADT"
row.names(pos@assays$ADT)
cluster_markers_ADT <- c("CD45RA:HI100-PTPRC-AHS0009-pAbO",
                         "CD45RO-PTPRC-AHS0036-pAbO",
                         "CD3:SK7-CD3E-AHS0033-pAbO",
                         "CD4:RPA-T4-CD4-AHS0227-pAbO",
                         "CD8:SK1-CD8A-AHS0228-pAbO",
                         "CD14:M5E2-CD14-AHS0173-pAbO",
                         "CD16:3G8-FCGR3A-AHS0053-pAbO",
                         "CD11b:M1-70-ITGAM-AHS0005-pAbO",
                         "CD11c:B-LY6-ITGAX-AHS0056-pAbO",
                         "CD1c-CD1C-AHS0088-pAbO",
                         "CD80-CD80-AHS0046-pAbO",
                         "CD86:2331-CD86-AHS0057-pAbO",
                         "HLA-DR-CD74-AHS0035-pAbO",
                         "CD117:104D2-KIT-AHS0165-pAbO",
                         "CD69-CD69-AHS0010-pAbO",
                         "CD138-SDC1-AHS0121-pAbO",
                         "CD19:SJ25C1-CD19-AHS0030-pAbO")
labels <- gsub(":.*","",cluster_markers_ADT)
labels <- gsub("HLA-DR-CD74-AHS0035-pAbO","CD74",labels)
labels <- gsub("CD80-CD80-AHS0046-pAbO","CD80",labels)
labels <- gsub("CD69-CD69-AHS0010-pAbO","CD69",labels)
labels <- gsub("CD138-SDC1-AHS0121-pAbO","CD138",labels)
labels <- gsub("CD45RO-PTPRC-AHS0036-pAbO","CD45RO",labels)
labels <- gsub("CD1c-CD1C-AHS0088-pAbO", "CD1c",labels)
labels
DotPlot(pos, features = cluster_markers_ADT, cols = c("lightgrey", "darkgreen")) +
  scale_x_discrete(labels = labels) +
  RotatedAxis() + ylab("") + xlab("") +
  ggtitle("CD45+ Marker Protein Expression")
ggsave("plots_analysis/dotplot_pos_ADT_sepMastClusters.png", bg = "white")
ggsave("plots_analysis/dotplot_pos_ADT_sepMastClusters.pdf", units = "in",
       width = 8, height = 5)

DefaultAssay(pos) <- "RNA"
dittoDimPlot(pos, var = "annot_clusters") +
  ggtitle("scRNA-seq of CD45+ cells") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size = 11),
        legend.position = "bottom") +
  guides(color=guide_legend(ncol=2, override.aes = list(size=5)))
ggsave("plots_analysis/dimPlot_pos_annot_sepMastClusters.png", bg = "white", 
       units = "in", width = 5, height = 6)
ggsave("plots_analysis/dimPlot_pos_annot_sepMastClusters.pdf")

dittoDimPlot(pos, var = "annot_clusters", split.by = "digest_method") +
  ggtitle("CD45+ cells split by digestion method") +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(rect = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.text=element_text(size=10)) +
  guides(color=guide_legend(ncol=3, override.aes = list(size=5)))
ggsave("plots_analysis/dimPlot_pos_annot_splitDigest_sepMastClusters.png", units = "in", width = 5, height = 5)
ggsave("plots_analysis/dimPlot_pos_annot_splitDigest_sepMastClusters.pdf")

#####################################
### CD45 pos: prop test all cells ###
#####################################

pos <- readRDS(file = "integratedAnalysis/data_pos_annotated_sepMastClusters.RDS")

# first visualize percentage of cells in each cluster
unique(pos$annot_clusters)
dittoBarPlot(pos, var = "annot_clusters", group.by = "digest_method",
             var.labels.reorder = c(8,6,4,5,3,7,2,9,1), x.labels.rotate = FALSE) +
  ggtitle("CD45+ cell type proportions") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.text=element_text(size = 10)) +
  xlab("") +
  RotatedAxis()
ggsave("plots_analysis/barplot_proportions_pos_mastSep.png")
ggsave("plots_analysis/barplot_proportions_pos_mastSep.pdf")

# create the analysis object
prop_test <- sc_utils(pos)

# run testing and bootstrapping-- sample1 ctrl, sample2 treatment
prop_test <- permutation_test(
  prop_test, cluster_identity = "annot_clusters",
  sample_1 = "Simultaneous", sample_2 = "Sequential",
  sample_identity = "digest_method"
)
saveRDS(prop_test, file = "data_pos_prop_test_results_sepMastClusters.RDS")

write.table(prop_test@results, quote = FALSE, sep = "\t", file = "integratedAnalysis/prop_test_results_sepMastClusters.txt")

head(prop_test@results)

### visualize results ###
# significantly up is higher in sequential
# significantly down is lower in sequential
permutation_plot(prop_test) +
  ggtitle("CD45+ Sequential vs. Simultaneous\n") +
  ylab("log2FD") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))
ggsave("plots_analysis/prop_plot_pos_sepMastClusters.png", bg = "white")
ggsave("plots_analysis/prop_plot_pos_sepMastClusters.pdf")

######################################
### analysis CD45 negative dataset ###
######################################

filtered_seurat <- readRDS(file = "integratedAnalysis/data_listOfSeuratObjects_filtered.RDS")

# get just CD45 negative dataset
neg <- filtered_seurat[["DE2_CD45neg"]]

#saveRDS(neg, file = "integratedAnalysis/data_neg.RDS")
neg <- readRDS(file = "integratedAnalysis/data_neg.RDS")

### standard analysis ###
neg <- NormalizeData(neg)
neg <- FindVariableFeatures(neg, selection.method = "vst", nfeatures = 2000)
neg <- ScaleData(neg)
neg <- RunPCA(neg)
ElbowPlot(neg, ndims = 20)
neg <- RunUMAP(neg, dims = 1:8, reduction = "pca")
neg <- FindNeighbors(neg, dims = 1:8)

res_range <- c(seq(0.1,1.0,0.1), seq(1.2,2.0,0.2))
res_range
neg <- FindClusters(neg, resolution = res_range)

head(neg[[]])
clustree(neg, prefix = "RNA_snn_res.", show_axis = TRUE) +
  theme(legend.position="none")
ggsave("plots_analysis/clustree_prebatch_fresh_neg.png", bg = "white", width = 7, height = 12, units = "in")
ggsave("plots_analysis/clustree_prebatch_fresh_neg.pdf", width = 7, height = 12, units = "in")

# set optimal resolution
Idents(neg) <- neg$RNA_snn_res.0.1

#saveRDS(neg, file = "integratedAnalysis/data_neg_analyzed.RDS")
neg <- readRDS(file = "integratedAnalysis/data_neg_analyzed.RDS")

dittoDimPlot(neg, var = Idents(neg), reduction.use = "umap") +
  ggtitle("CD45- Pre-batch correction")
ggsave("plots_analysis/dimPlot_seurat_clusters_fresh_neg.png", bg = "white")
ggsave("plots_analysis/dimPlot_seurat_clusters_fresh_neg.pdf")

# visualize metadata to see if batch correction is needed
dittoDimPlot(neg, var = "digest_method") +
  ggtitle("CD45- Pre-batch correction: digestion method")
ggsave("plots_analysis/dimPlot_digest_method_fresh_neg.png", bg = "white")
ggsave("plots_analysis/dimPlot_digest_method_fresh_neg.pdf")

dittoDimPlot(neg, var = "digest_method", split.by = "digest_method") +
  theme(legend.position = "none") +
  ggtitle("CD45- Pre-batch correction: digestion method")
ggsave("plots_analysis/dimPlot_splitDigestMethod_fresh_neg.png", bg = "white")
ggsave("plots_analysis/dimPlot_splitDigestMethod_fresh_neg.pdf")

dittoDimPlot(neg, var = "patientID") +
  ggtitle("CD45- Pre-batch correction: patient ID")
ggsave("plots_analysis/dimPlot_patientID_fresh_neg.png", bg = "white")
ggsave("plots_analysis/dimPlot_patientID_fresh_neg.pdf")

dittoDimPlot(neg, var = "patientID", split.by = "patientID") +
  ggtitle("CD45- Pre-batch correction: patient ID") +
  theme(legend.position = "none")
ggsave("plots_analysis/dimPlot_splitPatientID_fresh_neg.png", bg = "white")
ggsave("plots_analysis/dimPlot_splitPatientID_fresh_neg.pdf")

### batch correct by patient ID ###
neg <- RunHarmony(neg, 
                  group.by.vars = "patientID",
                  assay.use = "RNA",
                  reduction = "pca",
                  reduction.save = "harmony",
                  plot_convergence = TRUE)
ElbowPlot(neg, reduction = "harmony")
neg <- RunUMAP(neg, dims = 1:8, reduction = "harmony")
neg <- FindNeighbors(neg, dims = 1:8, reduction = "harmony",
                     graph.name = c("RNA.batch_nn","RNA.batch_snn"))

res_range <- c(seq(0.1,1.0,0.1), seq(1.2,2.0,0.2))
res_range
neg <- FindClusters(neg, resolution = res_range, graph.name = "RNA.batch_snn")

head(neg[[]])
clustree(neg, prefix = "RNA.batch_snn_res.", show_axis = TRUE) +
  theme(legend.position="none")
ggsave("plots_analysis/clustree_postbatch_fresh_neg.png", bg = "white", width = 7, height = 12, units = "in")
ggsave("plots_analysis/clustree_postbatch_fresh_neg.pdf", width = 7, height = 12, units = "in")

# set optimal resolution
Idents(neg) <- neg$RNA_snn_res.0.1

#saveRDS(neg, file = "integratedAnalysis/data_neg_analyzed_postbatch.RDS")
neg <- readRDS(file = "integratedAnalysis/data_neg_analyzed_postbatch.RDS")

### visualize post batch correction ###
dittoDimPlot(neg, var = Idents(neg), reduction.use = "umap") +
  ggtitle("CD45- Post-batch correction")
ggsave("plots_analysis/dimPlot_seurat_clusters_fresh_neg_postbatch.png", bg = "white")
ggsave("plots_analysis/dimPlot_seurat_clusters_fresh_neg_postbatch.pdf")

# visualize metadata to see if batch correction is needed
dittoDimPlot(neg, var = "digest_method") +
  ggtitle("CD45- Post-batch correction: digestion method")
ggsave("plots_analysis/dimPlot_digest_method_fresh_neg_postbatch.png", bg = "white")
ggsave("plots_analysis/dimPlot_digest_method_fresh_neg_postbatch.pdf")

dittoDimPlot(neg, var = "digest_method", split.by = "digest_method") +
  theme(legend.position = "none") +
  ggtitle("CD45- Post-batch correction: digestion method")
ggsave("plots_analysis/dimPlot_splitDigestMethod_fresh_neg_postbatch.png", bg = "white")
ggsave("plots_analysis/dimPlot_splitDigestMethod_fresh_neg_postbatch.pdf")

dittoDimPlot(neg, var = "patientID") +
  ggtitle("CD45- Post-batch correction: patient ID")
ggsave("plots_analysis/dimPlot_patientID_fresh_neg_postbatch.png", bg = "white")
ggsave("plots_analysis/dimPlot_patientID_fresh_neg_postbatch.pdf")

dittoDimPlot(neg, var = "patientID", split.by = "patientID") +
  ggtitle("CD45- Post-batch correction: patient ID") +
  theme(legend.position = "none")
ggsave("plots_analysis/dimPlot_splitPatientID_fresh_neg_postbatch.png", bg = "white")
ggsave("plots_analysis/dimPlot_splitPatientID_fresh_neg_postbatch.pdf")

### normalize ADT assay ###
DefaultAssay(neg) <- "ADT"
neg <- NormalizeData(neg, normalization.method = "CLR", margin = 2)

# check out a few markers for visualization
row.names(neg@assays$ADT)

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD138-SDC1-AHS0121-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD138 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "SDC1") + ggtitle("CD138 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD19:SJ25C1-CD19-AHS0030-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD19") + ggtitle("CD19 RNA")
p1 | p2

####################################
### CD45 neg: cluster annotation ###
####################################

unique(neg$seurat_clusters)
neg$seurat_clusters <- Idents(neg)
unique(neg$seurat_clusters)

### run find markers ###
DefaultAssay(neg) <- "RNA"
neg_markers <- FindAllMarkers(neg, assay = "RNA")
#saveRDS(neg_markers, file = "integratedAnalysis/data_neg_markers.RDS")
neg_markers <- readRDS("integratedAnalysis/data_neg_markers.RDS")

# view top genes by cluster number
head(neg_markers[neg_markers$cluster == "3",], n = 10)

# fibroblasts
FeaturePlot(neg, features = "COL1A1")
FeaturePlot(neg, features = "COL1A2")
FeaturePlot(neg, features = "LUM")
FeaturePlot(neg, features = "FBLN1")

# myeloid cells 
FeaturePlot(neg, features = "LYZ")
FeaturePlot(neg, features = "HLA-DRA")
FeaturePlot(neg, features = "HLA-DPA1")

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, features = "HLA-DR-CD74-AHS0035-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD74 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, features = "CD74") + ggtitle("CD74 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD15-FUT4-AHS0196-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD15 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "FUT4") + ggtitle("CD15 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD14:M5E2-CD14-AHS0173-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD14 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD14") + ggtitle("CD14 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD16:3G8-FCGR3A-AHS0053-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD16 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "FCGR3A") + ggtitle("CD16 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD11b:M1-70-ITGAM-AHS0005-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD11b protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "ITGAM") + ggtitle("CD11b RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD11c:B-LY6-ITGAX-AHS0056-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD11c protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "ITGAX") + ggtitle("CD11c RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD206-MRC1-AHS0072-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD206 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "MRC1") + ggtitle("CD206 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD80-CD80-AHS0046-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD80 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD80") + ggtitle("CD80 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD86:2331-CD86-AHS0057-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD86 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD86") + ggtitle("CD86 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD56:B159-NCAM1-AHS0257-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD56 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "NCAM1") + ggtitle("CD56 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "Tim3-HAVCR2-AHS0016-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("Tim3 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "HAVCR2") + ggtitle("Tim3 RNA")
p1 | p2

# CD1a Langerhans cells
DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD1a-CD1A-AHS0067-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD1a protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD1A") + ggtitle("CD1a RNA")
p1 | p2

FeaturePlot(neg, features = "CD207")

# endothelial cells
FeaturePlot(neg, features = "PECAM1")
FeaturePlot(neg, features = "CDH5")
FeaturePlot(neg, features = "KDR")
FeaturePlot(neg, features = "CD34")

# lymphoid cells
DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD3:SK7-CD3E-AHS0033-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD3 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD3E") + ggtitle("CD3 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD8:SK1-CD8A-AHS0228-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD8 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD8A") + ggtitle("CD8 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD4:RPA-T4-CD4-AHS0227-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD4 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD4") + ggtitle("CD4 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD279:MIH4-PDCD1-AHS0190-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD279 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "PDCD1") + ggtitle("CD279 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD19:SJ25C1-CD19-AHS0030-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD19") + ggtitle("CD19 RNA")
p1 | p2

FeaturePlot(neg, features = "MS4A1")

# mast cells and basophils
DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD117:104D2-KIT-AHS0165-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD117 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "KIT") + ggtitle("CD117 RNA")
p1 | p2

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD69-CD69-AHS0010-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD69 protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "CD69") + ggtitle("CD69 RNA")
p1 | p2

FeaturePlot(neg, features = "CPA3")
FeaturePlot(neg, features = "FCGR2A")
FeaturePlot(neg, features = "FCER1A")
FeaturePlot(neg, features = "IL13")

# mast cells
FeaturePlot(neg, features = "TPSAB1")
FeaturePlot(neg, features = "CMA1")
FeaturePlot(neg, features = "MITF")

# basophils
FeaturePlot(neg, features = "ENPP3") # CD203c
FeaturePlot(neg, features = "CD63")
FeaturePlot(neg, features = "IL3RA") # CD123
FeaturePlot(neg, features = "CCR3")
FeaturePlot(neg, features = "IL4")
FeaturePlot(neg, features = "CEBPA")

# vascular smooth muscle cells
FeaturePlot(neg, features = "ACTA2")
FeaturePlot(neg, features = "TAGLN")
FeaturePlot(neg, features = "MYH11")
FeaturePlot(neg, features = "CNN1")

# keratinocytes
FeaturePlot(neg, features = "KRT15")
FeaturePlot(neg, features = "IVL")
FeaturePlot(neg, features = "KRT14")
FeaturePlot(neg, features = "KRT16")
FeaturePlot(neg, features = "DMKN")
FeaturePlot(neg, features = "S100A14")
FeaturePlot(neg, features = "CASP14")

# melanocytes
FeaturePlot(neg, features = "TYRP1")
FeaturePlot(neg, features = "MLANA")
FeaturePlot(neg, features = "PMEL")
FeaturePlot(neg, features = "TYR")

# view CD45 protein
row.names(neg@assays$ADT)
DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD45RA:HI100-PTPRC-AHS0009-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD45A protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "PTPRC") + ggtitle("CD45 RNA")
p1 | p2
ggsave("featurePlot_neg_CD45RA_RNA_ADT.png")
ggsave("featurePlot_neg_CD45RA_RNA_ADT.pdf")

DefaultAssay(neg) <- "ADT"
p1 <- FeaturePlot(neg, "CD45RO-PTPRC-AHS0036-pAbO", cols = c("lightgrey", "darkgreen")) + ggtitle("CD45RO protein")
DefaultAssay(neg) <- "RNA"
p2 <- FeaturePlot(neg, "PTPRC") + ggtitle("CD45 RNA")
p1 | p2
ggsave("featurePlot_neg_CD45RO_RNA_ADT.png")
ggsave("featurePlot_neg_CD45RO_RNA_ADT.pdf")

# rename identities
neg <- RenameIdents(neg, 
                    `0` = "Fibroblasts",
                    `1` = "Endothelial cells",
                    `2` = "Myeloid cells",
                    `3` = "Lymphoid cells",
                    `4` = "Mast cells",
                    `5` = "Vascular smooth muscle cells",
                    `6` = "Keratinocytes",
                    `7` = "Melanocytes")

neg$annot_clusters <- Idents(neg)
unique(neg$annot_clusters)

#saveRDS(neg, file = "integratedAnalysis/data_neg_annotated.RDS")

### visualize ###
DefaultAssay(neg) <- "RNA"
cluster_markers <- c("PTPRC","COL1A1","FBLN1","LUM",
                     "PECAM1","CDH5","KDR","CD34",
                     "CD14","FCGR3A","ITGAM","ITGAX","MRC1",
                     "CD86","CD74","LYZ",
                     "CD3E","CD8A","MS4A1",
                     "KIT","CD69","CMA1","TPSAB1",
                     "ACTA2","CNN1","TAGLN","MYH11",
                     "KRT14","KRT15","DMKN","S100A14","CASP14",
                     "TYRP1","MLANA","PMEL","TYR")

legend <- get_legend(DotPlot(neg, features = cluster_markers) +
                       theme(legend.title=element_text(size=10)))
p <- DotPlot(neg, features = cluster_markers) +
  RotatedAxis() +
  ggtitle("CD45- Marker Gene Expression") +
  theme(legend.position = "",
        axis.title.y = element_blank())
plot_grid(p, legend, nrow = 1, ncol = 2, rel_widths = c(6,1))
ggsave("plots_analysis/dotplot_neg_RNA.png", bg = "white", units = "in", width = 10, height = 5)
ggsave("plots_analysis/dotplot_neg_RNA.pdf", units = "in", width = 10, height = 5)

DefaultAssay(neg) <- "ADT"
row.names(neg@assays$ADT)
cluster_markers_ADT <- c("CD3:SK7-CD3E-AHS0033-pAbO",
                         "CD4:RPA-T4-CD4-AHS0227-pAbO",
                         "CD8:SK1-CD8A-AHS0228-pAbO",
                         "CD19:SJ25C1-CD19-AHS0030-pAbO",
                         "CD14:M5E2-CD14-AHS0173-pAbO",
                         "CD16:3G8-FCGR3A-AHS0053-pAbO",
                         "CD11b:M1-70-ITGAM-AHS0005-pAbO",
                         "CD11c:B-LY6-ITGAX-AHS0056-pAbO",
                         "CD80-CD80-AHS0046-pAbO",
                         "CD86:2331-CD86-AHS0057-pAbO",
                         "HLA-DR-CD74-AHS0035-pAbO",
                         "CD117:104D2-KIT-AHS0165-pAbO",
                         "CD69-CD69-AHS0010-pAbO",
                         "CD138-SDC1-AHS0121-pAbO")
labels <- gsub(":.*","",cluster_markers_ADT)
labels <- gsub("HLA-DR-CD74-AHS0035-pAbO","CD74",labels)
labels <- gsub("CD80-CD80-AHS0046-pAbO","CD80",labels)
labels <- gsub("CD69-CD69-AHS0010-pAbO","CD69",labels)
labels <- gsub("CD138-SDC1-AHS0121-pAbO","CD138",labels)
labels
DotPlot(neg, features = cluster_markers_ADT, cols = c("lightgrey", "darkgreen")) +
  scale_x_discrete(labels = labels) +
  RotatedAxis() +
  ggtitle("CD45- Marker Protein Expression")
ggsave("plots_analysis/dotplot_neg_ADT.png", bg = "white")
ggsave("plots_analysis/dotplot_neg_ADT.pdf")

DefaultAssay(neg) <- "RNA"
dittoDimPlot(neg, var = "annot_clusters") +
  ggtitle("scRNA-seq of CD45- cells") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size = 12))
ggsave("plots_analysis/dimPlot_neg_annot.png")
ggsave("plots_analysis/dimPlot_neg_annot.pdf")

dittoDimPlot(neg, var = "annot_clusters", split.by = "digest_method") +
  ggtitle("CD45- cells split by digestion method") +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(rect = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.text=element_text(size=10)) +
  guides(color=guide_legend(ncol=2, override.aes = list(size=5)))
ggsave("plots_analysis/dimPlot_neg_annot_splitDigest.png", units = "in", width = 5, height = 5)
ggsave("plots_analysis/dimPlot_neg_annot_splitDigest.pdf")

#####################################
### CD45 neg: prop test all cells ###
#####################################

neg <- readRDS("integratedAnalysis/data_neg_annotated.RDS")

# first visualize percentage of cells in each cluster
unique(neg$annot_clusters)
dittoBarPlot(neg, var = "annot_clusters", group.by = "digest_method",
             var.labels.reorder = c(2,1,7,4,5,8,3,6), x.labels.rotate = FALSE) +
  ggtitle("CD45- cell type proportions") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  xlab("")
ggsave("plots_analysis/barplot_proportions_neg.png")
ggsave("plots_analysis/barplot_proportions_neg.pdf")

# create the analysis object
prop_test <- sc_utils(neg)

# run testing and bootstrapping
# sample1 control, sample2 treatment
prop_test <- permutation_test(
  prop_test, cluster_identity = "annot_clusters",
  sample_1 = "Simultaneous", sample_2 = "Sequential",
  sample_identity = "digest_method"
)
#saveRDS(prop_test, file = "integratedAnalysis/data_neg_prop_test_results.RDS")

write.table(prop_test@results, file = "integratedAnalysis/prop_test_results_neg.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

head(prop_test@results)

### visualize results ###
# significantly up is higher in sequential
# significantly down is lower in sequential (none)
permutation_plot(prop_test) +
  ggtitle("CD45- Sequential vs. Simultaneous") +
  ylab("log2FD") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))
ggsave("plots_analysis/prop_plot_neg.png", bg = "white")
ggsave("plots_analysis/prop_plot_neg.pdf")

#################################
### CD45 pos: T cell analysis ###
#################################

pos <- readRDS("integratedAnalysis/data_pos_annotated_sepMastClusters.RDS")

# subset T cells
unique(Idents(pos))
Tcells <- subset(pos, idents = "T cells")

# drop levels left over from full object
unique(Tcells$annot_clusters)
Tcells$annot_clusters <- droplevels(Tcells$annot_clusters)
unique(Tcells$annot_clusters)

# recluster T cells
DefaultAssay(Tcells) <- "integrated"
Tcells <- ScaleData(Tcells, vars.to.regress = "patientID", features = row.names(Tcells))
Tcells <- RunPCA(Tcells, features = row.names(Tcells))
ElbowPlot(Tcells)

Tcells <- RunUMAP(Tcells, dims = 1:10)
Tcells <- FindNeighbors(Tcells, dims = 1:10)
res_range <- c(seq(0.1,0.5,0.05), seq(0.6, 1.0, 0.1))
res_range
Tcells <- FindClusters(Tcells, res = res_range)

# run clustree to find optimal resolution
head(Tcells[[]]) # use to find prefix for resolutions
clustree(Tcells, prefix = "integrated_snn_res.") +
  guides(edge_colour = "none", edge_alpha = "none", size = "none")
ggsave("plots_analysis/clustree_Tcells.png", bg = "white")
ggsave("plots_analysis/clustree_Tcells.pdf")

saveRDS(Tcells, file = "integratedAnalysis/Tcells.RDS")

# set optimal ident
Idents(Tcells) <- Tcells$integrated_snn_res.0.15

dittoDimPlot(Tcells, var = Idents(Tcells), reduction = "umap")

DefaultAssay(Tcells) <- "RNA"
FeaturePlot(Tcells, features = "CD4")
FeaturePlot(Tcells, features = "CD8A")
FeaturePlot(Tcells, features = "RORC")
FeaturePlot(Tcells, features = "PTPRC")
FeaturePlot(Tcells, features = "ITGAE")
FeaturePlot(Tcells, features = "CD3E")
FeaturePlot(Tcells, features = "FOXP3")
FeaturePlot(Tcells, features = "IL17A")
FeaturePlot(Tcells, features = "IL2RG")

FeaturePlot(Tcells, features = "CCR8")
FeaturePlot(Tcells, features = "CCR7")
FeaturePlot(Tcells, features = "SELPLG")

gene_markers <- c("PTPRC","CD4", "CD8A","ITGAE","CD3E")
DotPlot(Tcells, features = gene_markers)

DefaultAssay(Tcells) <- "ADT"
row.names(Tcells@assays$ADT)
FeaturePlot(Tcells, features = "CD3:SK7-CD3E-AHS0033-pAbO",
            cols = c("lightgrey","darkgreen"))
FeaturePlot(Tcells, features = "CD4:RPA-T4-CD4-AHS0227-pAbO",
            cols = c("lightgrey","darkgreen"))
FeaturePlot(Tcells, features = "CD8:SK1-CD8A-AHS0228-pAbO",
            cols = c("lightgrey","darkgreen"))

FeaturePlot(Tcells, features = "TCR-gamma-delta:B1-TRD-TRG-AHS0015-pAbO",
            cols = c("lightgrey","darkgreen"))

FeaturePlot(Tcells, features = "CD161:HP-3G10-KLRB1-AHS0205-pAbO",
            cols = c("lightgrey","darkgreen"))

# memory T cells
FeaturePlot(Tcells, features = , "CD45RO-PTPRC-AHS0036-pAbO",
            cols = c("lightgrey","darkgreen"))

# TRMs
FeaturePlot(Tcells, features = , "CD69-CD69-AHS0010-pAbO",
            cols = c("lightgrey","darkgreen"))

# naive T cells
FeaturePlot(Tcells, features = , "CD45RA:HI100-PTPRC-AHS0009-pAbO",
            cols = c("lightgrey","darkgreen"))
FeaturePlot(Tcells, features = , "CD62L:SK11-SELL-AHS0265-pAbO",
            cols = c("lightgrey","darkgreen"))



DefaultAssay(Tcells) <- "RNA"
Tcell_markers <- FindAllMarkers(Tcells)
head(Tcell_markers[Tcell_markers$cluster == "1",], n = 10)

Tcell_markers[Tcell_markers$cluster == "0",] %>%
  arrange(-avg_log2FC) %>%
  slice(1:10)

Tcells <- RenameIdents(Tcells,
                       `0` = "ANXA1+ CD4 T cells",
                       `1` = "KYNU+ CD4 T cells",
                       `2` = "Regulatory T cells",
                       `3` = "CD8 T cells")

unique(Tcells$annot_clusters)  # all T cells, from before
Tcells$annot_clusters <- Idents(Tcells)
unique(Tcells$annot_clusters)

#saveRDS(Tcells, file = "integratedAnalysis/Tcells_annot.RDS")

dittoDimPlot(Tcells, var = "annot_clusters") +
  ggtitle("T cell clusters") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  guides(color = guide_legend(ncol=2, override.aes = list(size=5)))
ggsave("plots_analysis/dimplot_Tcells.png", units = "in", 
       width = 5, height = 6)

gene_markers <- c("PTPRC","CD3E","CD4","ANXA1","KYNU","FOXP3","CD8A")
DotPlot(Tcells, features = gene_markers) +
  RotatedAxis() +
  ggtitle("T cell Gene Expression Markers") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("") +
  xlab("Genes")
ggsave("plots_analysis/dotplot_Tcell_RNA.png")
ggsave("plots_analysis/dotplot_Tcell_RNA.pdf")

## prop test ##
Tcells <- readRDS(file = "integratedAnalysis/Tcells_annot.RDS")

# first visualize percentage of cells in each cluster
unique(Tcells$annot_clusters)
dittoBarPlot(Tcells, var = "annot_clusters", group.by = "digest_method",
             var.labels.reorder = c(1,3,4,2), x.labels.rotate = FALSE) +
  ggtitle("T cell type proportions") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.text=element_text(size = 10),
        legend.position = "bottom") +
  xlab("") +
  guides(fill = guide_legend(ncol=2, override.aes = list(size=5)))
ggsave("plots_analysis/barplot_proportions_Tcells.png", units = "in", height = 5, width = 4)
ggsave("plots_analysis/barplot_proportions_Tcells.pdf", units = "in", height = 5, width = 4)

# create the analysis object
prop_test <- sc_utils(Tcells)

# run testing and bootstrapping-- sample1 ctrl, sample2 treatment
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "annot_clusters",
  sample_1 = "Simultaneous", sample_2 = "Sequential",
  sample_identity = "digest_method"
)
write.table(prop_test@results, file = "integratedAnalysis/prop_test_results_Tcells.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
#saveRDS(prop_test, file = "data_Tcells_prop_test_results.RDS")

head(prop_test@results)

### visualize results ###
# significantly up is higher in sequential
# significantly down is lower in sequential
permutation_plot(prop_test) +
  ggtitle("T cells: Sequential vs. Simultaneous") +
  ylab("log2FD") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))
ggsave("plots_analysis/prop_plot_Tcells.png", bg = "white")
ggsave("plots_analysis/prop_plot_Tcells.pdf")

###########################
### QC by digest method ###
###########################

### CD45 positive ###
pos <- readRDS("integratedAnalysis/data_pos_annotated_sepMastClusters.RDS")

## overview plot
unique(Idents(pos))
Idents(pos) <- "digest_method"
feats <- c("nCount_RNA","nFeature_RNA","percent_mito")
titles <- c("Number of\nMolecules","Number of\nUnique Genes","Percent\nMitochondrial Genes")
plist <- list()
for (i in seq(1:length(feats))){
  p <- VlnPlot(pos, features = feats[i], pt.size = 0, cols = c("coral","turquoise")) +
    ggtitle(titles[i]) + xlab("") + NoLegend()
  plist[[i]] <- p
}
plot_grid(plotlist= plist, nrow = 1)
ggsave("plots_analysis/violin_QC_digestMethod_overview_pos.pdf")
Idents(pos) <- "annot_clusters"

## split by digest method
# number of UMIs
VlnPlot(pos, features = "nCount_RNA", split.by = "digest_method",
          split.plot = TRUE, pt.size = 0, cols = c("coral","turquoise")) +
  xlab("") +
  ggtitle("Number of Molecules")
ggsave("plots_analysis/violin_QC_digestMethod_nCount_pos.pdf")

# number of unique genes
VlnPlot(pos, features = "nFeature_RNA", split.by = "digest_method", 
        split.plot = TRUE, pt.size = 0, cols = c("coral","turquoise")) +
  xlab("") +
  ggtitle("Number of Unique Genes")
ggsave("plots_analysis/violin_QC_digestMethod_nFeat_pos.pdf")

# percent mitochondrial genes
VlnPlot(pos, features = "percent_mito", split.by = "digest_method", 
        split.plot = TRUE, pt.size = 0, cols = c("coral","turquoise")) +
  xlab("") +
  ggtitle("Percent Mitochondrial Genes")
ggsave("plots_analysis/violin_QC_digestMethod_mito_pos.pdf")

### CD45 negative ###
neg <- readRDS("integratedAnalysis/data_neg_annotated.RDS")

## overview plot
unique(Idents(neg))
Idents(neg) <- "digest_method"
feats <- c("nCount_RNA","nFeature_RNA","percent_mito")
titles <- c("Number of\nMolecules","Number of\nUnique Genes","Percent\nMitochondrial Genes")
plist <- list()
for (i in seq(1:length(feats))){
  p <- VlnPlot(neg, features = feats[i], pt.size = 0, cols = c("coral","turquoise")) +
    ggtitle(titles[i]) + xlab("") + NoLegend()
  plist[[i]] <- p
}
plot_grid(plotlist= plist, nrow = 1)
ggsave("plots_analysis/violin_QC_digestMethod_overview_neg.pdf")
Idents(neg) <- "annot_clusters"

### split by digest method ###
# number of UMIs
VlnPlot(neg, features = "nCount_RNA", split.by = "digest_method",
        split.plot = TRUE, pt.size = 0, cols = c("coral","turquoise")) +
  xlab("") +
  ggtitle("Number of Molecules")
ggsave("plots_analysis/violin_QC_digestMethod_nCount_neg.pdf")

# number of unique genes
VlnPlot(neg, features = "nFeature_RNA", split.by = "digest_method", 
        split.plot = TRUE, pt.size = 0, cols = c("coral","turquoise")) +
  xlab("") +
  ggtitle("Number of Unique Genes")
ggsave("plots_analysis/violin_QC_digestMethod_nFeat_neg.pdf")

# percent mitochondrial genes
VlnPlot(neg, features = "percent_mito", split.by = "digest_method", 
        split.plot = TRUE, pt.size = 0, cols = c("coral","turquoise")) +
  xlab("") +
  ggtitle("Percent Mitochondrial Genes")
ggsave("plots_analysis/violin_QC_digestMethod_mito_neg.pdf")


table(pos$digest_method)

