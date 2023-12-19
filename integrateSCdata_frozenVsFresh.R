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
    png(file = paste("plots_analysis/violinPlot_",outbase,"_",name,".png", sep = ""), bg = "white")
    p <- plot_grid(plotlist = plist, ncol = 2) + plot_annotation(title = name)
    print(p)
    dev.off()

    pdf(file = paste("plots_analysis/violinPlot_",outbase,"_",name,".pdf", sep = ""))
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

# read input files
fileList_RSEC <- "fileList.RSEC.txt"  # text file with name of each RSEC file on a separate line
fileList_SMK <- "fileList.SampleTag.txt"  # text file with name of each sampleTag file on a separate line
metadataFile <- "metadata.txt"
patientMetadataFile <- "metadata_patientLevel.txt"

RSEC_container <- readRSECfiles(fileList_RSEC) # this takes a long time
names(RSEC_container)
names(RSEC_container) <- gsub("-", "_", names(RSEC_container))
names(RSEC_container)

#saveRDS(RSEC_container, file = "integratedAnalysis/data_RSECdataframes_frozenFresh.RDS")
RSEC_container <- readRDS(file = "integratedAnalysis/data_RSECdataframes_frozenFresh.RDS")

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
patient_metadata <- patient_metadata %>% 
  mutate(Sample_Tag = case_when(
    as.numeric(SMK) < 10 ~ paste("SampleTag0", SMK, "_hs", sep = ""),
    as.numeric(SMK) >= 10 ~ paste("SampleTag", SMK, "_hs", sep = ""),
    SMK %in% c("Multiplet","Undetermined") ~ SMK))

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
  sample_tag_df <- sampleTag_container[[name]]
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

#saveRDS(seurat_container, file = "integratedAnalysis/data_listOfSeuratObjects_frozenFresh.RDS")
seurat_container <- readRDS("integratedAnalysis/data_listOfSeuratObjects_frozenFresh.RDS")

#######################
### quality control ###
#######################

seurat_container <- readRDS("integratedAnalysis/data_listOfSeuratObjects_frozenFresh.RDS")

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
ggsave(filename = "plots_analysis/barPlotCellCountPerDataset.pdf")
ggsave(filename = "plots_analysis/barPlotCellCountPerDataset.png", bg = "white")

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
  ggsave(filename = paste("plots_analysis/barPlotCellCountPerSMK_",name,".pdf",sep=""))
  ggsave(filename = paste("plots_analysis/barPlotCellCountPerSMK_",name,".png",sep=""), bg = "white")
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
ggsave("plots_analysis/violinPlot_prefilter_allDatasets_nFeature.png", height = 5, width = 8, bg = "white")
ggsave("plots_analysis/violinPlot_prefilter_allDatasets_nFeature.pdf")

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

# by digest method
unique(Idents(obj))
Idents(obj) <- "digest_method"
unique(Idents(obj))
feats = c("nCount_RNA","nFeature_RNA","percent_mito")
titles = c("Number of\nMolecules","Number of\nUnique Genes","Percent\nMitochondrial Genes")
plist <- list()
for (i in seq(1,length(feats))){
  p <- VlnPlot(obj, features = feats[[i]],
          pt.size = 0, cols = c("coral","turquoise")) + 
    NoLegend() + ggtitle(titles[[i]]) + xlab("") +
    theme(axis.text.x = element_blank())
  plist[[i]] <- p
}
legend <- get_legend(VlnPlot(obj, features = feats[[i]],
                             cols = c("coral","turquoise")))
blank_plot <- ggplot() + theme_void()
length(plist)
plist[[4]] <- blank_plot
plist[[5]] <- legend
plist[[6]] <- blank_plot
plot_grid(plotlist = plist, nrow = 2, rel_heights = c(6,1))
ggsave("plots_analysis/violin_QC_digestMethod_prefilter.pdf")

# by tissue state
unique(Idents(obj))
Idents(obj) <- "tissue_state"
unique(Idents(obj))
feats = c("nCount_RNA","nFeature_RNA","percent_mito")
titles = c("Number of\nMolecules","Number of\nUnique Genes","Percent\nMitochondrial Genes")
plist <- list()
for (i in seq(1,length(feats))){
  p <- VlnPlot(obj, features = feats[[i]],
               pt.size = 0, cols = c("coral","turquoise")) + 
    NoLegend() + ggtitle(titles[[i]]) + xlab("") +
    theme(axis.text.x = element_blank())
  plist[[i]] <- p
}
legend <- get_legend(VlnPlot(obj, features = feats[[i]],
                             cols = c("coral","turquoise")) +
                       guides(fill=guide_legend(ncol=2)))
blank_plot <- ggplot() + theme_void()
length(plist)
plist[[4]] <- blank_plot
plist[[5]] <- legend
plist[[6]] <- blank_plot
plot_grid(plotlist = plist, nrow = 2, rel_heights = c(6,1))
ggsave("plots_analysis/violin_QC_tissueState_prefilter.pdf", 
       units = "in", height = 5, width = 8)

# number of cells per digest method
obj@meta.data %>% group_by(digest_method) %>% summarise(n())

# number of cells per tissue state
obj@meta.data %>% group_by(tissue_state) %>% summarise(n())

# by tissue state -- CD45+ simultaneous only
simul <- subset(obj, digest_method == "Simultaneous")
simul <- subset(simul, CD45_isolation == "positive")
unique(simul$CD45_isolation)
unique(simul$digest_method)

unique(Idents(simul))
Idents(simul) <- "tissue_state"
unique(Idents(simul))
feats = c("nCount_RNA","nFeature_RNA","percent_mito")
titles = c("Number of\nMolecules","Number of\nUnique Genes","Percent\nMitochondrial Genes")
plist <- list()
for (i in seq(1,length(feats))){
  p <- VlnPlot(simul, features = feats[[i]],
               pt.size = 0, cols = c("coral","turquoise")) + 
    NoLegend() + ggtitle(titles[[i]]) + xlab("") +
    theme(axis.text.x = element_blank())
  plist[[i]] <- p
}
legend <- get_legend(VlnPlot(simul, features = feats[[i]],
                             cols = c("coral","turquoise")) +
                       guides(fill=guide_legend(ncol=2)))
blank_plot <- ggplot() + theme_void()
length(plist)
plist[[4]] <- blank_plot
plist[[5]] <- legend
plist[[6]] <- blank_plot
plot_grid(plotlist = plist, nrow = 2, rel_heights = c(6,1))
ggsave("plots_analysis/violin_QC_tissueState_simultaneous_CD45pos_prefilter.pdf", 
       units = "in", height = 5, width = 8)