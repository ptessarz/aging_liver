# use R4 and above for this script


# load libraries and other data  ----------------------

library(Seurat)
library(data.table)
library(ggplot2)
library(ggpointdensity)
library(viridisLite)
library(tidyverse)
library(broom)
library(inflection)


apply_seurat_workflow <- function(seurat_obj, label, minGenes = 100, minReads = 1000, nVar = 2000, res = 1, use_sctransform = FALSE){
  #library(future)
  #plan("multiprocess", workers = 16)
  bla <- seurat_obj
  bla$method <- label
  bla[["percent.mt"]] <- PercentageFeatureSet(bla, pattern = "^MT-")
  
  expr <- FetchData(object = bla, vars = c("nFeature_RNA","nCount_RNA"))
  bla <- bla[, which(x = expr$nFeature_RNA > minGenes & expr$nCount_RNA > minReads)]
  
  if(use_sctransform){
    bla <- SCTransform(bla, variable.features.n = nVar,return.only.var.genes = FALSE )
  }else{
    bla <- NormalizeData(bla, normalization.method = "LogNormalize", scale.factor = 10000)
    bla <- FindVariableFeatures(bla, selection.method = "vst", nfeatures = nVar)
    bla <- ScaleData(bla, vars.to.regress = c("nCount_RNA"))
  }
  
  bla <- RunPCA(object = bla, npcs = 40)
  
  bla <- FindNeighbors(object = bla, reduction = "pca", dims = 1:40)
  
  bla <- FindClusters(object = bla, resolution = res,group.singletons = FALSE,algorithm = 1)
  
  
  bla <- RunTSNE(object = bla)
  bla <- RunUMAP(object = bla, dims = 1:40 )
  return(bla)
}

#The following files are from zUMIs output
genenames <- fread("ss3nuclei.gene_names.txt")
rds <- readRDS("ss3nuclei.dgecounts.rds")
stats <- fread("ss3nuclei.readspercell.txt")
anno.ss3nuc <- read.table(file = "ss3_nuclei_annotation_barcodes.txt", header = T, sep = "\t", stringsAsFactors = F)

ss3nucleikept_barcodes <- read.csv("ss3nucleikept_barcodes.txt")


# basic descriptive stats  ------------------------------------------------


stats <- dcast(stats, RG~type, value.var = 'N', fill = 0)
stats[,User := NULL][,nReads := Exon+Intergenic+Intron+Unmapped][,pct_exon := (Exon/nReads)*100]
stats <- stats[RG != "bad"]


stats <- merge(stats, anno.ss3nuc, by.x = 'RG', by.y = 'XC_DNBPE')

umis <- as.matrix(rds$umicount$exon$all)
reads <- as.matrix(rds$readcount$exon$all)
stats <- stats[RG %in% colnames(umis)]

stats[, nGenes := colSums(umis[,RG]>0)][, nUMIs := colSums(umis[,RG])][, sparsityIndex := colSums(umis[,RG]>0)/colSums(umis[,RG])]

cells_good <- stats[nReads >= 1000 & pct_exon >= 50,]



# Analyse using Seurat workflow ----------------------------------------------------------------

reads_good <- reads[,cells_good$RG]
reads_good <- reads_good[which(rowSums(reads_good)>0),]

genenames_match <- genenames[match(row.names(reads_good),gene_id)]
genenames_match[duplicated(gene_name), gene_name:= paste0(gene_name,sample(genenames_match[duplicated(gene_name),.N], replace = F))]
row.names(reads_good) <- genenames_match$gene_name

ss3reads <- CreateSeuratObject(counts = reads_good, project = "ss3nuclei", min.cells = 0)

metadat <- cells_good %>% 
  dplyr::select(nReads,pct_exon,nuclei,Age,Replicate,RG) 
row.names(metadat) <- cells_good$RG

ss3reads <- AddMetaData(ss3reads, metadat)
ss3reads_x <- apply_seurat_workflow(ss3reads, label = "Smart-seq3", minGenes = 200, minReads = 1000, nVar = 6000, res = 0.8, use_sctransform = TRUE)
ss3meta_reads <- ss3reads_x@meta.data


clustermarkers_reads <- FindAllMarkers(ss3reads_x, only.pos = T)

topmarkers_reads <- clustermarkers_reads %>% group_by(cluster) %>% slice_max(n=2, order_by = avg_log2FC)


saveRDS(ss3reads_x, file = "ss3nuclei/polynucleiSS3_allcelltypes_reads_seurat.rds")

# The same seurat workflow was used to save objects without macrophages.

# check DGE between lowcov pericentral and pericentral --------------------

heps <- readRDS("ss3nuclei/polynucleiSS3_noMacrophages_reads_seurat.rds")

PClowmarkers <- FindMarkers(object = heps,ident.1 = c(3,4), ident.2 = c(0,1,2))
PClowmarkers.filt <- PClowmarkers %>% 
  filter( p_val_adj <= 0.05 & 
           abs(avg_log2FC) >= 1.5) %>% 
  mutate(GeneName = rownames(.))


avg.exp.dge <-AverageExpression(
  heps,
  assay="SCT",
  features = PClowmarkers.filt$GeneName,
  return.seurat = FALSE,
  group.by = c("ident"),
  slot = "data",
)
avg.exp.dge.sct <- data.frame(avg.exp.dge, GeneName = rownames(avg.exp.dge$SCT))
gather(avg.exp.dge.sct, key="Attribute", value="norm.exp",-GeneName) %>% 
  separate(.,col="Attribute", into=c("Assay","ident")) %>% 
  ggplot(., aes(ident,norm.exp,fill=ident))+
  geom_boxplot(notch = T,width=0.8)+
  scale_y_log10()+
  geom_signif(comparisons = list(c("3", "0"),c("3","1"),c("3","2"),c("4", "0"),c("4","1"),c("4","2"),c("3","4")),step_increase = 0.1)+
  ylab("Average expression")+
  theme_classic()+
  theme(text=element_text(size=18))
ggsave(filename = "hepatocytes_PClowcov_DEgenes_boxplot_celltype.pdf", height = 5, width = 7)
write.table(PClowmarkers.filt, file = "Pericentral_lowcov_marker_DGE.txt", sep = "\t", quote = F, row.names = F, col.names = T)

## Check general GO pathway enrichment for any genelist ##
library(ReactomePA)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

ensembl2eg<-AnnotationDbi::select(org.Mm.eg.db,keys=PClowmarkers.filt$GeneName,keytype="SYMBOL",columns=c("ENTREZID"))


gobpreduced <- simplify(enrichGO(ensembl2eg$ENTREZID,'org.Mm.eg.db',ont="BP",readable = T))

p <- pairwise_termsim(gobpreduced)
treeplot(p)
emapplot(p)
gobpbar <- gobpreduced@result
gobpbar %>% 
  dplyr::arrange(., desc(Count)) %>% 
  head(n=20) %>% 
  ggplot(.,aes(x=log10(p.adjust),y=reorder(Description,-p.adjust),fill=Count))+
  geom_col()+
  ylab("")+xlab("Adjusted p-value")+
  ggtitle("change in mean expression")+
  theme_bw()+
  theme(text = element_text(size = 18),
        legend.position = c(0.2, 0.2))


# Featureplots of some differentially dispersed/expressed genes -----------

ss3 <- readRDS("polynucleiSS3_allcelltypes_reads_seurat.rds")
xl <- readxl::read_xlsx(path = "featureplots.xlsx")

FeaturePlot(ss3, features = xl$overdispersed, order = T)
ggsave(filename="overdispersed_featureplots.pdf", width=12, height = 7)
RidgePlot(ss3, features = na.omit(xl$overdispersed), group.by = "Age", stack = F, cols = c("#800040","#1B81E5"))
ggsave(filename="overdispersed_ridgeplots.pdf", width=12, height = 7)

FeaturePlot(ss3, features = na.omit(xl$underdispersed), order = T)
ggsave(filename="underdispersed_featureplots.pdf", width=12, height = 7)
RidgePlot(ss3, features = na.omit(xl$underdispersed), group.by = "Age", stack = F, cols = c("#800040","#1B81E5"))
ggsave(filename="underdispersed_ridgeplots.pdf", width=12, height = 7)

FeaturePlot(ss3, features = na.omit(xl$overexpressed), order = T, ncol = 4)
ggsave(filename="overexpressed_featureplots.pdf", width=16, height = 7)
RidgePlot(ss3, features = na.omit(xl$overexpressed), group.by = "Age", stack = F, ncol=4, cols = c("#800040","#1B81E5"))
ggsave(filename="overexpressed_ridgeplots.pdf", width=16, height = 7)

FeaturePlot(ss3, features = na.omit(xl$underexpressed), order = T, ncol = 4)
ggsave(filename="underexpressed_featureplots.pdf", width=16, height = 9)
RidgePlot(ss3, features = na.omit(xl$underexpressed), group.by = "Age", stack = F, ncol=4, cols = c("#800040","#1B81E5"))
ggsave(filename="underexpressed_ridgeplots.pdf", width=16, height = 9)
