# loading required packages and data --------------------------------------

pkgs <- c("tidyverse","pheatmap","ggrepel","RColorBrewer","Seurat","ggsignif","ggVennDiagram","SummarizedExperiment","SingleCellExperiment")
paks<-lapply(pkgs, function(x) suppressMessages(require(x, character.only = TRUE)))
rm(paks)

# This file has the genes differentially expressed in old specific pericentral hepatocyte clusters compared to other cells
pcdge <- read.table(file="Pericentral_oldspecific_DGE.txt",sep="\t", header = T)

apply_seurat_workflow <- function(seurat_obj, label, minGenes = 100, minReads = 1000, nVar = 2000, res = 1, use_sctransform = FALSE){
  #library(future)
  #plan("multiprocess", workers = 16)
  tmp <- seurat_obj
  tmp$method <- label
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  
  expr <- FetchData(object = tmp, vars = c("nFeature_RNA","nCount_RNA"))
  tmp <- tmp[, which(x = expr$nFeature_RNA > minGenes & expr$nCount_RNA > minReads)]
  
  if(use_sctransform){
    tmp <- SCTransform(tmp, variable.features.n = nVar,return.only.var.genes = FALSE )
  }else{
    tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = nVar)
    tmp <- ScaleData(tmp, vars.to.regress = c("nCount_RNA"))
  }
  
  tmp <- RunPCA(object = tmp, npcs = 40)
  
  tmp <- FindNeighbors(object = tmp, reduction = "pca", dims = 1:40)
  
  tmp <- FindClusters(object = tmp, resolution = res,group.singletons = FALSE,algorithm = 1)
  
  
  tmp <- RunTSNE(object = tmp)
  tmp <- RunUMAP(object = tmp, dims = 1:40 )
  return(tmp)
}

tms.facs <- readRDS("seurat_objects/TM_senis_liver_facs_Seurat_raw.rds")
tms.facs.counts <- GetAssayData(tms.facs, slot = "counts", assay="RNA")
tms.meta <- tms.facs@meta.data
tms.meta$cellID <- rownames(tms.meta)
keep.tms <- 
  tms.meta %>% 
  filter(tissue == "Liver" & 
           sex == "male") 
keep.tms %>% 
  group_by(sex,age) %>% dplyr::count()

mat <- as.matrix(tms.facs.counts)[,keep.tms$cellID]
mat <- mat[rowSums(mat)>0,]
dim(mat)

tms<-Seurat::CreateSeuratObject(counts = mat, meta.data = keep.tms,min.cells=3)
tms.sct.male <- apply_seurat_workflow(tms, label = "TMS-male", minGenes = 200, minReads = 1000, nVar = 6000, res = 0.8, use_sctransform = TRUE)

DimPlot(tms.sct.male, group.by = c("cell.ontology.class","age"), reduction = "umap")

FeaturePlot(tms.sct.male, c(pcdge$GeneName[1:5],"Cyp2e1","Glul","Cyp2f2"), order = T, ncol = 4)
ggsave(filename="TMS_male_FACS_PClowDGE_featureplots.pdf", height = 7, width = 14)
