#This code was used to analyse the spatial data mentioned in the manuscript and it's corresponding figures.
#after general processing the code is structured to generate the figures in the manuscript in chronological order. 
#to increse readability to users code has been rephrased to use standard R style not tidyverse style.
library(Seurat)
library(ggrepel)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(svglite)
cbp2 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

## Loading data and separate normalization 
young1 <- Load10X_Spatial(data.dir = "/141895/outs",
                          slice  = "/141895/outs/spatial",
                          assay="RNA")

young2 <- Load10X_Spatial(data.dir = "/141897/outs",
                          slice  = "/141897/outs/spatial",
                          assay="RNA")

old1 <- Load10X_Spatial(data.dir = "141896/outs",
                        slice  = "/141896/outs/spatial",
                        assay="RNA")

old2<- Load10X_Spatial(data.dir = "/141898/outs",
                       slice  = "/141898/outs/spatial",
                       assay="RNA")

plan("multisession", workers = 4)
young1 <- NormalizeData(young1)
young1 <-  FindVariableFeatures(young1, selection.method = "vst", nfeatures = 2000)
plan("multisession", workers = 4)

young2 <- NormalizeData(young2)
young2 <-  FindVariableFeatures(young2, selection.method = "vst", nfeatures = 2000)
plan("multisession", workers = 4)

old1 <- NormalizeData(old1)
old1 <-  FindVariableFeatures(old1, selection.method = "vst", nfeatures = 2000)
plan("multisession", workers = 4)

old2 <- NormalizeData(old2)
old2 <-  FindVariableFeatures(old2, selection.method = "vst", nfeatures = 2000)

#Figure 1B 
a <- SpatialFeaturePlot(young1, features = c("Glul"),pt.size.factor = 1.8,ncol=1,alpha = 0) & scale_fill_gradientn(colours = c("grey","gold1","red"), breaks = c(0.0,1.0,2.0,3.0,4.0,5.0), limits=c(0,5)) & 
  theme( legend.title=element_text(size=13) )

b <- SpatialFeaturePlot(young1, features = c("Glul","Cyp2e1","Cyp2f2"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )
b+a

c <- SpatialFeaturePlot(old2, features = c("Glul"),pt.size.factor = 1.8,ncol=1,alpha = 0) & scale_fill_gradientn(colours = c("grey","gold1","red"), breaks = c(0.0,1.0,2.0,3.0,4.0,5.0), limits=c(0,5)) & 
  theme( legend.title=element_text(size=13) )

d <- SpatialFeaturePlot(old2, features = c("Glul","Cyp2e1","Cyp2f2"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )
d+c

#Figure 1C 

consensus_genes_y <- rownames(young1)[rownames(young1) %in% rownames(young2)]
consensus_genes_o <- rownames(old1)[rownames(old1) %in% rownames(old2)]
consensus_genes <- consensus_genes_y[consensus_genes_y %in% consensus_genes_o]

## integrating young
young.list <- list(young1,young2)

features.young <- SelectIntegrationFeatures(object.list =young.list , nfeatures = length(rownames(young1)))
young.anchors <- FindIntegrationAnchors(object.list = young.list, normalization.method = "LogNormalize", 
                                        anchor.features = features.young,reduction = "cca")
young.combined.sct <- IntegrateData(anchorset = young.anchors, normalization.method = "LogNormalize",features.to.integrate = consensus_genes)
DefaultAssay(young.combined.sct) <- "integrated"

young.combined.sct[["percent.mt"]] <- PercentageFeatureSet(young.combined.sct, pattern = "^mt-",assay = "RNA")


young.combined.sct <- ScaleData(young.combined.sct, verbose = FALSE,vars.to.regress = "percent.mt")
young.combined.sct <- RunPCA(young.combined.sct, verbose = FALSE,features = rownames(young.combined.sct))
set.seed(42069)
young.combined.sct <- RunUMAP(young.combined.sct, reduction = "pca", dims = 1:30)
young.combined.sct <- FindNeighbors(young.combined.sct, reduction = "pca", dims = 1:30)
young.combined.sct <- FindClusters(young.combined.sct, resolution = 0.4)

bcs <- colnames(young.combined.sct)
young11 <- grep(pattern = "_1$", x = bcs, value = TRUE)
young22 <- grep(pattern = "_2$", x = bcs, value = TRUE)
type <- data.frame(type=rep(NA,length(bcs)))
rownames(type) <- bcs
type[rownames(type) %in% young11,] = "Young1"
type[rownames(type) %in% young22,] = "Young2"
young.combined.sct[["type"]] <-as.vector(as.character(type$type))

x <- FeaturePlot(young.combined.sct, reduction = "umap",cols = c("grey",cbp2[1:2]),features = c("Cyp2e1","Cyp2f2"))&NoAxes()
k <- VlnPlot(young.combined.sct,cols = c("grey",cbp2[1:3]),features = c("Cyp2e1","Cyp2f2"))
c <- DimPlot(young.combined.sct, reduction = "umap",cols = cbp2)
j<- DimPlot(young.combined.sct, reduction = "umap",cols = cbp2[4:5],group.by = "type")&NoAxes()
x+k+c+j

j+x 

DefaultAssay(young.combined.sct) <- "RNA"
young.combined.markers <- list()
                    
new.cluster.ids <- c("Periportal young","Pericentral young","Periportal young","Pericentral young")
Idents(young.combined.sct) <- "seurat_clusters"
names(new.cluster.ids) <- levels(young.combined.sct)
young.combined.sct<- RenameIdents(young.combined.sct, new.cluster.ids)
young.combined.sct[["celltype"]] <- Idents(young.combined.sct)
VlnPlot(young.combined.sct,cols = c(cbp2[1:2]),features = c("Cyp2e1","Cyp2f2"),group.by = "celltype")

tt <- SpatialDimPlot(young.combined.sct,ncol=1)
zz <- SpatialFeaturePlot(young.combined.sct,features=c("Cyp2e1","Cyp2f2"))& scale_fill_gradientn(colours = c("grey","gold1","red"))
tt+zz

## Integrate old 
old.list <- list(old1,old2)

features.old <- SelectIntegrationFeatures(object.list =old.list , nfeatures = length(rownames(old1)))
old.anchors <- FindIntegrationAnchors(object.list = old.list, normalization.method = "LogNormalize", 
                                      anchor.features = features.old,reduction = "cca")
old.combined.sct <- IntegrateData(anchorset = old.anchors, normalization.method = "LogNormalize",features.to.integrate = consensus_genes)
DefaultAssay(old.combined.sct) <- "integrated"

old.combined.sct[["percent.mt"]] <- PercentageFeatureSet(old.combined.sct, pattern = "^mt-",assay = "RNA")

DefaultAssay(old.combined.sct) <- "integrated"
old.combined.sct <- ScaleData(old.combined.sct, verbose = FALSE,vars.to.regress = "percent.mt")
old.combined.sct <- RunPCA(old.combined.sct, verbose = FALSE,features = rownames(old.combined.sct))
set.seed(42069)
old.combined.sct <- RunUMAP(old.combined.sct, reduction = "pca", dims = 1:30)
old.combined.sct <- FindNeighbors(old.combined.sct, reduction = "pca", dims = 1:30)
old.combined.sct <- FindClusters(old.combined.sct, resolution = 0.5)

bcs <- colnames(old.combined.sct)
old11 <- grep(pattern = "_1$", x = bcs, value = TRUE)
old22 <- grep(pattern = "_2$", x = bcs, value = TRUE)
type <- data.frame(type=rep(NA,length(bcs)))
rownames(type) <- bcs
type[rownames(type) %in% old11,] = "old1"
type[rownames(type) %in% old22,] = "old2"
old.combined.sct[["type"]] <-as.vector(as.character(type$type))

x <- FeaturePlot(old.combined.sct, reduction = "umap",cols = c("grey",cbp2[1:2]),features = c("Cyp2e1","Cyp2f2"))&NoAxes()
k <- VlnPlot(old.combined.sct,cols = c(cbp2),features = c("Cyp2e1","Cyp2f2"))
c <- DimPlot(old.combined.sct, reduction = "umap",cols = cbp2)
j<- DimPlot(old.combined.sct, reduction = "umap",cols = cbp2[3:4],group.by = "type")&NoAxes()
x+k+c+j

j+x 

DefaultAssay(old.combined.sct) <- "RNA"
old.combined.markers <- list()
new.cluster.ids <- c("Pericentral old","Periportal old","Periportal old","Periportal old")
Idents(old.combined.sct) <- "seurat_clusters"
names(new.cluster.ids) <- levels(old.combined.sct)
old.combined.sct<- RenameIdents(old.combined.sct, new.cluster.ids)
old.combined.sct[["celltype"]] <- Idents(old.combined.sct)
VlnPlot(old.combined.sct,cols = c(cbp2[1:2]),features = c("Cyp2e1","Cyp2f2"),group.by = "celltype")

tt <- SpatialDimPlot(old.combined.sct,ncol=1)
zz <- SpatialFeaturePlot(old.combined.sct,features=c("Cyp2e1","Cyp2f2"))& scale_fill_gradientn(colours = c("grey","gold1","red"))
tt
zz

## merge
spatial.full <- merge(young.combined.sct, y = old.combined.sct, add.cell.ids = c("young","old"), project = "allspatialmerged",merge.data=T)
spatial.full
spatial.full <- subset(spatial.full,  subset = nFeature_RNA > 1250 & nFeature_RNA < 7000)

DefaultAssay(spatial.full) <- "integrated"

##getting tissue type meta data from the data.
bcs <- rownames(spatial.full@meta.data)
young11 <- grep(pattern = "^young_*", x = bcs, value = TRUE)
old11 <- grep(pattern = "^old_*", x = bcs, value = TRUE)

type <- data.frame(type=rep(NA,length(bcs)))
rownames(type) <- bcs
type[rownames(type) %in% young11,] = "Young"
type[rownames(type) %in% old11,] = "Old"
spatial.full[["type_merge"]] <-as.vector(as.character(type$type))
plot1 <- VlnPlot(spatial.full, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 2,group.by = "type",cols = cbp2)
plot2 <- FeatureScatter(spatial.full, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "type",cols = cbp2)
plot1+plot2

spatial.full <- ScaleData(spatial.full,assay = "integrated")
spatial.full <- RunPCA(spatial.full, features = rownames(spatial.full),assay = "integrated")
DimPlot(spatial.full, reduction = "pca",group.by = "type",cols = cbp2)
VizDimLoadings(spatial.full, dims = 1, reduction = "pca",nfeatures = 100) & theme(axis.text.y = element_text(size = 9))

spatial.full <- FindNeighbors(spatial.full, dims = 1:30)
spatial.full <- FindClusters(spatial.full)

set.seed(42069)
spatial.full <- RunUMAP(spatial.full, dims = 1:10)
c1 <- DimPlot(spatial.full, reduction = "umap",group.by = "celltype",cols = cbp2)+NoAxes()
c1



c2 <- FeaturePlot(spatial.full,features =c("Glul", "Cyp2e1", "Cyp2f2"),pt.size = .4,keep.scale = "all",cols=c("blue","grey","red"))  &theme(axis.title.x=element_blank(),
                                                                                                                                      axis.text.x=element_blank(),
                                                                                                                                      axis.ticks.x=element_blank(),
                                                                                                                                      axis.title.y=element_blank(),
                                                                                                                           axis.text.y=element_blank(),axis.ticks.y=element_blank())

##figure 1C
c1+c2

##figure 1D and F
DefaultAssay(spatial.full) <- "RNA"
plan("multisession", workers = 4)
spatial.full <- NormalizeData(spatial.full)
plan("multisession", workers = 4)
spatial.full <- ScaleData(spatial.full)
spatial.full <- FindVariableFeatures(spatial.full)
gc()

Idents(spatial.full) <- factor(x = Idents(spatial.full), levels = c("Pericentral young","Pericentral old","Periportal young","Periportal old"))

##figure 4B 
VlnPlot(spatial.full,features = c("Cidea","Cideb","Cidec"),cols = c(cbp2),assay = "RNA",group.by = "celltype")

spatial.full <- SetIdent(spatial.full, value = "celltype")

plan("multisession", workers = 6)
full.pp.y.vs.o <- FindMarkers(spatial.full,assay = "RNA",test.use = "MAST",only.pos = F,ident.1 = "Periportal young",ident.2 = "Periportal old",latent.vars = c("type","celltype"))
gc()
plan("multisession", workers = 6)
full.pc.y.vs.o <- FindMarkers(spatial.full,assay = "RNA",test.use = "MAST",only.pos = F,ident.1 = "Pericentral young",ident.2 = "Pericentral old",latent.vars = c("type","celltype"))
gc()


spatial.split_PC<- subset(spatial.full,  subset = nFeature_RNA > 1250 & nFeature_RNA < 7000,  idents = c("Pericentral young" , "Pericentral old"))

log2fcall_PC <- Seurat::FoldChange(spatial.split_PC,ident.1 = "Pericentral old",ident.2 = "Pericentral young",slot="data")

log2fcall_PC$p_val_adj <- rep(1,length(nrow(log2fcall_PC)))

avg_gene_expression <- log2(rowSums(GetAssayData(spatial.split_PC, "counts", "RNA"))+1)

avg_gene_expression <- apply(GetAssayData(spatial.split_PC, "data", "RNA"),MARGIN = 1,FUN = mean)

full.pc.y.vs.o$avg_gene_expression <-avg_gene_expression[rownames(full.pc.y.vs.o)]

log2fcall_PC$avg_gene_expression <- avg_gene_expression[rownames(log2fcall_PC)]

log2fcall_PC  <- log2fcall_PC[,c("avg_log2FC","p_val_adj","avg_gene_expression")]

full.pc.y.vs.o <- full.pc.y.vs.o[,c("avg_log2FC","p_val_adj","avg_gene_expression")]

MA_data_PC <- rbind(full.pc.y.vs.o,log2fcall_PC)
MA_data_PC$p_val_adj_t <- rep("Not significant",nrow(MA_data_PC))
MA_data_PC$p_val_adj_t <- ifelse(MA_data_PC$p_val_adj<0.05 & MA_data_PC$avg_log2FC >0,"Significantly upregulated",ifelse(MA_data_PC$p_val_adj<0.05 & MA_data_PC$avg_log2FC < 0,"Significantly downregulated","Not significant"))

MA_data_PC_sig <- MA_data_PC[which(MA_data_PC$p_val_adj<0.05),]
order <- sort(MA_data_PC_sig$avg_log2FC,decreasing = T,index.return=T)$ix
MA_data_PC$tops  <- rep(NA,nrow(MA_data_PC))

MA_data_PC[rownames(MA_data_PC_sig)[head(order,10)],"tops"] <- rownames(MA_data_PC_sig)[head(order,10)]
MA_data_PC[rownames(MA_data_PC_sig)[tail(order,10)],"tops"] <- rownames(MA_data_PC_sig)[tail(order,10)]

MA_data_PC$avg_log2FC <- MA_data_PC$avg_log2FC*-1
MA_data_PC$p_val_adj_t <- ifelse(MA_data_PC$p_val_adj<0.05 & MA_data_PC$avg_log2FC >0,"Significantly upregulated",ifelse(MA_data_PC$p_val_adj<0.05 & MA_data_PC$avg_log2FC < 0,"Significantly downregulated","Not significant"))

###Figure 1D
a <- MA_data_PC %>%
  ggplot() +
  aes(x=avg_gene_expression, y=avg_log2FC, col=p_val_adj_t,label=tops) +
  geom_point()+
  scale_color_manual(name = "",values=c("grey","dodgerblue2", "#E31A1C"))+
  geom_text_repel(min.segment.length = 0,box.padding = 0.9, max.overlaps = 90,
                  size=4,aes(x=avg_gene_expression, y=avg_log2FC, col=p_val_adj_t,label=tops),show.legend = F)+
  ggtitle("MA plot Pericentral young vs old")+
  scale_x_continuous(name="mean of normalized counts") +
  scale_y_continuous(name="log2 fold change")+
  theme_classic()
a 


spatial.split_PP<- subset(spatial.full,  subset = nFeature_RNA > 1250 & nFeature_RNA < 7000,  idents = c("Periportal young" , "Periportal old"))

log2fcall_PP <- Seurat::FoldChange(spatial.split_PP,ident.1 = "Periportal old",ident.2 = "Periportal young",slot="data")

log2fcall_PP$p_val_adj <- rep(1,length(nrow(log2fcall_PP)))

avg_gene_expression <- log2(rowSums(GetAssayData(spatial.split_PP, "counts", "RNA"))+1)

avg_gene_expression <- apply(GetAssayData(spatial.split_PP, "data", "RNA"),MARGIN = 1,FUN = mean)

full.pp.y.vs.o$avg_gene_expression <-avg_gene_expression[rownames(full.pp.y.vs.o)]

log2fcall_PP$avg_gene_expression <- avg_gene_expression[rownames(log2fcall_PP)]

log2fcall_PP  <- log2fcall_PP[,c("avg_log2FC","p_val_adj","avg_gene_expression")]

full.pp.y.vs.o <- full.pp.y.vs.o[,c("avg_log2FC","p_val_adj","avg_gene_expression")]

MA_data_PP <- rbind(full.pp.y.vs.o,log2fcall_PP)
MA_data_PP$p_val_adj_t <- rep("Not significant",nrow(MA_data_PP))
MA_data_PP$p_val_adj_t <- ifelse(MA_data_PP$p_val_adj<0.05 & MA_data_PP$avg_log2FC >0,"Significantly upregulated",ifelse(MA_data_PP$p_val_adj<0.05 & MA_data_PP$avg_log2FC < 0,"Significantly downregulated","Not significant"))

MA_data_PP_sig <- MA_data_PP[which(MA_data_PP$p_val_adj<0.05),]
order <- sort(MA_data_PP_sig$avg_log2FC,decreasing = T,index.return=T)$ix
MA_data_PP$tops  <- rep(NA,nrow(MA_data_PP))

MA_data_PP[rownames(MA_data_PP_sig)[head(order,10)],"tops"] <- rownames(MA_data_PP_sig)[head(order,10)]
MA_data_PP[rownames(MA_data_PP_sig)[tail(order,10)],"tops"] <- rownames(MA_data_PP_sig)[tail(order,10)]

MA_data_PP$avg_log2FC <- MA_data_PP$avg_log2FC*-1
MA_data_PP$p_val_adj_t <- ifelse(MA_data_PP$p_val_adj<0.05 & MA_data_PP$avg_log2FC >0,"Significantly upregulated",ifelse(MA_data_PP$p_val_adj<0.05 & MA_data_PP$avg_log2FC < 0,"Significantly downregulated","Not significant"))

###Figure 1F
a <- MA_data_PP %>%
  ggplot() +
  aes(x=avg_gene_expression, y=avg_log2FC, col=p_val_adj_t,label=tops) +
  geom_point()+
  scale_color_manual(name = "",values=c("grey","dodgerblue2", "#E31A1C"))+
  geom_text_repel(min.segment.length = 0,box.padding = 0.9, max.overlaps = 90,
                  size=4,aes(x=avg_gene_expression, y=avg_log2FC, col=p_val_adj_t,label=tops),show.legend = F)+
  ggtitle("MA plot Pericentral young vs old")+
  scale_x_continuous(name="mean of normalized counts") +
  scale_y_continuous(name="log2 fold change")+
  theme_classic()
a 

#### supplemental figures S1
###Figure S1E
a <- SpatialFeaturePlot(young2, features = c("Glul"),pt.size.factor = 1.8,ncol=1,alpha = 0) & scale_fill_gradientn(colours = c("grey","gold1","red"), breaks = c(0.0,1.0,2.0,3.0,4.0,5.0), limits=c(0,5)) & 
  theme( legend.title=element_text(size=13) )

b <- SpatialFeaturePlot(young2, features = c("Glul","Cyp2e1","Cyp2f2"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )
b+a

c <- SpatialFeaturePlot(old1, features = c("Glul"),pt.size.factor = 1.8,ncol=1,alpha = 0) & scale_fill_gradientn(colours = c("grey","gold1","red"), breaks = c(0.0,1.0,2.0,3.0,4.0,5.0), limits=c(0,5)) & 
  theme( legend.title=element_text(size=13) )

d <- SpatialFeaturePlot(old1, features = c("Glul","Cyp2e1","Cyp2f2"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )
d+c

###Figure S1F+G
all.merged <- merge(young1,c(young2,old1,old2), add.cell.ids = c("young_1","wawawa","old_1","wowowo"), project = "allspatialmerged",merge.data=T)
bcs <- rownames(all.merged@meta.data)
young11 <- grep(pattern = "^young_1*", x = bcs, value = TRUE)
young22 <- grep(pattern = "^wawawa", x = bcs, value = TRUE)
old11 <- grep(pattern = "^old_1*", x = bcs, value = TRUE)
old22 <- grep(pattern = "^wowowo", x = bcs, value = TRUE)

type <- data.frame(type=rep(NA,length(bcs)))
rownames(type) <- bcs
type[rownames(type) %in% young11,] = "Young 1"
type[rownames(type) %in% young22,] = "Young 2"
type[rownames(type) %in% old11,] = "Old 1"
type[rownames(type) %in% old22,] = "Old 2"
all.merged[["type"]] <-as.vector(as.character(type$type))
all.merged[["percent.mt"]] <- PercentageFeatureSet(all.merged, pattern = "^mt-")

all.merged <- subset(all.merged, subset = nFeature_RNA > 200  & percent.mt < 15)

plan("multisession", workers = 4)
all.merged <- NormalizeData(all.merged)
plan("multisession", workers = 4)
all.merged <- ScaleData(all.merged)
all.merged <- FindVariableFeatures(all.merged)
gc()

plan("multisession", workers = 4)
all.merged <- RunPCA(all.merged, verbose = FALSE,features = VariableFeatures(all.merged))  

DimPlot(all.merged,group.by = "type",cols=cbp2,reduction = "pca")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

##Figure S1G 
VizDimLoadings(all.merged,dims = 1,nfeatures = 50)+theme_classic(base_family = "sans")

##Figure S6A 
spatial.full <- SetIdent(spatial.full, value = "celltype")
Idents(spatial.full) <- factor(x = Idents(spatial.full), levels = c("Periportal young","Periportal old","Pericentral young","Pericentral old"))
cowplot::plot_grid(plotlist = VlnPlot(spatial.full,features = c("Cidea","Cideb","Cidec"),pt.size = .4, cols=c("#E31A1C","dodgerblue2","#6A3D9A","green4"),split.by = "type",combine = F,group.by = "celltype"),nrow = 1,ncol = 3)

##Figure S9
kupffer_young <- SpatialFeaturePlot(young1, features = c("Cec4f","rf7"),pt.size.factor = 1.8,ncol=1,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red"), breaks = c(0.0,1.0,2.0,3.0,4.0,5.0), limits=c(0,5)) & 
  theme( legend.title=element_text(size=13) )

endothelial_young <- SpatialFeaturePlot(young1, features = c("Dpp4","Ushbp1"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )

stellate_young <- SpatialFeaturePlot(young1, features = c("Igfbp3","Stx2"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )

kupffer_old <- SpatialFeaturePlot(old2, features = c("Cec4f","rf7"),pt.size.factor = 1.8,ncol=1,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red"), breaks = c(0.0,1.0,2.0,3.0,4.0,5.0), limits=c(0,5)) & 
  theme( legend.title=element_text(size=13) )

endothelial_old <- SpatialFeaturePlot(old2, features = c("Dpp4","Ushbp1"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )

stellate_old <- SpatialFeaturePlot(old2, features = c("Igfbp3","Stx2"),pt.size.factor = 1.8,ncol=5,alpha = 1) & scale_fill_gradientn(colours = c("grey","gold1","red")) & 
  theme( legend.title=element_text(size=13) )






