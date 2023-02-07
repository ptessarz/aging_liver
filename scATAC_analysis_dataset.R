#This code was used to analyse the scATAC data mentioned in the manuscript and it's corresponding figures.
#to increse readability to users code has been rephrased to use standard R style not tidyverse style.
#to analyse the II Dataset just adjust the folder paths accordingly. 
# 1.Library initialization -----------------------------------------------
library(cisTopic)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(R.utils)
library(rtracklayer)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(tidyverse)
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(grid)
library(pheatmap)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(Gviz)
library(monocle3)
library(cicero)
library(ggrepel)
set.seed(42069)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
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
# 2. Cistopic Object initialization and Analysis ------------------------------------------------
pathTo10x <-"/outs/" 
data_folder <- paste0(pathTo10x, 'filtered_peak_bc_matrix')
metrics <- paste0(pathTo10x, 'singlecell.csv')

cisTopicObject <- createcisTopicObjectFrom10Xmatrix(data_folder, metrics, project.name='young_old_tissue_datasetI')
bcs <- colnames(cisTopicObject@binary.count.matrix)
old <- grep(pattern = "-1$", x = bcs, value = TRUE)
young <- grep(pattern = "-2$", x = bcs, value = TRUE)
type <- data.frame(type=rep(NA,length(bcs)))
rownames(type) <- bcs
type[rownames(type) %in% young,] = "Young"
type[rownames(type) %in% old,] = "Old"
cisTopicObject <- addCellMetadata(cisTopicObject,cell.data = type)

cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:50), seed=42069, nCores=52, addModels=FALSE, tmp="/outs/cistopictemps")

#saveRDS(cisTopicObject,file = "cistopicfitcistopic32.RDS")

cisTopicObject <- selectModel(cisTopicObject, select=32)

cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=4269, method='Probability')

UMAP_cis<- as.data.frame(cisTopicObject@dr$cell$Umap)
UMAP_cis <- merge(UMAP_cis,type,by="row.names")
rownames(UMAP_cis) <- UMAP_cis[,"Row.names"]
UMAP_cis <- UMAP_cis[,-1]
UMAP_cis <- as.data.frame(UMAP_cis)

##Figure 3A
ggplot(UMAP_cis, aes(x=UMAP1, y=UMAP2, color=type))+ geom_point(size=.75,)+scale_color_manual(values=c("#1B81E5","#800040"))+
  theme_classic()+theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())

##Figure 3E
par(mar = c(1, 1, 1, 1))
par(mfrow=c(6,6))
plotFeatures(cisTopicObject, method='Umap', target='cell',
             topic_contr='Probability', colorBy=NULL, cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE,col.low = "white",col.mid = "darkblue", col.high = "darkred")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
cisTopicObject <- annotateRegions(cisTopicObject, txdb=txdb, annoDb='org.Mm.eg.db')

cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=FALSE)

cisTopicObject <- GREAT(cisTopicObject, genome='mm10', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)

##Figure3F
ontologyDotPlot(cisTopicObject, var.y='name', order.by='Binom_Adjp_BH')
#ggsave("3F_ontologydotplot_32topic.pdf", width = 297, height = 210, units = "mm")
pathToFeather <- "mm9-regions-9species.all_regions.mc9nr.feather" 
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz"
mm10Tomm9.chain <- "mm10Tomm9.over.chain"
download.file(url, destfile = paste0(mm10Tomm9.chain, ".gz"))
gunzip(paste0(mm10Tomm9.chain, ".gz"))
mm10Tomm9.chain  <- import.chain(mm10Tomm9.chain)

# Obtain liftOver dictionary (as list)
mm10_coord <- cisTopicObject@region.ranges
mm10_to_mm9_list <- liftOver(mm10_coord, mm10Tomm9.chain)

cisTopicObject <- binarizedcisTopicsToCtx(cisTopicObject, liftOver=mm10_to_mm9_list, genome='mm9')
cisTopicObject <- scoredRegionsToCtx(cisTopicObject, liftOver=mm10_to_mm9_list, genome='mm9')
cisTopicObject <- topicsRcisTarget(cisTopicObject, genome='mm9', pathToFeather, reduced_database=FALSE, nesThreshold=3, rocthr=0.005, maxRank=20000, nCores=1)

##generate tables for Figure 3G 
Topic10_motif_enr <- cisTopicObject@binarized.RcisTarget[[5]]
DT::datatable(Topic10_motif_enr[,-c("enrichedRegions", "TF_lowConf"), with=FALSE], escape = FALSE, filter="top", options=list(pageLength=5))

# 3. complementary Signac analysis --------------------------------------
counts <- Read10X_h5("/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

tissue_dat <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/outs/fragments.tsv.gz', min.cells = 1)

tissue_dat <- CreateSeuratObject(
  counts = tissue_dat,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(tissue_dat) <- annotations

# add age information into metadata
bcs <- colnames(tissue_dat)
young <- grep(pattern = "-1$", x = bcs, value = TRUE)
old <- grep(pattern = "-2$", x = bcs, value = TRUE)
type <- data.frame(type=rep(NA,length(bcs)))
rownames(type) <- bcs
type[rownames(type) %in% young,] = "Young"
type[rownames(type) %in% old,] = "Old"
tissue_dat[["type"]] <-as.vector(as.character(type$type))

tissue_dat <- NucleosomeSignal(object = tissue_dat)

tissue_dat$nucleosome_group <- ifelse(tissue_dat$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = tissue_dat, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

tissue_dat <- TSSEnrichment(tissue_dat, fast = FALSE)

tissue_dat$high.tss <- ifelse(tissue_dat$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(tissue_dat, group.by = 'high.tss') + NoLegend()

tissue_dat$pct_reads_in_peaks <- tissue_dat$peak_region_fragments / tissue_dat$passed_filters * 100
tissue_dat$blacklist_ratio <- tissue_dat$blacklist_region_fragments / tissue_dat$peak_region_fragments

VlnPlot(
  object = tissue_dat,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1, cols = cbp2,group.by = "type",
  ncol = 4
)

tissue_dat <- subset(tissue_dat, subset = peak_region_fragments > 1000 & peak_region_fragments < 35000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)

#remove cells from dataset which couldn't be unambigously identified in the downstream cell identification process.
toremove <- readRDS("outs/HM_integration_filtered_cells/remove.RDS")
tissue_dat <- subset(tissue_dat, cells= colnames(tissue_dat)[!colnames(tissue_dat) %in% toremove] )
tissue_dat

## Normalization and Dimred 
tissue_dat <- RunTFIDF(tissue_dat)
tissue_dat <- FindTopFeatures(tissue_dat, min.cutoff = 'q0')
tissue_dat <- RunSVD(
  object = tissue_dat,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

DepthCor(tissue_dat)

#non-linear Dimred
set.seed(1234)
DefaultAssay(tissue_dat) <- "peaks"

tissue_dat <- RunUMAP(object = tissue_dat, reduction = 'lsi', dims = 2:30)
##figure S3a
DimPlot(object = tissue_dat, label = F,group.by = "type",cols=cbp2,reduction = "umap") &theme(axis.title.x=element_blank(),
                                                                                                       axis.text.x=element_blank(),
                                                                                                       axis.ticks.x=element_blank(),
                                                                                                       axis.title.y=element_blank(),
                                                                                                       axis.text.y=element_blank(),
                                                                                                       axis.ticks.y=element_blank())


##import cistopic object to transfer topic model for dimensionality reduction. 
cisTopicObject <-  readRDS("cistopicfitcistopic32.RDS")

umap_cis <- cisTopicObject@dr$cell$Umap
umap_cis <-umap_cis[rownames(umap_cis)%in% colnames(tissue_dat),]

tissue_dat[["umap_cis"]] <- CreateDimReducObject(embeddings = umap_cis, key = "UMAP", assay = "peaks")


#generate ggplot object for celltype selection
tissue_dat <- FindNeighbors(object = tissue_dat, reduction = 'lsi', dims = 2:30)
tissue_dat <- FindClusters(object = tissue_dat, verbose = FALSE,resolution = c(0.5,0.6,0.7,0.8,0.9,1))
a <- DimPlot(object = tissue_dat, label = F,group.by = "type",cols=cbp2,reduction = "umap_cis") &theme(axis.title.x=element_blank(),
                                                                                                            axis.text.x=element_blank(),
                                                                                                            axis.ticks.x=element_blank(),
                                                                                                            axis.title.y=element_blank(),
                                                                                                            axis.text.y=element_blank(),
                                                                                                            axis.ticks.y=element_blank())
# compute gene activities
gene.activities <- GeneActivity(tissue_dat)

# add the gene activity matrix to the Seurat object as a new assay
tissue_dat[['RNA']] <- CreateAssayObject(counts = gene.activities)
tissue_dat <- NormalizeData(
  object = tissue_dat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(tissue_dat$nCount_RNA)
)

#generate marker gene lists
hepatocytes <- c("Afp","Acly","Alb","Apoa1","Asl","Ass1","Cyp2e1","Cyp2f2","G6pc","Glul","Mup3","Pck1")
kupffer <- c("Clec4f","Cd68","Adgre1","Irf7","Spic","Mnda","Ear2")
endothelial <- c("Bmp2", "C1qtnf1", "CD32b", "Dpp4", "F8", "Il1a", "LYVE-1", "Mmrn2", "Myf6", "Oit3", "Pcdh12", "Ushbp1","Pecam","Cldn5")
itocell <- c("Des","Acta2","Stx2","Gfap","Igfbp3","Syp")
liverprogenitor <- c("Ascl2","Krt19","Dlk1","Epcam","Klf5","Abcc1","Aldh18a1","Capg","Cdc20","Enah","Fhl2","Gnb5")
poly <- c("Bub1", "Ccnb1", "Ccnb2", "Cdk1")
liverstem <- c("Cd44","Cd105","Cd29","Gata4","Epcam")
bileduct <- c("Prom1", "Cdh1", "Spp1", "Sox9", "Hes1","Epcam","Cd24")
bcells <- c("Cd22","Iglc")
tcells <- c("Cd3d","Trac","Ltb","Cxcr6")

###Figure 3d
DefaultAssay(tissue_dat) <- "RNA"
FeaturePlot( 
  object = tissue_dat,
  features = hepatocytes,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  combine = T,
  cols=c("grey","blue","red"),
  reduction = "umap_cis"
)&theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#select cells per celltype individually based on feature plot gene activities 
endothelial_cells <- CellSelector(a)
immune_cells <- CellSelector(a)
hepatocytes_I<- CellSelector(a)
hepatocytes_II<- CellSelector(a)
ito<- CellSelector(a)
kupffer<- CellSelector(a)
liverprogenitor_I <- CellSelector(a)
liverprogenitor_II <- CellSelector(a)

#assign celltypes in Meta data
tissue_dat@meta.data[endothelial_cells,"celltypes"] <- "endothelial cells"
tissue_dat@meta.data[immune_cells,"celltypes"] <- "immune cells"
tissue_dat@meta.data[hepatocytes_I,"celltypes"] <- "hepatocytes"
tissue_dat@meta.data[hepatocytes_II,"celltypes"] <- "hepatocytes"
tissue_dat@meta.data[ito,"celltypes"] <- "Ito cells"
tissue_dat@meta.data[kupffer,"celltypes"] <- "Kupffer cells"
tissue_dat@meta.data[liverprogenitor_I,"celltypes"] <- "liver progenitor cells"
tissue_dat@meta.data[liverprogenitor_II,"celltypes"] <- "liver progenitor cells"

##figure 3b
DimPlot(object = tissue_dat, label = F,group.by = "celltypes",cols=cbp2,reduction = "umap_cis") &theme(axis.title.x=element_blank(),
                                                                                                  axis.text.x=element_blank(),
                                                                                                  axis.ticks.x=element_blank(),
                                                                                                  axis.title.y=element_blank(),
                                                                                                  axis.text.y=element_blank(),
                                                                                                  axis.ticks.y=element_blank())


##figure S3b
DimPlot(object = tissue_dat, label = F,group.by = "celltypes",cols=cbp2,reduction = "umap") &theme(axis.title.x=element_blank(),
                                                                                                       axis.text.x=element_blank(),
                                                                                                       axis.ticks.x=element_blank(),
                                                                                                       axis.title.y=element_blank(),
                                                                                                       axis.text.y=element_blank(),
                                                                                                       axis.ticks.y=element_blank())

##figure S5b
gradient_umap <- data.frame(
  umap1=embeddings_umap$UMAP_1,
  umap2=embeddings_umap$UMAP_2,
  ncount=tissue_dat@meta.data$nCount_peaks
)
ggplot(gradient_umap,pt.size = 1.25, aes(x=umap1, y=umap2, color=ncount))+scale_color_gradientn(colours = rainbow(5))+
  theme_classic()+ geom_point()



hepatocytes_young_pc<- CellSelector(a)
hepatocytes_old_pc<- CellSelector(a)
hepatocytes_young_pp<- CellSelector(a)
hepatocytes_old_pp<- CellSelector(a)

tissue_dat@meta.data[hepatocytes_young_pc,"celltypes"] <- "hepatocytes young pericentral"
tissue_dat@meta.data[hepatocytes_old_pc,"celltypes"] <- "hepatocytes old pericentral"
tissue_dat@meta.data[hepatocytes_young_pp,"celltypes"] <- "hepatocytes young periportal"
tissue_dat@meta.data[hepatocytes_old_pp,"celltypes"] <- "hepatocytes old periportal"

tissue_dat_hep<- subset(tissue_dat, idents = c("hepatocytes young pericentral" , "hepatocytes old pericentral",
                                               "hepatocytes young periportal","hepatocytes old periportal"))
DefaultAssay(tissue_dat) <- "peaks"
Idents(tissue_dat) <- "type"
da_peaks <- FindMarkers(
  object = tissue_dat,
  ident.1 ="Young", 
  ident.2 = "Old",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

##figure 3H
DefaultAssay(tissue_dat) <- "peaks"
CoveragePlot(tissue_dat,region = "chr4-60656466-60664395",group.by = "celltypes")
CoveragePlot(tissue_dat,region = "chr7-30940368-30946017",group.by = "celltypes")
CoveragePlot(tissue_dat,region = "chr12-104336486104347739",group.by = "celltypes")
CoveragePlot(tissue_dat,region = "chr9-106231455-106249954",group.by = "celltypes")

####figure 3C
DefaultAssay(tissue_dat) <- "RNA"
meta_data <-tissue_dat@meta.data[colnames(tissue_dat),]

# order of annotations/colors are defined here
ordered_meta_data <- meta_data[order(meta_data$celltypes), ]
#meta_data <-data.frame(type=tissue_dat@meta.data[,"type"],row.names = rownames(tissue_dat@meta.data))
ordered_meta_data <- data.frame(type=ordered_meta_data$celltypes,row.names = rownames(ordered_meta_data))

annotation_colors <- list("type"= c( "hepatocytes" = "#A52A2A",
                                     "Kupffer cells" = "#6A3D9A",
                                     "endothelial cells"= "#FF7F00",
                                     "immune cells"= "#E31A1C",
                                     "Ito cells"= "#006600",
                                     "liver progenitor cells" = "#B4FBB4"
))
# Expression data
genes_to_use <- c(hepatocytes,kupffer,endothelial,itocell,liverprogenitor,bcells,tcells)
'%ni%' <- Negate('%in%')

my_data <- as.matrix(tissue_dat[["RNA"]]@data[genes_to_use[genes_to_use %in% rownames(tissue_dat)],])

my_data <- my_data[, rownames(ordered_meta_data)]

col_fun = colorRampPalette(c("white","blue1","red"))(256)
rowind <- cumsum(c(length(hepatocytes),length(kupffer),length(endothelial),length(itocell),length(liverprogenitor),length(bcells)+length(tcells)))
colind <-  head(as.numeric(cumsum(table(ordered_meta_data$type))))  

markermap <- pheatmap(my_data,cluster_rows = F,
                      color = colorRampPalette(c("grey95","blue1","red"))(256)
                      ,cluster_cols = F,
                      annotation_col = ordered_meta_data,
                      annotation_colors =annotation_colors,
                      legend_labels="Normalized accessibility",
                      annotation_names_row=T,
                      annotation_names_col=F,
                      show_rownames=T,
                      show_colnames =F, 
                      gaps_row=rowind	,
                      gaps_col=colind	,
                      row_gap = unit(5000000, "points") , 
                      border_color="black",
                      row_title = c("hepatocytes","Kupffer cells", "endothelial cells","immune cells","Ito cells","liver progenitor cells")
)
markermap

save_pheatmap_pdf <- function(x, filename, width=21, height=21) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


# 4.Cicero CCann analysis --------------------------------------
tissue_dat_hep <- subset(tissue_dat,subset = celltypes == "Hepatocytes")

hep_young <- subset(tissue_dat_hep,subset = type == "Young")
hep_old <-  subset(tissue_dat_hep,subset = type == "Old")

# convert to CellDataSet format and make the cicero object
hep.young.cds <-  as.cell_data_set(x = hep_young)
hep.old.cds <- as.cell_data_set(x = hep_old)

young.cicero <- make_cicero_cds(hep.young.cds, reduced_coordinates = reducedDims(hep.young.cds)$UMAP)
old.cicero <- make_cicero_cds(hep.old.cds, reduced_coordinates = reducedDims(hep.old.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(tissue_dat_hep)

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns.young <- run_cicero(young.cicero, genomic_coords = genome.df, sample_num = 100)
conns.old <- run_cicero(old.cicero, genomic_coords = genome.df, sample_num = 100)

ccans.young <- generate_ccans(conns.young)
ccans.old <- generate_ccans(conns.old)

links.young <- ConnectionsToLinks(conns = conns.young, ccans = ccans.young)
Links(hep_young) <- links.young

links.old <- ConnectionsToLinks(conns = conns.old, ccans = ccans.old)
Links(hep_old) <- links.old

##generate figure 4D with Gviz
# download genome track and unzip
temp <- tempfile()
download.file("http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)

unlink(temp)
###import chipseq data 
y1 <-  import.bw("ChIP_H3K27ac-Y-1.norm.bw")
y2 <-  import.bw("ChIP_H3K27ac-Y-2.norm.bw")
o1 <-  import.bw("ChIP_H3K27ac-O-1.norm.bw")
o2 <-  import.bw("ChIP_H3K27ac-O-2.norm.bw")

y1_string <- GRangesToString(y1, sep = c("_", "_"))
y2_string <- GRangesToString(y2, sep = c("_", "_"))
o1_string <- GRangesToString(o1, sep = c("_", "_"))
o2_string <- GRangesToString(o2, sep = c("_", "_"))

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

y_conns <- data.frame(Peak1=y1_string,Peak2=y1_string,coaccess=rep(.001,length(y1_string)))
o_conns <- data.frame(Peak1=o1_string,Peak2=o1_string,coaccess=rep(.001,length(o1_string)))

chr <- "chr18"
cis_start <-67341564-10000 #add space left and right from gene
cis_end <-67369794+5000 

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gtTrack<- GeneRegionTrack(txdb,chromosome = chr,start= cis_start, end=cis_end,
                          transcriptAnnotation="GENEID", # symbol is the gene symbol
                          fontsize.group=20,col="#E1AD01" # free to adjust font size
                          
)

convertensembl <- function(x = gtTrack){
  require(biomaRt)
  require(org.Mm.eg.db)
  
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  convertedgene <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
                         filters = "entrezgene_id",
                         values = x@range@elementMetadata@listData$gene,
                         mart = mouse)
  
  
  x@range@elementMetadata@listData$gene<- convertedgene[ match(x@range@elementMetadata@listData$gene, convertedgene$entrezgene_id ) ,"external_gene_name"]
  
  return(x)
}
# convert ensembl id to gene name
gtTrack <- convertensembl(gtTrack)


y1_dat <- DataTrack(range=y1, genome="mm10", name="chIP-seq_3m",col="#1B81E5",col.line=c("#1B81E5"),lwd=4,chromosome = chr,start= cis_start, end=cis_end,ylim=c(0,200)) #ylim=c(0,200)
y2_dat <- DataTrack(range=y2, genome="mm10", name="chIP-seq_3m",col="#1B81E5",fill="#1B81E5",lwd=4,chromosome = chr,start= cis_start, end=cis_end,ylim=c(0,200))
o1_dat <- DataTrack(range=o1, genome="mm10", name="chIP-seq_12m",col="#800040",fill="#800040",col.line="#800040",lwd=4,chromosome = chr,start= cis_start, end=cis_end,ylim=c(0,200))
o2_dat <- DataTrack(range=o2, genome="mm10", name="chIP-seq_12m",col="#800040",fill="#800040",lwd=4,chromosome = chr,start= cis_start, end=cis_end,ylim=c(0,200))

plot_list_old <- plot_connections(conns.old,chr, cis_start, cis_end,
                                  gene_model = NULL,
                                  collapseTranscripts = "longest", 
                                  return_as_list = T,
                                  coaccess_cutoff = 0.15,include_axis_track = T,connection_width = .5,connection_color = "#800040",connection_ymax = 1)


plot_list_young <- plot_connections(conns.young,chr, cis_start, cis_end,
                                    gene_model = NULL,
                                    collapseTranscripts = "longest", 
                                    return_as_list = TRUE,coaccess_cutoff = 0.15,include_axis_track = F,connection_width = .5,connection_color = "#1B81E5",connection_ymax = 1)


plot_list_merged <- list(plot_list_young[[1]],plot_list_young[[2]],plot_list_old[[1]],plot_list_old[[2]])


###figure 4D for the full figure dataset II is to be processed the same way. 

plotTracks(c(plot_list_merged,y1_dat,o1_dat,gtTrack),type = c("l"), from =cis_start, to =cis_end, chromosome = chr,collapseTranscripts = "longest",
           transcriptAnnotation="gene",
           #type="hist",
           background.title = "white",
           sizes = c(1,.25,1,.25,1.5,1.5,.5),
           fontcolor = "black",
           col.axis="black",
           fontsize=15,
           showTitle=T,
           margin=40,
           innerMargin = 5
)


##figure 4E
refTSS_v3.3_mouse_coordinate.ann.mm10 <-  readRDS("/refTSS_v3.3_mouse_coordinate.ann.mm10.rds")
DE_genes_all_spatial <- read_excel("/DE_genes_all_spatial.xlsx")
names_de <- DE_genes_all_spatial$...1
DefaultAssay(tissue_dat) <- "peaks"
tissue_dat_hep <- subset(tissue_dat,subset = celltypes == "hepatocytes")
hep_young <- subset(tissue_dat_hep,subset = type == "Young")
hep_old <-  subset(tissue_dat_hep,subset = type == "Old")

#import merged spatial dataset
full.y.vs.o <- readRDS("~/Liver_Organoid_project/scRNA_spatial_tissue/ccans/full.y.vs.o.RDS")

DE_TSS_reg <- refTSS_v3.3_mouse_coordinate.ann.mm10[(elementMetadata(refTSS_v3.3_mouse_coordinate.ann.mm10)[,"Gene_symbol"] %in% names_de)]

conns.old.filt <- na.omit(conns.old)
conns.young.filt <- na.omit(conns.young)

reg.young <-  StringToGRanges(unique(conns.young.filt$Peak1),sep = c("-","-"))
reg.old <- StringToGRanges(unique(conns.old.filt$Peak1),sep = c("-","-"))

intersect.y <- IRanges::subsetByOverlaps(x=reg.young,ranges=DE_TSS_reg)
intersect.o <- IRanges::subsetByOverlaps(x=reg.old,ranges=DE_TSS_reg)

TSS_inters_DE_spatial_y <- GRangesToString(intersect.y)
TSS_inters_DE_spatial_o <- GRangesToString(intersect.o)

hits_y <- conns.young.filt[conns.young.filt$Peak1 %in% TSS_inters_DE_spatial_y,]
hits_o <- conns.old.filt[conns.old.filt$Peak1 %in% TSS_inters_DE_spatial_o,]

hits_y_filt <- hits_y[hits_y$coaccess>0.25,]
hits_o_filt <- hits_o[hits_o$coaccess>0.25,]

hits_y_count <- as.data.frame(table(hits_y_filt$Peak1))
names(hits_y_count)[1:2] <- c("region","count_young")

hits_o_count <- as.data.frame(table(hits_o_filt$Peak1))
names(hits_o_count)[1:2] <- c("region","count_old")

countmat <- merge(hits_y_count,hits_o_count,by.x="region",by.y="region",all.x=T,all.y=T)

countmat <- countmat %>% replace(is.na(.), 0)

countmat_diff <- abs(countmat$count_old-countmat$count_young)
names(countmat_diff) <- countmat$region

countmat_diff[order(countmat_diff,decreasing=T)]

countmat_diff_neg <- countmat$count_young-countmat$count_old
names(countmat_diff_neg) <- countmat$region

countmat_diff_neg_sort <- as.data.frame(countmat_diff_neg[order(countmat_diff_neg,decreasing=T)])
names(countmat_diff_neg_sort)[1] <- "connections"
countmat_diff_neg_sort$order <- c(1:length(countmat_diff_neg))
countmat_diff_neg_sort$regions <- names(countmat_diff_neg)[order(countmat_diff_neg,decreasing=T)]
gg <- ClosestFeature(tissue_dat_hep,regions=countmat_diff_neg_sort$regions)
countmat_diff_neg_sort$gene <- gg[ match(countmat_diff_neg_sort$regions, gg$query_region ) ,"gene_name"]

countmat_diff_neg_sort <-   countmat_diff_neg_sort[countmat_diff_neg_sort$gene %in% rownames(full.y.vs.o),]
countmat_diff_neg_sort$order <- c(1:length(countmat_diff_neg_sort$regions))
upslightly <- rownames(full.y.vs.o)[which(full.y.vs.o$avg_log2FC > 0 & full.y.vs.o$avg_log2FC < 0.5 )]
upmore  <- rownames(full.y.vs.o)[which(full.y.vs.o$avg_log2FC > 0.5 )]

downslightly <- rownames(full.y.vs.o)[which(full.y.vs.o$avg_log2FC < 0 & full.y.vs.o$avg_log2FC > -0.5 )]
downmore <- rownames(full.y.vs.o)[which(full.y.vs.o$avg_log2FC < -0.5)]
countmat_diff_neg_sort$logfc [countmat_diff_neg_sort$gene %in% upmore] <- "up > 0.5"
countmat_diff_neg_sort$logfc [countmat_diff_neg_sort$gene %in% upslightly] <- "up < 0.5"

countmat_diff_neg_sort$logfc [countmat_diff_neg_sort$gene %in% downmore] <- "down > -0.5"
countmat_diff_neg_sort$logfc [countmat_diff_neg_sort$gene %in% downslightly] <- "down < -0.5"

ggplot(data=countmat_diff_neg_sort, aes(x=order, y=connections, group=1)) +
  geom_segment( aes(x=order, xend=order, y=0, yend=connections, color=logfc)) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()) +
  geom_label_repel(aes(label = gene),
                   box.padding   = 0.45, 
                   point.padding = 0.5,max.overlaps = 15,
                   segment.color = "grey50")+
  xlab("") +
  ylab("number of connections")+
  scale_color_manual(values = c("deeppink4","brown1","aquamarine","blue1"))+
  theme_classic()














