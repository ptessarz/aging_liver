
setwd("~/Desktop//scAgeing/")
source("scripts/functions.R")
pkgs <- c("reshape2","tidyverse","DESeq2","limma","pheatmap","Seurat","ggthemes","ggsignif","scater","scran","Matrix")
paks<-lapply(pkgs, function(x) suppressMessages(require(x, character.only = TRUE)))
rm(paks)

ginfo <- read.table("../Machi/UMIseq_Machi/geneid_symbol.93.txt", header = F, stringsAsFactors = F, sep="\t")
colnames(ginfo) <- c("GeneID","GeneName")
ginfo$GeneName <- trimws((ginfo$GeneName))
ginfo$GeneID <- trimws((ginfo$GeneID))
ginfo <- ginfo %>% dplyr::group_by(GeneName) %>% dplyr::mutate(n = paste(toupper(GeneName),seq(1:length(GeneID)),sep="-"))

mtgenes <- read.table("MTgenes.txt",header=F)
mtgenes <- ginfo %>% dplyr::filter(GeneID %in% mtgenes$V1)

# Load dge data ---------------------------------------------------------------
y <- readRDS("final_scripts/Young_organoids.dgecounts.rds")
o <- readRDS("final_scripts/Old_organoids.dgecounts.rds")

inexp.y <- y$umicount$inex$all
inexp.o <- o$umicount$inex$all

lcs.y <- (head(sort(Matrix::colSums(inexp.y),decreasing=T),n=5000))
lcs.o <- (head(sort(Matrix::colSums(inexp.o),decreasing=T),n=5000))

inexp.y <- inexp.y[,names(lcs.y)]
inexp.o <- inexp.o[,names(lcs.o)]

colnames(inexp.y) <- paste("run2a3.Young_organoids",seq(1:ncol(inexp.y)),colnames(inexp.y),sep="_")
colnames(inexp.o) <- paste("run2a3.Old_organoids",seq(1:ncol(inexp.o)),colnames(inexp.o),sep="_")

inexp <- Seurat::RowMergeSparseMatrices(inexp.y,inexp.o)


# sparsity, mitopercent genes and UMIs per cell ------------------------------
# inexp1 <- readRDS("seurat_objects/run2and3merged.YoungOldorganoids.Seurat3.rds")
 inexp <- GetAssayData(org.seurat, slot = "counts")

g <- countGenes(inexp)
u <- countUMIs(inexp)
colnames(g) <- c("Counts","SampleID")
colnames(u) <- c("Counts","SampleID")

percent.mito <- data.frame(Counts=colSums(inexp[intersect(rownames(inexp),mtgenes$GeneID),])/colSums(inexp), SampleID=colnames(inexp))


gumtcounts <- bind_rows(g%>%dplyr::mutate(Type="Gene"),
                      u%>%dplyr::mutate(Type="UMI"),
                      percent.mito%>%dplyr::mutate(Type="PercentMito"))

gumtcounts <- tidyr::separate(gumtcounts,col=SampleID,into=c("Age"),sep="_",remove=F)
gumtspread <- gumtcounts %>% tidyr::spread(.,key=Type,value=Counts) %>% dplyr::mutate(SparsityIndex = Gene/UMI)

gumtspread %>% 
  tidyr::gather(.,key="Type",value="Counts",-SampleID, -Age) %>%
  ggplot(aes(x=Type,y=Counts,color=Age))+
  geom_violin()+geom_boxplot(width=0.2)+
  facet_wrap(~Type,scales="free",nrow=1)

ggplot(gumtspread, aes(x=UMI,y=SparsityIndex,color=Age))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  facet_wrap(~Age,scales="free")

ggplot(gumtspread, aes(x=SparsityIndex,color=Age))+geom_density()+geom_vline(xintercept = 0.7,linetype="dashed",col="blue")+geom_vline(xintercept = 0.6,linetype="dashed",col="blue")

threshold <- 0.5
gumtspread %>% dplyr::group_by(Age) %>% dplyr::filter(SparsityIndex<=threshold) %>% dplyr::count()
selectcells <- gumtspread %>% 
  dplyr::group_by(Age) %>% 
  dplyr::filter(SparsityIndex<=threshold & PercentMito < 0.1)

selectcells %>% 
  tidyr::gather(.,key="Type",value="Counts",-SampleID, -Age) %>%
  ggplot(aes(x=Type,y=Counts,color=Age))+
  geom_violin()+geom_boxplot(width=0.2)+
  facet_wrap(~Type,scales="free",nrow=1)

inexp <- inexp[,selectcells$SampleID]
inexp <- inexp[Matrix::rowSums(inexp) > 0,selectcells$SampleID]

dim(inexp)


# change rownames to gene names instead of ensgids ------------------------

r <- ginfo %>% dplyr::filter(GeneID %in% rownames(inexp))
inexp <- inexp[r$GeneID,]
rownames(inexp) <- r$GeneName

# QC check cells and features ----------------------------------------------------------
#R3.6.3 scater 1.14.6
#
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(inexp)))
# apply quality filters
sce <- scater::calculateQCMetrics(sce)
colnames(colData(sce))
colnames(rowData(sce))

#remove cells that are outliers in QC metrics
sce <- runPCA(sce, use_coldata = TRUE,
                      detect_outliers = TRUE)
table(sce$outlier)

libsize.drop <- scater::isOutlier(sce$total_counts, nmads=3, type="both", log=TRUE)
feature.drop <- scater::isOutlier(sce$total_features_by_counts, nmads=3, type="both", log=TRUE)

table(libsize.drop)
table(feature.drop)


## Some diagnostic plots to look at the data before filtering
plotReducedDim(sce, use_dimred="PCA_coldata")
plotHighestExprs(sce, exprs_values = "counts",n = 20)
plotExprsFreqVsMean(sce)
plotColData(sce, x = "total_counts", y = "total_features_by_counts")
plotRowData(sce, x = "n_cells_by_counts", y = "mean_counts")
sce <- normalize(sce)
plotExplanatoryVariables(sce)

sce.filt <- sce[,!(libsize.drop | sce$outlier | feature.drop)]

keep_feature <- nexprs(sce.filt, byrow=TRUE) >= 5
table(keep_feature)

remove_feature <- c("MALAT1-1",mtgenes$n)
keep_feature[intersect(names(keep_feature),remove_feature)] <- FALSE

# kick out features that do not pass thresholds
#remove all the mt genes and MALAT1 which is unecessarily highly expressed in all cells
dim(sce.filt)
sce.filt <- sce.filt[keep_feature,]
dim(sce.filt)

selectcells %>% 
  dplyr::filter(SampleID %in% colnames(sce.filt)) %>%
  tidyr::gather(.,key="Type",value="Counts",-SampleID, -Age) %>%
  ggplot(aes(x=Type,y=Counts,color=Age))+
  geom_violin()+geom_boxplot(width=0.2)+
  facet_wrap(~Type,scales="free",nrow=1)

## Some diagnostic plots to look at the data after filtering
plotReducedDim(sce.filt, use_dimred="PCA_coldata")
plotHighestExprs(sce.filt, exprs_values = "counts")
plotExprsFreqVsMean(sce.filt)
plotColData(sce.filt, x = "total_counts", y = "total_features_by_counts")
plotRowData(sce.filt, x = "n_cells_by_counts", y = "mean_counts")
sce.filt <- normalize(sce.filt)
plotExplanatoryVariables(sce.filt)


# save the filtered count matrix ------------------------------------------


inexp <- BiocGenerics::counts(sce.filt)
dim(inexp)

saveRDS(inexp, file = "seurat_objects/Organoids_cellbarcoderetained_YoungOld.DGEtable.rds")
