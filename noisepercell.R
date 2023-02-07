pkgs <- c("tidyverse","ggsignif","Seurat","BASiCS","SingleCellExperiment","ggridges")


paks<-lapply(pkgs, function(x) suppressMessages(require(x, character.only = TRUE)))
rm(paks)


reads.heps <- readRDS("seuratobj.rds")
meta.heps <- reads.heps@meta.data
meta.heps$celltype <- gsub("_lowcov","",meta.heps$celltype)
meta.heps$celltype.age <- gsub("_lowcov","",meta.heps$celltype.age)


pc.y.cells <- meta.heps %>% 
  filter(celltype.age == "Pericentral_hepatocytes_Young") %>% 
  arrange(.,Replicate)
pc.o.cells <- meta.heps %>% 
  filter(celltype.age == "Pericentral_hepatocytes_Old") %>% 
  arrange(.,Replicate)
pp.y.cells <- meta.heps %>% 
  filter(celltype.age == "Periportal_hepatocytes_Young") %>% 
  arrange(.,Replicate)
pp.o.cells <- meta.heps %>% 
  filter(celltype.age == "Periportal_hepatocytes_Old") %>% 
  arrange(.,Replicate)
y.cells <- meta.heps %>% 
  filter(Age == "Young" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)
o.cells <- meta.heps %>% 
  filter(Age == "Old" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)
y.2n.cells <- meta.heps %>% 
  filter(Age == "Young" & nuclei == "2n" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)
o.2n.cells <- meta.heps %>% 
  filter(Age == "Old" & nuclei == "2n" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)
y.4n.cells <- meta.heps %>% 
  filter(Age == "Young" & nuclei == "4n" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)
o.4n.cells <- meta.heps %>% 
  filter(Age == "Old" & nuclei == "4n" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)
y.8n.cells <- meta.heps %>% 
  filter(Age == "Young" & nuclei == "8n" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)
o.8n.cells <- meta.heps %>% 
  filter(Age == "Old" & nuclei == "8n" & celltype != "Macrophages" & !grepl("lowcov",celltype)) %>% 
  arrange(.,Replicate)

pc.y <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,pc.y.cells$RG]
pc.o <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,pc.o.cells$RG]
pp.y <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,pp.y.cells$RG]
pp.o <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,pp.o.cells$RG]

y <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,y.cells$RG]
o <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,o.cells$RG]

y.2n <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,y.2n.cells$RG]
o.2n <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,o.2n.cells$RG]
y.4n <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,y.4n.cells$RG]
o.4n <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,o.4n.cells$RG]
y.8n <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,y.8n.cells$RG]
o.8n <- GetAssayData(reads.heps, slot = "data", assay="RNA")[,o.8n.cells$RG]


keep_feature.yo <- intersect(rownames(y[rowSums(y)>0,]),rownames(o[rowSums(o)>0,]))
keep_feature.pc <- intersect(rownames(pc.y[rowSums(pc.y)>0,]),rownames(pc.o[rowSums(pc.o)>0,]))
keep_feature.pp <- intersect(rownames(pp.y[rowSums(pp.y)>0,]),rownames(pp.o[rowSums(pp.o)>0,]))
keep_feature.2n <- intersect(rownames(y.2n[rowSums(y.2n)>0,]),rownames(o.2n[rowSums(o.2n)>0,]))
keep_feature.4n <- intersect(rownames(y.4n[rowSums(y.4n)>0,]),rownames(o.4n[rowSums(o.4n)>0,]))
keep_feature.8n <- intersect(rownames(y.8n[rowSums(y.8n)>0,]),rownames(o.8n[rowSums(o.8n)>0,]))


o.filt <- o[keep_feature.yo,]
y.filt <- y[keep_feature.yo,]
pc.o.filt <- pc.o[keep_feature.pc,]
pc.y.filt <- pc.y[keep_feature.pc,]
pp.o.filt <- pp.o[keep_feature.pp,]
pp.y.filt <- pp.y[keep_feature.pp,]
o.2n.filt <- o.2n[keep_feature.2n,]
y.2n.filt <- y.2n[keep_feature.2n,]
o.4n.filt <- o.4n[keep_feature.4n,]
y.4n.filt <- y.4n[keep_feature.4n,]
o.8n.filt <- o.8n[keep_feature.8n,]
y.8n.filt <- y.8n[keep_feature.8n,]

cor.df <- NULL
cv.df <- NULL
norm.counts <- list()
for(i in ls(pattern=".filt")){
  a <- get(i)
  a <- NormalizeData(a, normalization.method = "LogNormalize", scale.factor = 10000)
  norm.counts[[i]] <- a
  cv <- rowSds(a)/rowMeans2(a)
  cv.df <- bind_rows(cv.df, data.frame(r = cv, GeneName = rownames(a), n = i))
  
  rowmeds <- rowMedians(as.matrix(a))
  cor.p <- apply(a, MARGIN = 2, function(x) cor(rowmeds, x, method = "pearson"))
  
  cor.df <- bind_rows(cor.df, data.frame(r = cor.p, n = i, cellID = names(cor.p)))
  
}
head(cor.df)
head(cv.df)
gene.exp.o <- data.frame(glul = norm.counts$o.filt["Glul",],
                       cyp2e1 = norm.counts$o.filt["Cyp2e1",],
                       cyp2f2 = norm.counts$o.filt["Cyp2f2",],
                       cellID = colnames(norm.counts$o.filt),
                       Age="Old")
gene.exp.y <- data.frame(glul = norm.counts$y.filt["Glul",],
                         cyp2e1 = norm.counts$y.filt["Cyp2e1",],
                         cyp2f2 = norm.counts$y.filt["Cyp2f2",],
                         cellID = colnames(norm.counts$y.filt),
                         Age="Young")
gene.exp <- bind_rows(gene.exp.o,gene.exp.y)

exp.cor <- cor.df %>% 
  filter(n %in% c("o.filt","y.filt")) %>% 
  mutate(Age = ifelse(grepl("o",n), "Old", "Young")) %>% 
  left_join(.,gene.exp,by=c("Age","cellID")) %>% 
  gather(.,key="Genes", value = "norm.exp", -r,-n,-cellID,-Age)

ggplot(exp.cor,aes(norm.exp,r,color=Age))+
  geom_point()+geom_smooth()+
  facet_wrap(~Genes,scales="fixed",ncol=1)

cor.df <- cor.df %>% 
  mutate(Age = ifelse(grepl("o",n), "Old", "Young"),
         nuclei = ifelse(grepl("[0-9]n",n), substr(n,3,4), NA),
         Zonation = ifelse(grepl("pc", n), "Pericentral", 
                           ifelse(grepl("pp",n),"Periportal", NA)))
table(cor.df$n,cor.df$Age)
table(cor.df$n,cor.df$Zonation)
table(cor.df$n,cor.df$nuclei)
cor.df %>% 
  filter(n %in% c("o.filt", "y.filt")) %>% 
  ggplot(.,aes(Age,r,fill=Age))+
  geom_boxplot(position=position_dodge(width=0.9),notch=T)+
  xlab("")+ylab("Pearson's R")+
  geom_signif(comparisons = list(c("Old","Young")))+
  scale_fill_manual(values = c("#800040","#1B81E5"))+
  theme_bw()+
  theme(text = element_text(size =18))

cor.df %>% 
  filter(!is.na(nuclei)) %>%
  ggplot(.,aes(Age,r,fill=Age))+
  geom_boxplot(position=position_dodge(width=0.9),notch=T)+
  xlab("")+ylab("Pearson's R")+
  facet_wrap(~nuclei,scales="free")+
  geom_signif(comparisons = list(c("Old","Young")))+
  scale_fill_manual(values = c("#800040","#1B81E5"))+
  theme_bw()+
  theme(text = element_text(size =18))

cor.df %>% 
  filter(!is.na(Zonation)) %>% 
  ggplot(.,aes(Age,r,fill=Age))+
  geom_boxplot(position=position_dodge(width=0.9),notch=T)+
  xlab("")+ylab("Pearson's R")+
  facet_wrap(~Zonation,scales="free")+
  geom_signif(comparisons = list(c("Old","Young")))+
  scale_fill_manual(values = c("#800040","#1B81E5"))+
  theme_bw()+
  theme(text = element_text(size =18),
        legend.position = "top")


cv.df <- cv.df %>% 
  mutate(Age = ifelse(grepl("o",n), "Old", "Young"),
         nuclei = ifelse(grepl("[0-9]n",n), substr(n,3,4), NA),
         Zonation = ifelse(grepl("pc", n), "Pericentral", 
                           ifelse(grepl("pp",n),"Periportal", NA)))

table(cv.df$n,cv.df$Age)
table(cv.df$n,cv.df$Zonation)
table(cv.df$n,cv.df$nuclei)

cv.df %>% 
  filter(n %in% c("o.filt", "y.filt")) %>% 
  ggplot(.,aes(Age,r,fill=Age))+
  geom_boxplot(position=position_dodge(width=0.9),notch=T)+
  xlab("")+ylab("Coefficient of variation")+
  geom_signif(comparisons = list(c("Old","Young")))+
  scale_fill_manual(values = c("#800040","#1B81E5"))+
  theme_bw()+
  theme(text = element_text(size =18))


cv.df %>% 
  filter(n %in% c("o.filt", "y.filt")) %>%
  select(GeneName,Age,r) %>% 
  spread(., key = "Age", value = "r") %>%
  ggplot(.,aes(Old,Young))+
  ggpointdensity::geom_pointdensity()+
  geom_abline()+
  theme_bw()+
  theme(text = element_text(size =18))

cv.df %>% 
  filter(!is.na(nuclei)) %>%
  ggplot(.,aes(Age,r,fill=Age))+
  geom_boxplot(position=position_dodge(width=0.9),notch=T)+
  xlab("")+ylab("Coefficient of variation")+
  facet_wrap(~nuclei,scales="free")+
  geom_signif(comparisons = list(c("Old","Young")))+
  scale_fill_manual(values = c("#800040","#1B81E5"))+
  theme_bw()+
  theme(text = element_text(size =18))

cv.df %>% 
  filter(!is.na(Zonation)) %>% 
  ggplot(.,aes(Age,r,fill=Age))+
  geom_boxplot(position=position_dodge(width=0.9),notch=T)+
  xlab("")+ylab("Coefficient of variation")+
  facet_wrap(~Zonation,scales="free")+
  geom_signif(comparisons = list(c("Old","Young")))+
  scale_fill_manual(values = c("#800040","#1B81E5"))+
  theme_bw()+
  theme(text = element_text(size =18),
        legend.position = "top")

cv.df %>% 
  filter(!is.na(nuclei)) %>%   
  ggplot(.,aes(x=r, y=Age,fill=Age))+
  geom_density_ridges(stat = "binline", bins = 10, scale=3, alpha=0.7)+
  stat_density_ridges() +
  xlab("")+ylab("Coefficient of variation")+
  facet_wrap(~nuclei,scales="fixed",ncol=1)+
  theme_bw()+
  theme(text = element_text(size =18))
cv.df %>% 
  filter(!is.na(Zonation)) %>%   
  ggplot(.,aes(x=r, y=Age,fill=Age))+
  geom_density_ridges(stat = "binline", bins = 10, scale=3, alpha=0.7)+
  #  stat_density_ridges() +
  xlab("")+ylab("Coefficient of variation")+
  facet_wrap(~Zonation,scales="fixed",ncol=1)+
  theme_bw()+
  theme(text = element_text(size =18))

cv.df %>% 
  filter(!is.na(nuclei)) %>% 
  select(GeneName,Age,nuclei,r) %>% 
  spread(., key = "Age", value = "r") %>%
  ggplot(.,aes(Old,Young))+
  ggpointdensity::geom_pointdensity()+
  geom_abline()+
  ggtitle("Coefficient of variation")+
  facet_wrap(~nuclei,ncol=1,scales="free")+
  theme_bw()+
  theme(text = element_text(size =18))

cv.df %>%
  filter(!is.na(Zonation)) %>% 
  select(GeneName,Zonation,Age,r) %>% 
  spread(., key = "Age", value = "r") %>%
  ggplot(.,aes(Old,Young))+
  ggpointdensity::geom_pointdensity()+
  geom_abline()+
  facet_wrap(~Zonation,ncol=1)+
  theme_bw()+
  theme(text = element_text(size =18))


