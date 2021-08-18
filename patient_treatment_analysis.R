library(data.table)
library(ggplot2)
library(reshape)
library(corrplot)
library(Matrix)
library(Seurat)
library(SeuratData)
library(dplyr)
library(tidyverse)
library(cowplot)
library(biomaRt)
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(patchwork)
library(ggpubr)
library(comprehenr)

project<-"AML-Cell"
Sys.setenv(language="en")
if (project == "AML-Cell"){
  working_dir<-"C:/Users/jtao/work_dir/AML/data/AMLPaper"
  graph_dir<- "C:/Users/jtao/work_dir/AML/data/images/"
  setwd(working_dir)
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}


targets <- c("CD70", "CD33", "IL3RA","CD34", "CD38",
             "CD123", "GP67", "CD90", "THY1", "CD45RA")
filename <- "SCT_AML_328_707_bm_norm.rds"
# filename <- "SCT_AML_328.rds"
AML<- readRDS(paste(working_dir, filename, sep = "/"))




data<- c()
Idents(AML) <- AML$CellType
heatmap <- DoHeatmap(AML, 
                     features = targets,   
                     size = 3.5, angle = 30, slot = "data")+
  labs(fills = "Cluster", colours = "Expression") + 
  theme(legend.position= "bottom", 
        legend.text = element_text(size = 9))+ 
  scale_fill_distiller(palette = "RdBu")
heatmap


Seurat_splt <- SplitObject(AML, split.by = "orig.ident")
Seurat_splt
sobjD <- Seurat_splt$AML707B.D18
DefaultAssay(sobjD) <- "SCT"
Idents(sobjD) <- sobjD$CellType
heatmap <- DoHeatmap(sobjD, 
                     features = targets,   
                     size = 3.5, angle = 30, slot = "data")+
  labs(fills = "Cluster", colours = "Expression") + 
  theme(legend.position= "bottom", 
        legend.text = element_text(size = 9))+ 
  scale_fill_distiller(palette = "RdBu")


heatmap

AML <- FindNeighbors(AML)
AML <- FindClusters(AML, resolution = 0.6)

DimPlot(AML)


Seurat_splt$AML328.D0

#visualization by D0-D171
# orig_ident <- c("AML707B.D0", "AML707B.D18",
#                 "AML707B.D41", "AML707B.D97", "AML707B.D113")
orig_ident <- c("AML328.D0", "AML328.D29", "AML328.D113", "AML328.D171")
table(HSC$orig.ident)

Idents(AML) <- AML$CellType
HSC <- subset(AML, idents = c("HSC-like"))
Idents(HSC) <- HSC$orig.ident
table(HSC$orig.ident)
HSC$orig.ident <- factor(x = HSC$orig.ident, levels = orig_ident)
pv <- VlnPlot(HSC, "CD33", group.by = "orig.ident",
        ident = orig_ident, slot = "data") + 
  scale_x_discrete(labels=c("D0", "D29", "D131", "D171")) +
  labs(x = "Days after treatment")
pv
DefaultAssay(HSC) <- "SCT"
meta <- AML@meta.data
HSC_counts <- c()
HSC_counts
for (i in orig_ident){
  HSC_counts <- append(HSC_counts
    , nrow(meta %>% filter(orig.ident == i & CellType == "HSC-like")))
  
}
#=================================================================
#statistical analysis 
#=================================================================
marker <- "CD70"
D0 <- Seurat_splt$AML328.D0
D2 <- Seurat_splt$AML328.D29
Idents(D0) <- D0$CellType
Idents(D2) <- D2$CellType
HSCD0 <- subset(D0, ident= "HSC-like")
HSCD2 <- subset(D2, ident = "HSC-like")
cd70D0 <- as.data.frame(HSCD0@assays[["SCT"]]@counts[marker, ])
cd70D2 <- as.data.frame(HSCD2@assays[["SCT"]]@counts[marker, ])
names(cd70D0) <- "counts"
names(cd70D2) <- "counts"
df <- rbind(cd70D0, cd70D2)
# rownames <- unlist(lapply(rownames(df), substr(5, 13)))
df$day_ident <- to_vec(for (i in rownames(df)) if(str_detect(i, "D0")) substr(i, 5, 13) else substr(i, 5, 14))

p1 <- ggboxplot(df, x = "day_ident", y = "counts",
          color = "day_ident", palette=  c("#00AFBB", "#E7B800"),
          ylab = "RNA expression", xlab = "Days of treatment",
          main = marker) +
  stat_compare_means()
p1
plot_grid(p1)
wc <- wilcox.test(cd70D0$counts, cd70D2$counts, alternative = "less")
wc

# meta %>% filter(orig.ident ==i) %>% filter(CellType == "HSC-like")
HSC_counts
ident_cnt <- as.data.frame(table(AML$orig.ident))
ident_cnt
Cell_counts <- ident_cnt[ident_cnt$Var1 %in% orig_ident, "Freq"]
Cell_counts
df <- data.frame(HSC_counts, Cell_counts)
colnames(df) <- c("HSC_counts", "Cell_counts")
rownames(df) <- orig_ident

png(paste(graph_dir, "cell summary-328.png", sep = ''), width = 400,
          height = 350)
barplot(as.matrix(t(df)), names = orig_ident, col = terrain.colors(5), ylab = "cell counts",
        xlab = "Days of treatment", width=0.3, font = 1,)
legend("topright", legend =  c("HSC-like cells", "All cells"), 
       fill = terrain.colors(5))
dev.off()
# ggplot(df, aes(x = rownames, y = HSC_counts)) + 
  # geom_bar(position = "dodge", stat = "identity")

table(AML$predicted.celltype.l2)






      