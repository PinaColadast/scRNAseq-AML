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


#=====================================================================
#helper functions
#=====================================================================
my_DimPlot <- function(SeuratObject){
  pdim<- DimPlot(SeuratObject, reduction = "umap",
                 label = TRUE, 
                 # cols = MyColorsCellTypes,
                 pt.size = 1, 
                 repel = TRUE,
                 label.size = 4)+theme_void()+theme(legend.position="top")
  return (pdim)}

save_png <- function(plot, pngname, png_width, png_height){
  png(file = paste(graph_dir, pngname, ".png", sep=""),
      width=png_width, height=png_height, units="in", res=200)
  print(plot)
  dev.off()
}


DGE_pvalue_plot <- function(data2, group1, group2){
  main = paste("differential expression of ", group1, " vs ", group2,
               sep = "")
  
  png(paste(graph_dir, group1," ", group2,".wilcox-cluster.png", sep= ""), 
      width = 12, height = 8, units = "in", res = 500)
  
  p <- ggplot(data2,aes(x=gene,y=avg_log2FC, fill=significance)) + 
    geom_bar( stat = "identity", size=0.4, width = 0.9, position="dodge", show.legend = T) +
    labs(y="average LOG2FC", x = "Genes", title =main) +
    theme_classic()+ theme(axis.text.x = element_text(angle = 90, size = 12),
                           axis.title.y = element_text(size = 14),
                           legend.title = element_text(size = 14),
                           legend.text = element_text(size = 13) )+ 
    scale_fill_discrete() 
  
  print(p)
  dev.off()
  
}

RECODE_Ontology <- function(Seuratobj){
  
  Seuratobj$predicted.celltype.l2 <- recode(Seuratobj$predicted.celltype.l2, "CD14 Mono"="CD14-positive, CD16-negative classical monocyte",
                                            "CD16 Mono"="CD14-low, CD16-positive monocyte", 
                                            "CD4 Memory"="effector memory CD4 T cell", "CD56 bright NK"="natural killer cell",
                                            "CD8 Effector_1"="effector CD8 T cell", "CD8 Effector_2"="effector CD8 T cell", "CD8 Memory_2"="effector memory CD8 T cell", 
                                            "CD8 Memory_1"="effector memory CD8 T cell", "cDC2"="CD1c-positive myeloid dendritic cell",
                                            "GMP" ="granulocyte monocyte progenitor", "HSC"="hematopoietic stem cell","LMPP"="lymphoid-primed multipotent progenitor" , 
                                            "Memory B"="memory B cell", "NK"="natural killer cell", "pDC"="plasmacytoid dendritic cell", "Plasmablast"="plasmablast", 
                                            "Prog_B 1"="pro-B cell", "Prog_B 2"="pro-B cell","Prog_DC"="common dendritic progenitor", "Prog_Mk"="megakaryocyte progenitor cell", "Prog_RBC"="erythroid progenitor cell",
                                            "Naive B"="naive B cell","NULL"="NULL")
  
  return(Seuratobj)
}

add.mut.CellType <- function(SeuratObj, reference = TRUE){
  cell.ident <- c()
  
  if (reference == TRUE){
    for (i in 1:length(SeuratObj$predicted.celltype.l2)){
      
      CellType <- as.vector(SeuratObj$predicted.celltype.l2)
      mut <- SeuratObj$mutation_annotation
      
      if (mut[i] == "mutation detected"){
        cell.ident<- append(cell.ident, 
                            paste(CellType[i], "malignant",
                                  sep =  "-"))
      }else{cell.ident<-append(cell.ident, CellType[i])}
      
    }
    SeuratObj$l2CellType.mut <- cell.ident
    
    
  } else{
    for (i in 1:length(SeuratObj$CellType)){
      
      CellType <- as.vector(SeuratObj$CellType)
      mut <- SeuratObj$mutation_annotation
      
      if (mut[i] == "mutation detected"){
        cell.ident<- append(cell.ident, 
                            paste(CellType[i], "malignant",
                                  sep =  "-"))
      }else{cell.ident<-append(cell.ident, CellType[i])}
      
    }
    SeuratObj$CellType.mut <- cell.ident}
  
  return (SeuratObj)
}


#======================================================================
#analysis 
#======================================================================



AML_norm<- readRDS(paste(working_dir, "/SCT_AML_D0_bm_norm.rds", sep = ""))
table(AML_norm@meta.data[["predicted_celltype_l2"]])

AML_norm <- FindNeighbors(AML_norm)
AML_norm <- FindClusters(AML_norm, resolution = 0.6)

Idents(AML_norm) <- AML_norm$CellType
my_DimPlot(AML_norm)



table(AML_norm$CellType)
aa <- AML_norm@meta.data %>% filter(CellType== "cDC-like")
table(aa$orig.ident)

#differential gene expression analysis#



pair_norm <- c("HSC", "Prog", "Mono", "ProMono", "cDC", "GMP")
pair_tumor <- c("HSC-like", "Prog-like", "Mono-like", "ProMono-like",
                "cDC-like", "GMP-like")

length(pair_norm) == length(pair_tumor)


Idents(AML_norm) <-AML_norm$orig.ident
Idents(AML_norm) <- AML_norm$CellType

celltype.markers.list <- list()
my_DimPlot(AML_norm)
for (i in 1:(length(pair_norm))){
  celltype.markers.list[[i]] <- FindMarkers(AML_norm, ident.1 = pair_norm[i],
                                           ident.2 = pair_tumor[i],
                                           logfc.threshold = 0.25, test.use = "negbinom", only.pos = TRUE)
}
 


celltype.markers.list[[1]]

data <- t(as.matrix(celltype.markers.list[[1]][0:20, c(1,2)]))
data
for (i in 1:length(pair_norm)){
  data1 <- t(as.matrix(celltype.markers.list[[i]][0:20, c(1,2)]))
  DGE_barplot(data1, pair_norm[i], pair_tumor[i])
}

DGE_barplot <- function(data, group1, group2){
  main = paste("differential expression of ", group1, " vs ", group2,
               sep = "")
  png(paste(graph_dir, group1," ", group2,".wilcox-cluster.png", sep= ""), 
      width = 500, height = 300)
  
  barplot(data, names = colnames(data), main = main, col=colors()[c(23,33)],
          border = "white", font.lab = 2,
          beside = T, font.axis = 1, las =2, srt = 45, 
         )
  legend("topright", legend = c("p_val", "avg log2FC"), ncol =1, cex = 0.8,
         fill = colors()[c(23,33)])
  dev.off()

}


dev.off()

my_DimPlot(AML_norm)


#====================================================================
#expression of groups with markers on AML nature dataset 
#=====================================================================
AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_508084.anno.rds")
AML_nature2 <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_548327.anno.rds")

AML_nature <- FindNeighbors(AML_nature)
AML_nature <- FindClusters(AML_nature, resolution = 0.5, graph.name = "SCT_snn")

Idents(AML_nature) <- AML_nature$predicted.celltype.l2


max(c(AML_nature$seurat_clusters))


my_DimPlot(AML_nature)
Idents(AML_nature) <- AML_nature$seurat_clusters
cluster_no <- length(unique(AML_nature$seurat_clusters))

cluster.markers.list <- list()
for (i in 0:(cluster_no-1)){
  j = i+1
  cluster.markers.list[[j]] <- FindMarkers(AML_nature, ident.1 = i,
                                            ident.2 = NULL,
                                            logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
}


cluster.markers.list[[1]][1:20, ]




AML_nature  <- add.mut.CellType(AML_nature)
AML_nature2 <- add.mut.CellType(AML_nature2)

Idents(AML_nature) <- AML_nature$CellType.mut
FindMarkers(AML_nature, ident.1 = "HSC",
            ident.2 = "HSC-malignant",
            logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

Idents(AML_nature2)<- AML_nature2$predicted.celltype.l2
DefaultAssay(AML_nature) <- "RNA"
target_marker <- append(targets , c("NPM1", "THY1"))
heatmap <- DoHeatmap(AML_sub, 
                     features = Tmarkers,   
                     size = 3.5, angle = 30, slot = "data", combine = FALSE,
                     label = FALSE)+
  labs(fills = "Cluster", colours = "Expression") + 
  theme(legend.position= "bottom", 
        legend.text = element_text(size = 9))+ 
  scale_fill_distiller(palette = "RdBu") 
# "THY1" %in% rownames(AML_nature@assays[["RNA"]]@counts)
heatmap
trace(DoHeatmap, edit = TRUE)

cluster0 <- AML_nature@meta.data %>% filter(seurat_clusters == 0)
table(cluster0$predicted.celltype.l2)
# df_meta <-AML_nature2@meta.data
# 
# table(AML_nature$predicted.celltype.l2)
# 
# 
# hsc_id <- rownames(df_meta[df_meta[["CellType.mut"]] == "HSC", ])
# hscmag_id <- rownames(df_meta[df_meta[["CellType.mut"]]== "HSC-malignant", ])
# 
# hsc_npm1 <- as.vector(AML_nature2@assays[["SCT"]]@counts["CD33", hsc_id])
# hscmag_npm1 <- as.vector(AML_nature2@assays[["SCT"]]@counts["CD33", hscmag_id])
# 
# wilcox.test(hsc_npm1, hscmag_npm1, alternative = "less")
# 
# 
# df1 <- as.data.frame(hsc_npm1,
#                      col.names = c("counts"), make.names = TRUE)
# df1$ident <- rep(c("HSC"), length(hsc_npm1))
# names(df1) <- c("counts", "ident")
# df2 <- as.data.frame(hscmag_npm1,
#                      col.names = c("counts"), make.names = TRUE) 
# 
# df2$ident <- rep(c("HSC-malignant"),
#                              length(hscmag_npm1))
# names(df2) <- c("counts", "ident")
# 
# df_NPM1 <- rbind(df1, df2)
# 
# 
# p1 <- ggboxplot(df_NPM1, x = "ident", y = "counts",
#           color = "ident", palette=  c("#00AFBB", "#E7B800"),
#           ylab = "RNA expression", xlab = "Cell Type",
#           main = "NPM1") +
#   stat_compare_means()

table(AML_nature$mutated_gene)

cluster.markers.list[[5]][0:30, ]

for (i in (1:length(cluster.markers.list))){
  data_a <- cluster.markers.list[[i]]
  if (dim(data_a)[1] <= 20){
    data2 <- cluster.markers.list[[i]][, c(1,2)]
  }else{
    data2 <-cluster.markers.list[[i]][0:25, c(1,2)] }
  
  data2$significance <- to_vec(for (i in as.vector(data2$p_val)) if(i<=0.01) "p-value <= 0.01" else "p-value > 0.01")
  data2$gene <- rownames(data2)
  DGE_pvalue_plot(data2, paste0("cluster",i-1, sep = ""),
  "other clusters")
}

table(AML_nature$predicted.celltype.l2)



# AML_nature <- add.mut.CellType(AML_nature, reference = FALSE)

LMPP <- AML_nature@meta.data %>% filter(CellType == "LMPP" | CellType == "HSC")
table(LMPP$predicted.celltype.l2)  
LMPP
Idents(AML_nature) <- AML_nature$predicted.celltype.l2
AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
              "granulocyte monocyte progenitor", "hematopoietic stem cell",
              "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")
AML_sub <- subset(AML_nature, ident = AML_category)  
table(AML_nature$predicted.celltype.l2)


Idents(AML_sub) <- AML_sub$l2CellType.mut




plot_stat_gene <- function(SeuratObj, gene_name){
  
  
  AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                    "granulocyte monocyte progenitor", "hematopoietic stem cell",
                    "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")
  
  df_meta <- SeuratObj@meta.data
  plots <- list()
  
  for (i in 1:length(AML_category)){
    celltype <- AML_category[i]
    muttype <- paste(celltype, "malignant", sep = "-")
    hsc_id <- rownames(df_meta[df_meta[["l2CellType.mut"]] == celltype, ])
    hscmag_id <- rownames(df_meta[df_meta[["l2CellType.mut"]]== muttype, ])
    
    hsc_npm1 <- as.vector(SeuratObj@assays[["SCT"]]@counts[gene_name, hsc_id])
    hscmag_npm1 <- as.vector(SeuratObj@assays[["SCT"]]@counts[gene_name, hscmag_id])
    
    if (length(hsc_npm1) == 0|length(hscmag_npm1) ==0){
      next
    }
    # wilcox.test(hsc_npm1, hscmag_npm1, alternative = "two.sided")
    
    df1 <- as.data.frame(hsc_npm1,
                         col.names = c("counts"), make.names = TRUE)
    
    df1$ident <- celltype
    names(df1) <- c("counts", "ident")
    
    df2 <- as.data.frame(hscmag_npm1,
                         col.names = c("counts"), make.names = TRUE) 
    
    df2$ident <- muttype
    # print(df2)
    names(df2) <- c("counts", "ident")
  
    df_gene<- rbind(df1, df2)
    
    
    p1 <- ggboxplot(df_gene, x = "ident", y = "counts",
                    color = "ident", palette=  c("#00AFBB", "#E7B800"),
                    ylab = "RNA expression", xlab = "Cell Type",
                    main = gene_name, ) +
      stat_compare_means(family = "Times New Roman",
                         size = 2) + 
        scale_x_discrete(breaks = c(celltype, muttype),
                         label=addline_format(c(celltype, muttype)))
                    
    
    plots[[i]] <- p1
  }
  
  return(plots)
  
}

p.cd70 <- plot_stat_gene(AML_sub, "CD70")

plot_grid(plotlist = p.cd70[[1:2]])


df_counts <- as.data.frame(table(AML_sub$l2CellType.mut))
names <- df_counts %>% filter(Freq >3)
as.vector(names$Var1)

#boxplot

boxplot_stat_gene <- function(SeuratObj, gene_name){
  
  subject <- strsplit(as.vector(unique(SeuratObj$orig.ident)), "[.]")[[1]][1]
  fontsize <- 4
  
  #helper functions
  give.n <- function(x){
    return(data.frame(y = median(x)+0.15, label = paste0("n = ", length(x))))
  }
  
  
  AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                    "granulocyte monocyte progenitor", "hematopoietic stem cell",
                    "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")
  
  df_meta <- SeuratObj@meta.data
  plots <- list()
  gene_exp <- as.vector(SeuratObj@assays[["SCT"]]@counts[gene_name, ])
  
  
  
    # wilcox.test(hsc_npm1, hscmag_npm1, alternative = "two.sided")
    
  df1 <- as.data.frame(gene_exp,
                         col.names = c("counts"), make.names = TRUE)
  df1$ident <- df_meta$predicted.celltype.l2
  df1$mut <- to_vec(for(i in df_meta$l2CellType.mut) if (str_detect(i, "malig")) "malignant" else "normal")
  df1$l2CellType.mut <- df_meta$l2CellType.mut
  df_cnt <- as.data.frame(table(df1$l2CellType.mut))
  
  #select the number of cells more than 3 
  qualify_cell <- df_cnt %>% filter(Freq >3)
  cellname <- as.vector(qualify_cell$Var1)
  
  # df1 <- df1 %>% filter(l2CellType.mut %in% cellname)
  # df1[[df1$l2CellType.mut %in% cellname, ]]

  names(df1) <- c("counts", "Cell", "malignancy")
    
  p1 <- ggboxplot(df1, x = "malignancy", y = "counts",
                    color = "malignancy", size  = 0.1,
                  add = "jitter", palette=  c("#00AFBB", "#E7B800"),
                    ylab = "RNA expression", xlab = "",
                    main = paste("Subject:", subject, "\n", gene_name, sep = ""), 
                  facet.by ="Cell", short.panel.labs = FALSE) +
    stat_compare_means(label = "p.format", size = fontsize,
                         family = "Times New Roman", 
                      label.x = 1.3, label.y = 4 )+
    stat_summary(fun.data = give.n, geom = "text")
    
  p1 + facet_wrap( ~ variable, scales="free")
  return (p1)
  
}

boxplot_stat_gene(AML_sub, "CD70")
boxplot_stat_gene(AML_sub, "CD33")
boxplot_stat_gene(AML_sub, "IL3RA")
factor(AML_sub$CellType)
celltype.markers.list <- list()
my_DimPlot(AML_sub)

AML_cat <- unique(AML_sub$l2CellType.mut)
df_AML <- as.data.frame(table(AML_sub$l2CellType.mut))
df_AML[df_AML$Var1=="common dendritic progenitor", ]$Freq

for (i in AML_cat){
  
  if (df_AML[df_AML$Var1==i, ]$Freq <=3){
    
    next
  
  }
  celltype.markers.list[[i]] <- FindMarkers(AML_sub, ident.1 = i,
                                            ident.2 = NULL,
                                            logfc.threshold = 0.25, 
                                            test.use = "wilcox", 
                                            only.pos = TRUE)
}

for (i in (1:length(celltype.markers.list))){
  data_a <- celltype.markers.list[[i]]
  if (dim(data_a)[1] <= 25){
    data2 <- celltype.markers.list[[i]][, c(1,2)]
  }else{
    data2 <-celltype.markers.list[[i]][0:25, c(1,2)] }
  
  data2$significance <- to_vec(for (i in as.vector(data2$p_val)) if(i<=0.01) "p-value <= 0.01" else "p-value > 0.01")
  data2$gene <- rownames(data2)
  DGE_pvalue_plot(data2, paste0(names(celltype.markers.list)[[i]], sep = ""),
                  "other celltypes")
}

aa <- FindMarkers(AML_sub, ident.2 = "lymphoid-primed multipotent progenitor-malignant",
                  ident.1  = "lymphoid-primed multipotent progenitor",
                  logfc.threshold = 0.25, 
                  test.use = "wilcox", 
                  only.pos = TRUE)

data2 <- aa[, c(1,2)]
data2$significance <- to_vec(for (i in as.vector(data2$p_val)) if(i<=0.01) "p-value <= 0.01" else "p-value > 0.01")
data2$gene <- rownames(data2)

DGE_pvalue_plot(data2, "LMPP", "LMPP-malignant")

dev.off()
barplot(table(AML_nature$mutated_gene))


