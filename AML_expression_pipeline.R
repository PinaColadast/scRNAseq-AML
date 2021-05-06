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

addline_format <- function(x,...){
  gsub('\\s','\n',x)
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

boxplot_stat_gene <- function(SeuratObj, gene_name){
  
  subject <- strsplit(as.vector(unique(SeuratObj$orig.ident)), "[.]")[[1]][1]
  fontsize <- 4
  
  #helper functions
  give.n <- function(x){
    return(data.frame(y = max(x)-0.5, label = paste0("n = ", length(x))))
  }
  
  
  AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                    "granulocyte monocyte progenitor", "hematopoietic stem cell",
                    "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")
  
  df_meta <- SeuratObj@meta.data
  plots <- list()
  gene_exp <- as.vector(SeuratObj@assays[["SCT"]]@counts[gene_name, ])
  
  lab_y <-  max(gene_exp) * (max(gene_exp)/(max(gene_exp) - min(gene_exp)) - 0.05)
  
  
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
                       family = "mono", 
                       label.x = 1.3, label.y = lab_y )+
    stat_summary(fun.data = give.n, geom = "text")
  
  p1 + facet_wrap( ~ variable, scales="free")
  return (p1)
  
}

#================================================================================
# analysis 
#================================================================================
AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                  "granulocyte monocyte progenitor", "hematopoietic stem cell",
                  "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")


# AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_508084.anno.rds")
AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_548327.anno.rds")
AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_721214.anno.rds")


AML_nature <- RECODE_Ontology(AML_nature)
AML_nature <- add.mut.CellType(AML_nature, reference = TRUE)

Idents(AML_nature) <- AML_nature$predicted.celltype.l2

mutual_celltype <- intersect(AML_category, unique(AML_nature$predicted.celltype.l2))
AML_sub <- subset(AML_nature, ident = mutual_celltype)

table(AML_nature$predicted.celltype.l2)

boxplot_stat_gene(AML_sub, "CD70")
boxplot_stat_gene(AML_sub, "CD33")
boxplot_stat_gene(AML_sub, "IL3RA")
