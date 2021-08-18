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
# library(biomaRt)
# library(SingleCellExperiment)
library(scater)
library(SingleR)
library(patchwork)
library(ggpubr)
library(comprehenr)
library(EnhancedVolcano)

project<-"AML-Nature"
Sys.setenv(language="en")
if (project == "AML-Nature"){
  working_dir<-"C:/Users/jtao/work_dir/AML/data/AMLnature"
  graph_dir<- "C:/Users/jtao/work_dir/AML/data/AMLnature/images/"
  setwd(working_dir)
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}



#===============================================================================
#helper functions
#===============================================================================
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
      width=png_width, height=png_height, units="in", res=300)
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

RECODE_Ontology_mix <- function(Seuratobj){
  
  Seuratobj$celltype.bmcite <- recode(Seuratobj$celltype.bmcite, "CD14 Mono"="CD14-positive, CD16-negative classical monocyte",
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

add.mut.CellType.mix <- function(SeuratObj, reference = TRUE){
  cell.ident <- c()
  
  if (reference == TRUE){
    for (i in 1:length(SeuratObj$celltype.bmcite)){
      
      CellType <- as.vector(SeuratObj$celltype.bmcite)
      mut <- SeuratObj$mutation_annotation
      
      if (mut[i] == "mutation detected" & is.na(mut[i])==FALSE ){
        cell.ident <- append(cell.ident, 
                            paste(CellType[i], "malignant",
                                  sep =  "-"))
      } else if(is.na(mut[i]) == TRUE){
        cell.ident <- append(cell.ident, paste(CellType[i], "bmcite",
                                            sep = "-"))
        
      }
      else{cell.ident <- append(cell.ident, CellType[i])}
      
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
  fontsize <- 7
  
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
  
  lab_y <-  max(gene_exp) * (max(gene_exp)/(max(gene_exp) - min(gene_exp)) - 0.04)
  
  
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
                  facet.by ="Cell", short.panel.labs = FALSE,
                  font.label = list(size = 14, face = "bold")) +
    stat_compare_means(label = "p.format", size = fontsize,
                       family = "mono", 
                       label.x = 1.3, label.y = lab_y)+
    stat_summary(fun.data = give.n, geom = "text", size = fontsize) + 
    theme(text = element_text(size = fontsize + 3),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"),
                                      face = "bold", size = fontsize ),
          axis.text.x  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize),
          aspect.ratio = 1/1) + 
    font("legend.title",  face = "bold", size = fontsize)+
    font("legend.text", size = fontsize +6)+
    font("xy.text", size = fontsize +4 )+
    font("title", size = fontsize +4 )
    
  
  p2 <- facet(p1, facet.by = "Cell", 
             panel.labs.font = list(size = fontsize +3))
    # + 
    # facet_wrap( ~ variable, scales="free")
  return (p1)
  
}

boxplot_stat_gene_bmcite <- function(SeuratObj, gene_name, subject){
  
  # subject <- strsplit(as.vector(unique(SeuratObj$orig.ident)), "[.]")[[1]][1]
  fontsize <- 7
  
  #helper functions
  give.n <- function(x){
    return(data.frame(y = 1.7, label = paste0("n = ", length(x))))
  }
  
  
  AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                    "granulocyte monocyte progenitor", "hematopoietic stem cell",
                    "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")
  
  df_meta <- SeuratObj@meta.data
  plots <- list()
  gene_exp <- as.vector(SeuratObj@assays[["integrated"]]@data[gene_name, ])
  
  lab_y <-  max(gene_exp) * (max(gene_exp)/(max(gene_exp) - min(gene_exp)) - 0.04)
  
  
  # wilcox.test(hsc_npm1, hscmag_npm1, alternative = "two.sided")
  
  df1 <- as.data.frame(gene_exp,
                       col.names = c("counts"), make.names = TRUE)
  df1$ident <- df_meta$celltype.bmcite
  df1$mut <- to_vec(for(i in df_meta$l2CellType.mut) if (str_detect(i, "malig")) "malignant" else if(str_detect(i, "bm")) "bmcite" else "normal")
  df1$l2CellType.mut <- df_meta$l2CellType.mut
  # df_cnt <- as.data.frame(table(df1$l2CellType.mut))
  
  #select the number of cells more than 3 
  # qualify_cell <- df_cnt %>% filter(Freq >3)
  # cellname <- as.vector(qualify_cell$Var1)
  
  # df1 <- df1 %>% filter(l2CellType.mut %in% cellname)
  # df1[[df1$l2CellType.mut %in% cellname, ]]
  
  names(df1) <- c("counts", "Cell", "malignancy", "original_name")
  # print(head(df1, 10, 10))
  var_width <- 30
  df1 <- df1 %>% mutate(cell = str_wrap(Cell, width = var_width))
  
  p1 <- ggviolin(df1, x = "malignancy", y = "counts",
                  color = "malignancy", size  = 0.5,
                  add = "none", palette=  c("#00AFBB", "#E7B800", "#2D2D2D"),
                  ylab = "RNA expression", xlab = "",
                  main = paste("Subject:", subject, "\n", gene_name, sep = ""), 
                  short.panel.labs = FALSE,
                  font.label = list(size = 7, face = "bold")) +
    theme(aspect.ratio = 1/1)
   
  
  
  p2<- facet(p1, facet.by = "cell", nrow =1, ncol=4, 
             panel.labs.font.x = list(size= 15, face = "bold"))+ 
    stat_compare_means(label = "p.format", size = 4,family = "mono", 
                       label.x = 1.3, label.y = lab_y,
                       ref.group = "bmcite", 
                       method.args = list(alternative = "less"))+
    stat_summary(fun.data = give.n, geom = "text", size = 4) + 
    theme(text = element_text(size = fontsize +2),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"),
                                      face = "bold", size = fontsize ),
          axis.text.x  = element_text(size = fontsize+2),
          legend.text = element_text(size = fontsize+2)) + 
    font("legend.title",  face = "bold", size = fontsize+2)+
    font("legend.text", size = fontsize+2)+
    font("xy.text", size = fontsize+2)+
    font("title", size = fontsize +2)
  return (p2)
}

#================================================================================
# analysis 
#================================================================================

#1. target visualization 

AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                  "granulocyte monocyte progenitor", "hematopoietic stem cell",
                  "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")


AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_508084.anno.rds")
AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_548327.anno.rds")
AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_721214.anno.rds")
AML_nature <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_782328.anno.rds")

# AML_nature <- Sobj
AML_nature <- RECODE_Ontology(AML_nature)
AML_nature <- add.mut.CellType(AML_nature, reference = TRUE)

Idents(AML_nature) <- AML_nature$predicted.celltype.l2

mutual_celltype <- intersect(AML_category, unique(AML_nature$predicted.celltype.l2))
AML_sub <- subset(AML_nature, ident = mutual_celltype)

table(AML_nature$predicted.celltype.l2)
# data<- c()
boxplot_stat_gene(AML_nature, "CD70")
boxplot_stat_gene(AML_sub, "CD33")
boxplot_stat_gene(AML_sub, "IL3RA")


pil <- boxplot_stat_gene(AML_sub, "IL3RA")
save_png(pil, "80_il3ra", 15, 10)
p70 <- boxplot_stat_gene(AML_sub, "CD70")
save_png(p70, "80_cd70", 15, 10 )
p33 <- boxplot_stat_gene(AML_sub, "CD70")
save_png(p33, "80_cd33", 15, 10 )

AML_nature$predicted.celltype.l1
cnt_RNA <- AML_nature@meta.data[, c("predicted.celltype.l2", "nCount_SCT")]
feature_n <- AML_nature@meta.data[, c("predicted.celltype.l2", "nFeature_SCT")]
p <- ggplot(cnt_RNA, aes(x=predicted.celltype.l2, y=nCount_SCT)) + 
  geom_violin()+ theme(axis.text.x = element_text(angle = 90)) + 
  geom_boxplot(width=0.1)

p

p2 <- ggplot(feature_n, aes(x=predicted.celltype.l2, y=nFeature_SCT)) + 
  geom_violin()+ theme(axis.text.x = element_text(angle = 90)) + 
  geom_boxplot(width=0.1)


#===============================================================================
#data preprocessing
#===============================================================================
#1. bmcite and AML data

bmcites <- c("hematopoietic stem cell-bmcite", "megakaryocyte progenitor cell-bmcite",
             "lymphoid-primed multipotent progenitor-bmcite", "granulocyte monocyte progenitor-bmcite"
             )
malignancy <- c("hematopoietic stem cell-malignant", "megakaryocyte progenitor cell-malignant",
                "lymphoid-primed multipotent progenitor-malignant", 
                "granulocyte monocyte progenitor-malignant"
                )

Sobj <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/721214.bmcite.rds")

bmcite <- Sobj@meta.data["celltype.l2"] %>% drop_na()
subject <- Sobj@meta.data["predicted.celltype.l2"] %>% drop_na()
colnames(bmcite) <- "celltype.bmcite"
colnames(subject) <- "celltype.bmcite"

celltype.bmcite <- rbind(bmcite, subject)

if (all.equal(rownames(Sobj@meta.data), rownames(celltype.bmcite))){
  Sobj$celltype.bmcite <- celltype.bmcite$celltype.bmcite
}


Sobj <- RECODE_Ontology_mix(Sobj)
Sobj <- add.mut.CellType.mix(Sobj, reference = TRUE)


table(Sobj$l2CellType.mut)
Idents(Sobj) <- Sobj$celltype.bmcite
# unique(Sobj$celltype.bmcite)
Idents(Sobj) <- Sobj$l2CellType.mut
mutual_celltype <- intersect(AML_category, unique(Sobj$celltype.bmcite))
AML_sub <- subset(Sobj, ident = c(bmcites, malignancy))



cd33 <- boxplot_stat_gene_bmcite(AML_sub, "CD33", "721214")
cd70 <- boxplot_stat_gene_bmcite(AML_sub, "CD70", "721214")
IL3RA <- boxplot_stat_gene_bmcite(AML_sub, "IL3RA", "721214")
save_png(cd33, "CD33", 8,5)
save_png(cd70, "CD70", 8,5)
save_png(IL3RA, "IL3RA", 8,5)
AML_sub$celltype.bmcite


unique(AML_sub$orig.ident)

#===============================================================================
#differential gene analysis
#===============================================================================
AML_cat <- unique(AML_sub$l2CellType.mut)
df_AML <- as.data.frame(table(AML_sub$l2CellType.mut))
AML_cat
celltype.markers.list <- list()
Idents(AML_sub)<- AML_sub$l2CellType.mut

bmcites <- c("hematopoietic stem cell-bmcite", "megakaryocyte progenitor cell-bmcite",
             "lymphoid-primed multipotent progenitor-bmcite", "granulocyte monocyte progenitor-bmcite",
             "erythroid progenitor cell-bmcite")
malignancy <- c("hematopoietic stem cell-malignant", "megakaryocyte progenitor cell-malignant",
                "lymphoid-primed multipotent progenitor-malignant", 
                "granulocyte monocyte progenitor-malignant",
                "erythroid progenitor cell-malignant")
df_AML
# vocalno_plot_list <- list()
for (i in (1:length(malignancy))){
  mag <- malignancy[i]
  bmci <- bmcites[i]
  
  if (df_AML[df_AML$Var1==mag, ]$Freq <=3|df_AML[df_AML$Var==bmci, ]$Freq <=3){
    
    next
    
  }
  celltype.markers.list[[mag]] <- FindMarkers(AML_sub, ident.1 = mag,
                                            ident.2 = bmci,
                                            logfc.threshold = 0.25, 
                                            test.use = "wilcox", 
                                            only.pos = FALSE)
  
  
  
}

# data<- c()

#plot volcano plot

i = 4
# DefaultAssay(AML_sub)
# hsc <- FindMarkers(AML_sub, ident.1 = "hematopoietic stem cell",
#             ident.2 = "hematopoietic stem cell-malignant",
#             logfc.threshold = 0.25, 
#             test.use = "wilcox", 
#             only.pos = FALSE)

P.vol <- EnhancedVolcano(celltype.markers.list[[i]],
                          lab = rownames(celltype.markers.list[[i]]),
                          x = "avg_log2FC",
                          y = "p_val", FCcutoff = 0.5,
                         title = paste0(names(celltype.markers.list[i]), " vs \n",
                                                      bmcites[i]))

P.vol
table(AML_sub$l2CellType.mut)



#===============================================================================
#co-expression analysis 
#===============================================================================


plot_coexpression <- function(seurat_object_sample, celltype) {
  co_markers <- c("CD33", "IL3RA")
  if (length(unique(seurat_object_sample$orig.ident))==1){
    sample <- unique(seurat_object_sample$orig.ident)
  }else if (str_detect(celltype, "bmcite")){
   sample<- "bmcite"
    
  }else if (str_detect(celltype, "malig")){
    sample<- unique(seurat_object_sample$orig.ident)[[2]]
  }
  title <- paste0("sample:",
                  sample, "\n", celltype)
  # title <- ("CTCL sample")
  df_meta <- seurat_object_sample@meta.data
  cell_UMI <- rownames(df_meta %>% filter(l2CellType.mut == celltype))
  df<- seurat_object_sample@assays[["SCT"]]@data[co_markers, cell_UMI ]
  
  df_aa <- data.frame("CD33" = as.vector(df["CD33", ]),
                      "IL3RA" = as.vector(df["IL3RA", ]))
  # df_aa <- df_aa[df_aa[, "CD33"] >0 & df_aa[, "IL3RA"] >0, ]
  df_aa$color <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "orchid",
                        ifelse(df_aa$CD33>0.5, "skyblue3"
                               , ifelse(df_aa$IL3RA>0.5, "indianred3", "grey59")) )
  df_aa$label <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "CD33 > 0.5, IL3RA > 0.5",
                        ifelse(df_aa$CD33>0.5 & df_aa$IL3RA<0.5, "CD33 > 0.5, IL3RA < 0.5"
                               , ifelse(df_aa$IL3RA>0.5, "CD33 < 0.5, IL3RA > 0.5", "CD33 < 0.5, IL3RA < 0.5")) )
  
  sp <- ggplot(df_aa) +
    geom_point(aes(x = IL3RA, y = CD33, colour = label), size  = 2) + 
    scale_color_manual(values =c("grey59", "indianred3","skyblue3","orchid")) + 
    geom_hline(yintercept= 0.5, linetype = "dotted", color = "red", 
               size = 1)  +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size = 1) +
    labs(x = 'IL3RA Expression',y = "CD33 Expression", colour = "Expression Group", title = "") + 
    theme(axis.title= element_text(size = 13),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 11),
          aspect.ratio = 1/1) +coord_fixed()
  
  sp
  bp <- ggplot(df_aa, aes(x=label, fill = color)) + geom_bar(width = 0.5)+ coord_flip() +
    scale_fill_manual(values = c("grey59", "indianred3", "orchid","skyblue3")) +
    theme (legend.position = "none", aspect.ratio = 1.5/1, axis.text = element_text(size = 10),
           axis.title.y = element_blank()
    ) + labs(y = "Cell Counts")
  
  plot_row <- plot_grid(sp,bp, align = "h", rel_widths = c(2,1.2)) 
  ptitle <- ggdraw() + draw_label(title, fontface='bold', size = 9)
  plot_row <- plot_grid(ptitle, plot_row , ncol=1, rel_heights=c(0.05, 1)) # rel_heights values control title margins
  
  return (plot_row)
  
}
# malignancy <- c("hematopoietic stem cell-malignant", "megakaryocyte progenitor cell-malignant",
#                 "lymphoid-primed multipotent progenitor-malignant", 
#                 "granulocyte monocyte progenitor-malignant",
#                 "erythroid progenitor cell-malignant")

coexp_plots <- list()
for (i in malignancy){
  coexp_plots[[i]] <- plot_coexpression(AML_nature,i) 
}


coexp_plots[[5]]
plot_coexpression(AML_nature, "common dendritic progenitor-malignant")

DoHeatmap(AML_nature, features = c("CD33", "IL3RA"))


table(AML_nature$l2CellType.mut)
Idents(AML_nature) <- AML_nature$l2CellType.mut
Idents(AML_nature) <- AML_nature$predicted.celltype.l2
aa <- my_DimPlot(AML_nature)
aa
save_png(aa, "dimplot_78_l2", 12, 8)
DefaultAssay(AML_nature) <- "RNA"
DoHeatmap(AML_nature, features =, slot = "data")




#===============================================================================
# T cell expression 
#===============================================================================
bm <- readRDS("C:/Users/jtao/work_dir/AML/data/SCT_bm_8000.rds")
hcabm40k <- readRDS("C:/Users/jtao/work_dir/AML/data/SCT_hcabm40k.rds")
S_Tcell <- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/Tcell_SCT_72_78.rds")
S_Tcell <- RECODE_Ontology(S_Tcell)

bm <- RECODE_Ontology(bm)
Tmarkers <- c("CD3D", "CD3E", "CD8A","CD8B", "CD4", "PD-1",
                 "CD279", "PDCD1", "TRGC1", "TRDC1","TRGC2")
T_reg_marker <- c("CCR8", "CTLA4", "IKZF2",	"IKZF4", "FOXP3", "IL10", "TNFRSF18")
cytotoxi_marker <- c("CCL3", "IFNG",	"GZMB", "CCL4",	"GZMA",	"PRF1",	"CST7", "NKG7")
ckpoint_marker <- c("CTLA4",	"LAG3",	"HAVCR2",	"PDCD1", "TIGIT")
naive <- c("CCR7",	"LEF1", "SELL",	"TCF7")
GE8 <- c("IFNG",	"CXCL9", "CD8A",	"GZMA",	"GZMB",	"CXCL10",	"PRF1",	"TBX21")
markers<- c(Tmarkers, T_reg_marker,cytotoxi_marker,ckpoint_marker, naive,GE8
                 )

DefaultAssay(S_Tcell) <- "SCT"
Idents(S_Tcell) <- S_Tcell$predicted.celltype.l2
Tcell<- subset(S_Tcell,
               ident = grep("T", unique(S_Tcell$predicted.celltype.l2), value = TRUE))
Idents(Tcell) <- Tcell$predicted.celltype.l2
DefaultAssay(Tcell) <- "SCT"
heatmap <- DoHeatmap(Tcell, 
                     # title = unique(Tcell$orig.ident),
                     features = markers,   
                     size = 4, slot = "data", combine = TRUE,
                     label = TRUE)+ 
  theme(legend.position= "bottom", 
        legend.text = element_text(size = 11))+ 
  scale_fill_distiller(palette = "RdBu") 

save_png(heatmap, "heatmap_tcell_2AML_markers", 18, 12) 
heatmap
"TRGC1" %in% rownames(Tcell@assays[["RNA"]])

trace(DoHeatmap, edit = TRUE)
Idents(AML_nature) <- AML_nature$predicted.celltype.l2
DefaultAssay(S_Tcell) <- "integrated"
VariableFeatures(S_Tcell)

S_Tcell<- RunUMAP(S_Tcell, slot = "data", features = VariableFeatures(S_Tcell))
S_Tcell <-RunPCA(S_Tcell, slot = "data", features = VariableFeatures(S_Tcell))
S_Tcell <-FindNeighbors(S_Tcell, features = VariableFeatures(S_Tcell))
S_Tcell <- FindClusters(S_Tcell, algorithm = 3, resolution= 0.5)
Idents(S_Tcell) <- S_Tcell$seurat_clusters
my_DimPlot(S_Tcell)


AvailableData()





