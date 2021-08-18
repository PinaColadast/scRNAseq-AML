library(remoter)
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
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(patchwork)
library(ggpubr)
library(comprehenr)
library(EnhancedVolcano)
library(CiteFuse)
library(comprehenr)


project<-"protein-genomics"
Sys.setenv(language="en")
if (project == "protein-genomics"){
  working_dir<-"C:/Users/jtao/work_dir/AML/data/pro-genomics"
  graph_dir<- "C:/Users/jtao/work_dir/AML/data/pro-genomics/images/"
  setwd(working_dir)
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}

# load("healthy.rda")
# data <- WTA
# data <- UpdateSeuratObject(object=data)

BM_All <- readRDS("healthy.rds")
BM_All <- UpdateSeuratObject(object = BM_All)
#===============================================================================
#helper function 
#===============================================================================
plot_coexpression_ab <- function(seurat_object_sample, celltype) {
  co_markers <- c("CD33-AB", "CD123-AB")
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
  cell_UMI <- rownames(df_meta %>% filter(Prediction_Ind == celltype))
  df<- seurat_object_sample@assays[["AB"]]@data[co_markers, cell_UMI]
  
  df_aa <- data.frame("CD33-AB" = as.vector(df["CD33-AB", ]),
                      "CD123-AB" = as.vector(df["CD123-AB", ]))
  # df_aa <- df_aa[df_aa[, "CD33"] >0 & df_aa[, "IL3RA"] >0, ]
  df_aa$color <- ifelse(df_aa$CD33.AB>0.5 & df_aa$CD123.AB >0.5, "orchid",
                        ifelse(df_aa$CD33.AB>0.5, "skyblue3"
                               , ifelse(df_aa$CD123.AB>0.5, "indianred3", "grey59")) )
  df_aa$label <- ifelse(df_aa$CD33.AB>0.5 & df_aa$CD123.AB >0.5, "CD33 > 0.5, IL3RA > 0.5",
                        ifelse(df_aa$CD33.AB>0.5 & df_aa$CD123.AB<0.5, "CD33 > 0.5, IL3RA < 0.5"
                               , ifelse(df_aa$CD123.AB>0.5, "CD33 < 0.5, IL3RA > 0.5", "CD33 < 0.5, IL3RA < 0.5")) )
  
  sp <- ggplot(df_aa) +
    geom_point(aes(x = CD123.AB, y = CD33.AB, colour = label), size  = 2) + 
    scale_color_manual(values =c("grey59", "indianred3","skyblue3","orchid")) + 
    geom_hline(yintercept= 0.5, linetype = "dotted", color = "red", 
               size = 1)  +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size = 1) +
    labs(x = 'IL3RA Expression',y = "CD33 Expression", colour = "Expression Group", title = "") + 
    theme(axis.title= element_text(size = 13),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 11),
          aspect.ratio = 1/3) +coord_fixed()
  
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


boxplot_stat_gene_mix <- function(SeuratObj, gene_name){
  
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
  gene_exp <- as.vector(SeuratObj@assays[["RNA"]]@counts[gene_name, ])
  
  lab_y <-  max(gene_exp) * (max(gene_exp)/(max(gene_exp) - min(gene_exp)) - 0.04)
  
  
  # wilcox.test(hsc_npm1, hscmag_npm1, alternative = "two.sided")
  
  df1 <- as.data.frame(gene_exp,
                       col.names = c("counts"), make.names = TRUE)
  df1$ident <- df_meta$Cell_Type4
  df1$mut <- to_vec(for(i in df_meta$Batch) if (str_detect(i, "AML")) "malignant" else if(str_detect(i, "BM")) "normal" else "normal")
  
 
  
  names(df1) <- c("counts", "Cell", "malignancy")
  # print(head(df1, 10, 10))
  var_width <- 30
  df1 <- df1 %>% mutate(cell = str_wrap(Cell, width = var_width))
  
  p1 <- ggviolin(df1, x = "malignancy", y = "counts",
                  color = "malignancy", size  = 0.5,
                  add = "none", palette=  c("#E7B800", "#00AFBB", "#2D2D2D"),
                  ylab = "RNA expression", xlab = "",
                  main = paste("\n", gene_name, sep = ""), 
                  short.panel.labs = FALSE,
                  font.label = list(size = 7, face = "bold")) +
    theme(aspect.ratio = 1/1)
  
  
  
  p2<- facet(p1, facet.by = "cell", nrow =1, ncol=4, 
             panel.labs.font.x = list(size= 15, face = "bold"))+ 
    stat_compare_means(label = "p.format", size = 4,family = "mono", 
                       label.x = 1.3, label.y = lab_y )+
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


RECODE_Ontology <- function(Seuratobj){
  Seuratobj$cell_type <- recode(Seuratobj$Prediction_Ind, 'B cell progenitor'='B cell progenitor',
  'CD11c+ memory B cells'='Memory B cell',
  'CD141+ cDCs'='cDCs',
  'CD1c cDCs'='cDCs',
  'CXCR5-neg memory B cells?'="memory B cell",
  'Central memory / effector memory CD4 T cell'='CD4 T cell',
  'Central memory CD8 T cell'="CD8 T cell",
  'Central memory CD8 T cells'="CD8 T cell",
  'Class switched B cell'="B cell",
  'Class-switched memory B cells'="memory B cell",
  'Conventional dendritic cells (cDCs)'="cDCs",
  'Cycling MPPs'='MPP',
  'Cytotoxic CD4 T cells'='CD4 T cell',
  'Doublet T and Neutrophils'='Doublet T and Neutrophils',
  'Doublet and Triplets'='Doublet and Triplets',
  'Early EMP'="EMP",
  'Early GMP'="GMP",
  'Early MPPs'="MPP",
  'Early erythroid progenitor'='early erythroid progenitor',
  'Effector / memory CD4 T cells'='memory CD4 T cell',
  'Effector / memory CD8 T cells'='memory CD8 T cell',
  'Effector CD8 T cells'='CD8 T cell',
  'EoBasoMast'="EoBasoMast",
  'Eosinophil / Basophil progenitors'='eosinophil/baslphil progenitors',
  'Erythroid / Megakaryocyte progenitors'='MEP',
  'Erythroid progenitors'='erythroid progenitor',
  'GMP'='GMP',
  'Gamma delta T cells'='gamma delta T cell',
  'HSC'='HSC',
  'Hematopoietic stem cells (HSCs)'='HSC',
  'Immature B cells'='B cell',
  'LMPP'='LMPP',
  'LMPPs'="LMPP",
  'Late erythroid progenitor'='erythroid progenitor',
  'MEP'='MEP',
  'MPP'='MPP',
  'Mature naive B cells'='B cell',
  'Megakaryocyte progenitor'='megakaryocyte progenitor',
  'Megakaryocyte progenitors'="megakaryocyte progenitor",
  'Metaphase MPP' ="MPP",
  'Metaphase MPPs' ="MPP",
  'Monocytes'="monocyte",
  'Multipotent progenitors (MPPs)'='MPP',
  'Myelocyte'='myelocyte',
  'NA'='NA',
  'NK cell progenitor'='NK cell progenitor',
  'NK cells'='NK cell',
  'Naive CD4 T cells'='CD4 T cell',
  'Naive CD8 T cells'='CD8 T cell',
  'Neutrophil'='Neutrophil',
  'Neutrophils'='Neutrophil',
  'Plasma Cells'='plasma cell',
  'Plasma cells'='plasma cell',
  'Plasmacytoid DC progenitors'='plasmacytoid DC progenitors',
  'Plasmacytoid DCs'="Plasmacytoid DCs",
  'Plasmacytoid dendritic cells (pDCs)'='plasmacytoid DCs',
  'Pre-B cells'='pre-B cell',
  'Pro-B cell'='pro-B cell',
  'Promyelocytes'='promeylocyte',
  'Promyelocytes / GMPs' ="promyelocyte",
  'T cell progenitors'='T cell progenitor',
  'Tissue resident memory CD8 T cells'='memory CD8 T cell',
  'Transitional (T1) B cells'='B cell',
  'Transitional (T2) B cell'='B cell',
  'Transitional (T2) B cells'='B cell',
  'Unswitched memory B cells'='B cell',
  'bottom: cytotoxic CD4 T cells'='CD4 T cell',
  'doublets from Erys, T cell'='doublets from Erys, T cell',
  'early MEP'="MEP",
  'pDC progenitors'='plasmacytoid DC progenitors',
  'pro-B cells'='pro-B cell')
  return(Seuratobj)}

#renameing 


#===============================================================================
#data processing 
#===============================================================================
s97 <- NormalizeData(BM_All, normalization.method = "CLR", margin = 2, assay = "AB")
BM_All<-RECODE_Ontology(BM_All)
#===============================================================================
#coexpression analysis 
#===============================================================================
DefaultAssay(s97) <- "AB"

df <- as.data.frame(table(s97
                          @meta.data[["Prediction_Ind"]]))
cell.list<- c("GMP", "HSC", "Early erythroid progenitor", "Cycling MPPs", "early MEP",
             "Late erythroid progenitor", "LMPP", "Megakaryocyte progenitor", 
             "MEP", "MPP", "Plasmacytoid DC progenitors")
cell.list
  
coexp_plots <- list()
for (i in cell.list){
  coexp_plots[[i]] <- plot_coexpression_ab(s97,i) 
}

coexp_plots[2]
Idents(s200) <- s200$Prediction_Ind
markers <- c("CD3-AB", "CD20-AB")
heatmap <- DoHeatmap(s200, 
                     features = markers,   
                     size = 3.5, slot = "data", label = TRUE)+
  labs(fills = "Cluster", colours = "Expression") + 
  theme(legend.position= "bottom", 
        legend.text = element_text(size = 9))+ 
  scale_fill_distiller(palette = "RdBu") 
heatmap
# trace(DoHeatmap, edit=TRUE)

#===============================================================================
#Citefuse 
#===============================================================================
#pbmc
pbmc <- readRDS("C:/Users/jtao/work_dir/CTCL/data/pbmc10k-multimodal.rds")
DefaultAssay(pbmc) <- "ADT"
Idents(pbmc)<- pbmc$cell.type

sobj_tcell <- subset(pbmc, ident = c("CD4 T cell, naive", "CD4 T cell, memory"))
sobj_bcell <- subset(pbmc, ident = "B cell")

SceObj_pbmc <- as.SingleCellExperiment(pbmc)
SceObj_tcell <-  as.SingleCellExperiment(sobj_tcell)
SceObj_bcell <-  as.SingleCellExperiment(sobj_bcell)

p_pbmc <- visualiseExprs(SceObj_tcell,
                     plot = "pairwise", exprs_value = "logcounts",
                     feature_subset = c("CD4", "CD3"), threshold= c(0.2, 0.4))
p_pbmc <- visualiseExprs(SceObj_tcell,
                         plot = "pairwise", exprs_value = "logcounts",
                         feature_subset = c("CD4", "CD3"))

cell.list<- c("GMP", "HSC", "MPP", "MEP",
              "LMPP", "Plasmacytoid DC progenitors")

Idents(BM_All) <- BM_All$cell_type
sobj_subset <- list()
for (i in cell.list){
  obj <- subset(BM_All, ident = i)
  sobj_subset[[i]] <- obj
}

i = 3
sobj_subset[[i]]
SceObj <-  as.SingleCellExperiment(sobj_subset[[i]])
p4 <- visualiseExprs(SceObj,
                    plot = "pairwise", exprs_value = "logcounts",
                    feature_subset = c("CD33-AB", "CD123-AB"),
                    threshold=c(0.5,1))

# trace(visualiseExprs, edit = TRUE)
sce <- as.SingleCellExperiment(BM_All)
# trace(as.SingleCellExperiment, edit = TRUE)

?visualiseExprs
#===============================================================================
#co-expression on RNA level
#===============================================================================
unique(BM_All$orig.ident)
unique(BM_All$Cell_Type4)
table(BM_All$Batch)

cell.types <- c("HSCs & MPPs", "Lymphoid-primed multipotent progenitors",
                "Megakaryocyte progenitors", "Early erythroid progenitor")
Idents(BM_All) <- BM_All$Batch
AML_BM <- subset(BM_All, ident = c("AML1", "AML3", "BM1", "BM2", "BM3"))
Idents(AML_BM) <- AML_BM$Cell_Type4
AML_BM <- subset(AML_BM, ident = cell.types)

df_meta <- AML_BM@meta.data
AML_BM$malignancy <- to_vec(for(i in df_meta$Batch) if (str_detect(i, "AML")) "malignant" else if(str_detect(i, "BM")) "normal")

#differential analysis of respective marker

"IL3RA" %in% rownames(AML_BM@assays[["RNA"]]@counts)
CD33 <- boxplot_stat_gene_mix(AML_BM, "CD33")
IL3RA <- boxplot_stat_gene_mix(AML_BM, "IL3RA")

save_png(CD33, "prot-geno-CD33", 8,5)
save_png(IL3RA, "prot-geno-IL3RA", 8,5)

#co-expression by scattering plots 

cell.order <-  c("Early erythroid progenitor",
                 "HSCs & MPPs", "Lymphoid-primed multipotent progenitors",
                              "Megakaryocyte progenitors" )
cell.order.malignant <- to_vec(for (i in cell.order) paste(i, "-malignant", sep = ""))
unique(AML_BM$orig.ident)
Idents(AML_BM) <- AML_BM$malignancy
AML_BM_mag <- subset(AML_BM, ident= "malignant")
AML_BM_bm <- subset(AML_BM, ident = "normal") 

df_meta <- AML_BM_bm@meta.data

co_expression_matrix <- as.data.frame(t(AML_BM_bm@assays[["RNA"]]@counts[co_markers,]))
rownames(co_expression_matrix) <-c()
co_expression_matrix$Cell <- df_meta$Cell_Type4
co_expression_matrix$mut <- df_meta$malignancy
co_expression_matrix <- mutate(co_expression_matrix, Cell = paste(Cell,"-", mut, sep = ""))

df_aa<-co_expression_matrix
df_aa$color <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "orchid",
                      ifelse(df_aa$CD33>0.5, "skyblue3"
                             , ifelse(df_aa$IL3RA>0.5, "indianred3", "grey59")) )
df_aa$label <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "CD33 > 0.5, IL3RA > 0.5",
                      ifelse(df_aa$CD33>0.5 & df_aa$IL3RA<0.5, "CD33 > 0.5, IL3RA < 0.5"
                             , ifelse(df_aa$IL3RA>0.5, "CD33 < 0.5, IL3RA > 0.5", "CD33 < 0.5, IL3RA < 0.5")) )
var_width_bm <- 20
df_aa <- mutate(df_aa, Cell.Name = str_wrap(Cell, width = var_width_bm))

proportion_bm <- c()
n_cells_bm <-c()
for (i in cell.order){
  sample_cell <- df_aa %>% filter(Cell==i)
  if (nrow(sample_cell %>% filter(label == "CD33 > 0.5, IL3RA > 0.5")) < 1){
    freq <- 0
  }else{
    freq <- prop.table(table(sample_cell$label))[["CD33 > 0.5, IL3RA > 0.5"]]
  }
  n_cells_bm <- append(n_cells_bm, nrow(sample_cell))
  proportion_bm <- append(proportion_bm, freq)
}


dat_sum_bm <- data.frame(
  num_of_cells = n_cells_bm,
  proportion = proportion_bm,
  Cell.Name = cell.order
)

dat_sum_bm <- mutate(dat_sum_bm, Cell.Name = str_wrap(Cell.Name, width = var_width_bm))

dat_sum_bm$text <- sprintf("number of cells = %s\nCD33++IL3RA++ %s", format(dat_sum_bm$num_of_cells, nsmall = 0), 
                           percent(dat_sum_bm$proportion))



gbmp <- ggplot(df_aa)  + 
  geom_point(aes(x = IL3RA, y = CD33, colour = label), size  = 2) + 
  scale_color_manual(values =c("grey59", "indianred3","skyblue3","orchid")) + 
  geom_hline(yintercept= 0.5, linetype = "dotted", color = "red", 
             size = 1)  +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size = 1) +
  labs(x = 'IL3RA Expression',y = "CD33 Expression", colour = "Expression Group", title = "") + 
  theme(axis.title= element_text(size = 8+4),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8+4),
        aspect.ratio = 1/1.2) +coord_fixed()+
  facet_grid(cols = vars(Cell.Name), margins = FALSE)+
  facet_wrap(~ Cell.Name, ncol=2)+
  geom_text(check_overlap = TRUE,
            size    = 4.5,
            data    = dat_sum_bm,
            mapping = aes(x = Inf, y = Inf, label = text)
            ,hjust   = 1.05,
            vjust   = 1.5
  )+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 12, face = "bold"),
  )
gbmp
save_png(gbmp, "targeted-panel-AMLcoexpress", 10, 10)

cd123_ab <- as.vector(BM_All@assays[["AB"]]@counts["CD123-AB", ])
cd123_RNA <- as.vector(BM_All@assays[["RNA"]]@counts["IL3RA", ])

cd123_ab
cd33_ab <- as.vector(BM_All@assays[["AB"]]@counts["CD33-AB", ])
cd33_RNA <-as.vector(BM_All@assays[["RNA"]]@counts["CD33", ])

cd4_ab <- as.vector(BM_All@assays[["AB"]]@counts["CD4-AB", ])
cd4_RNA <-as.vector(BM_All@assays[["RNA"]]@counts["CD4", ])

ab_rna <- data.frame(
  cd33_ab = cd33_ab,
  cd33_RNA = cd33_RNA,
  cd123_ab = cd123_ab,
  cd123_RNA = cd123_RNA,
  cd4_ab = cd4_ab, 
  cd4_RNA = cd4_RNA
  
)

plot(cd33_ab, cd33_RNA)
ggscatter(ab_rna, x= "cd4_ab", y= "cd4_RNA",
          conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4 antibody expression", ylab = "CD4 RNA expression") 

rownames(BM_All@assays[["AB"]]@counts)


sceasy::convertFormat(bm, from="seurat", to="anndata", assay = "RNA", main_layer = "counts",
                      outFile='pbmc97ab_RNA.h5ad', drop_single_values = FALSE)

sceasy::convertFormat(bm, from="seurat", to="anndata", assay = "ADT", main_layer = "counts",
                      outFile='pbmc97ab_AB.h5ad', drop_single_values = FALSE)

sceasy::convertFormat(bm, from="seurat", to="anndata", assay = "BOTH", main_layer = "counts",
                      outFile='pbmc97ab.h5ad', drop_single_values = FALSE)



#===============================================================================
#from Scanpy to Seurat
#===============================================================================
library(reticulate)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert(paste(working_dir, "TotalVI_denoised_200ab.h5ad", sep = "/"), dest= "h5seurat")
bm <- LoadH5Seurat(paste(working_dir, "TotalVI_denoised_200ab.h5ad", sep = "/"))

