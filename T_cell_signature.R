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
library(RColorBrewer)


Sys.setenv(language="en")


project<-"AML-nature"
Sys.setenv(language="en")
if (project == "AML-nature"){
  working_dir<-"C:/Users/jtao/work_dir/AML/data/AMLnature"
  graph_dir<- "C:/Users/jtao/work_dir/AML/data/images/"
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
      width=png_width, height=png_height, units="in", res=200)
  print(plot)
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
                                            "Naive B"="naive B cell","NULL"="NULL",
                                            "CD4 Naive" = "CD4 naive T cell", "CD8 Naive" = "CD8 Naive T cell")
  
  return(Seuratobj)
}


RECODE_Ontology_bm <- function(Seuratobj){
  
  Seuratobj$celltype.l2 <- recode(Seuratobj$celltype.l2, "CD14 Mono"="CD14-positive, CD16-negative classical monocyte",
                                  "CD16 Mono"="CD14-low, CD16-positive monocyte", 
                                  "CD4 Memory"="effector memory CD4 T cell", "CD56 bright NK"="natural killer cell",
                                  "CD8 Effector_1"="effector CD8 T cell", "CD8 Effector_2"="effector CD8 T cell", "CD8 Memory_2"="effector memory CD8 T cell", 
                                  "CD8 Memory_1"="effector memory CD8 T cell", "cDC2"="CD1c-positive myeloid dendritic cell",
                                  "GMP" ="granulocyte monocyte progenitor", "HSC"="hematopoietic stem cell","LMPP"="lymphoid-primed multipotent progenitor" , 
                                  "Memory B"="memory B cell", "NK"="natural killer cell", "pDC"="plasmacytoid dendritic cell", "Plasmablast"="plasmablast", 
                                  "Prog_B 1"="pro-B cell", "Prog_B 2"="pro-B cell","Prog_DC"="common dendritic progenitor", "Prog_Mk"="megakaryocyte progenitor cell", "Prog_RBC"="erythroid progenitor cell",
                                  "Naive B"="naive B cell","NULL"="NULL",
                                  "CD4 Naive" = "CD4 naive T cell", "CD8 Naive" = "CD8 Naive T cell")
  return(Seuratobj)
  }
#===============================================================================
#analysis 
#===============================================================================

TH1_differ <- c("CD28","CD40","CD40LG","CD86", "HLA-DRA", "HLA-DRB1")
IFNG_production <- c("IDO1", "CXCL9", "CXCL10", "STAT1", "HLA-DRA", "IFNG")

cytotoxicity <- c("CCL3", "IFNG", "GZMB", "CCL4", "GZMA", "PRF1", "CST7",
                  "NKG7")

Checkpoint <- c("CTLA4", "LAG3", "HAVCR2", "PDCD1", "TIGIT")
Naive_T    <- c("CCR7", "LEF1", "SELL", "TCF7")



Sobj_Tcell <- readRDS(paste0(getwd(), "/Tcell_SCT_72_78.rds"))
Sobj_Tcell <- RECODE_Ontology(Sobj_Tcell)
Sobj_Tcell <- RunUMAP(Sobj_Tcell, features = VariableFeatures(object = Sobj_Tcell))


bmcite<- readRDS("C:/Users/jtao/work_dir/AML/data/SCT_bm_8000.rds")
bmcite <- RECODE_Ontology_bm(bmcite)

saveRDS(tcell, "C:/Users/jtao/work_dir/AML/data/SCT_BM_Tcell_8000.rds")

Idents(Sobj_Tcell) <- Sobj_Tcell$predicted.celltype.l2
Idents(bmcite) <- bmcite$celltype.l2
tcell <- subset(Sobj_Tcell, ident = grep("T", unique(Sobj_Tcell$predicted.celltype.l2), value = TRUE))
tcell <- subset(bmcite, ident = grep("T", unique(bmcite$celltype.l2), value = TRUE))
# table(tcell$predicted.celltype.l2)
# tcell <- runUMAP(tcell, features = VariableFeatures(object = tcell))
t_markers <- c("CD3D", "CD3E", "CD3F", "CD3G", "CD4", "CD8A", "CD8B")
function_marker <- Checkpoint
heatmap <- DoHeatmap(tcell, 
                    features = append(t_markers, function_marker),
                    size = 5, slot = "data",
                    label = TRUE)+
  labs(fills = "Cluster", colours = "Expression") + 
  theme(legend.position= "bottom", 
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 13))+ 
  scale_fill_distiller(palette = "RdBu")
heatmap
save_png(heatmap, "checkpoint_bm", 14, 10)
# DefaultAssay(tcell) <- "SCT"

tcell <- AddModuleScore(tcell, features = TH1_differ, 
                        ctrl = 5, assay = "SCT", name = "TH1_differ")
tcell <- AddModuleScore(tcell, features = intersect(IFNG_production, rownames(tcell)), 
                        ctrl = 5, assay = "SCT", name = "IFNG")
tcell <- AddModuleScore(tcell, features = cytotoxicity, 
                        ctrl = 5, assay = "SCT", name = "cytotoxic")
tcell <- AddModuleScore(tcell, features = Checkpoint, 
                        ctrl = 5, assay = "SCT", name = "checkpoint")
tcell <- AddModuleScore(tcell, features = Naive_T, 
                        ctrl = 5, assay = "SCT", name = "Naive_T")


tcell@meta.data["TH1_differ1"]

th1 <- FeaturePlot(tcell,
            features = "TH1_differ1", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
  labs(title = "TH1/2 differentiation")

IFNG <- FeaturePlot(tcell,
                    features = "IFNG1", label = FALSE , repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = 'IFNG production')

cytotoxic <- FeaturePlot(tcell, features = "cytotoxic1", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = "Cytotoxicity")

checkpoint <- FeaturePlot(tcell, features = "checkpoint1", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = "Checkpoint")

naive_T <- FeaturePlot(tcell, features = "Naive_T1", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = "Naive T cell markers")

cell.annot <- DimPlot(tcell, label.size = 3) 
cell.annot


plot_grid(cell.annot, th1, cytotoxic, IFNG, checkpoint, naive_T, 
             nrow = 3, rel_widths = c(1.5,1,1,1,1,1))



#===============================================================================
#GSEA analysis
#===============================================================================
sobj_gsea <- readRDS("C:/Users/jtao/work_dir/AML/data/GSEA/GSEA_SCT_BM_Tcell_8000.rds")
aml_gsea <- readRDS("C:/Users/jtao/work_dir/AML/data/GSEA/GSVA_Tcell_SCT_72_78.rds")

aml_gsea <- RECODE_Ontology(aml_gsea)
sobj_gsea <- RECODE_Ontology_bm(sobj_gsea)
Idents(sobj_gsea) <- sobj_gsea$celltype.l2
sobj_gsea <- RunUMAP(sobj_gsea, features = VariableFeatures(object = sobj_gsea))
sobj_gsea <- RECODE_Ontology_bm(sobj_gsea)

DefaultAssay(sobj_gsea) <- "SCT"
DefaultAssay(aml_gsea) <- "SCT"
sobj_gsea <- subset(sobj_gsea, ident = grep("T", unique(sobj_gsea$celltype.l2),
                                            value = TRUE))
aml_gsea <- subset(aml_gsea, ident =grep("T", unique(aml_gsea$predicted.celltype.l2),
                                         value = TRUE)) 

th1 <- FeaturePlot(sobj_gsea,
                   features = "BIOCARTA_TH1TH2_PATHWAY", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
  labs(title = "TH1/2 markers")

IFNG <- FeaturePlot(sobj_gsea,
                    features = "IFNG.6.Merk", label = FALSE , repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = 'IFNG production')

cytotoxic <- FeaturePlot(sobj_gsea, features = "Cytotoxicity", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = "Cytotoxicity")

checkpoint <- FeaturePlot(sobj_gsea, features = "Checkpoint.expression", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = "Checkpoint")

naive_T <- FeaturePlot(sobj_gsea, features = "Naive", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = "Naive T cell markers")

cell.annot <- DimPlot(sobj_gsea, label.size = 3)



plot_grid(cell.annot, th1, cytotoxic, IFNG, checkpoint, naive_T, 
          nrow = 3, rel_widths = c(1.5,1,1,1,1,1))

Idents(aml_gsea) <- aml_gsea$predicted.celltype.l2  
vln_aml <- Seurat::VlnPlot(aml_gsea, features = "Checkpoint.expression")
vln_bm <- VlnPlot(sobj_gsea, features = "Checkpoint.expression")
save_png(vln_bm, "TH1_bm_vln", 8, 6)

vln_aml
vln_bm
ident.order <- c("CD4 naive T cell", "CD8 Naive T cell", "effector CD8 T cell",
                 "effector memory CD4 T cell", "effector memory CD8 T cell", 
                 "gdT", "MAIT", "Treg")
aml_gsea$predicted.celltype.l2 <- factor (x = aml_gsea$predicted.celltype.l2,
                                          levels = ident.order)
sobj_gsea$celltype.l2 <- factor (x = sobj_gsea$celltype.l2,
                                          levels = ident.order)


