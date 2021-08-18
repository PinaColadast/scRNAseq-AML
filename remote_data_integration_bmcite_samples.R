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
# library(SingleCellExperiment)
library(scater)
library(SingleR)
library(patchwork)
library(ggpubr)
library(comprehenr)
# library(EnhancedVolcano)

# remoter::client(port = 55555, addr = "10.100.0.50")
project<-"remote-AML"
Sys.setenv(language="en")
if (project == "remote-AML"){
  working_dir<-"/home/tjinyan/work_dir/AML/data/"
  graph_dir<- "/home/tjinyan/work_dir/AML/data/images"
  setwd(working_dir)
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}


getwd()
Sobj <- readRDS("/home/tjinyan/work_dir/AML/data/classifier/AMLnature/SCT_721214.anno.rds")
Sobj1 <- readRDS("/home/tjinyan/work_dir/AML/data/classifier/AMLnature/SCT_782328.anno.rds")
#============================================================
#bmcite =======================================================
#=============================================================
  
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
# bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
#               reduction.key = "wnnUMAP_", return.model = TRUE)
# bm <- ScaleData(bm, assay = 'RNA')
# # bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')
# bm <- RunPCA(bm, assay = "RNA", graph = "wsnn")
S_list <- list(bm)

# bm <- FindNeighbors(
#   object = bm,
#   reduction = "spca",
#   dims = 1:50,
#   graph.name = "spca.annoy.neighbors", 
#   k.param = 50,
#   cache.index = TRUE,
#   return.neighbor = TRUE,
#   l2.norm = TRUE
# )

Seurat.STnorm.pca <- function(SeuratObj){
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-")
  SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SeuratObj.ST <- SCTransform(SeuratObj, assay = "RNA",vars.to.regress =c( "percent.mt", "S.Score", "G2M.Score" ), return.only.var.genes = FALSE)
  SeuratObj.ST <- FindVariableFeatures(SeuratObj.ST, nfeatures = 8000) 
  SeuratObj.ST <- ScaleData(SeuratObj.ST)
  SeuratObj.ST <- RunPCA(SeuratObj.ST)
  SeuratObj.ST <- RunUMAP(SeuratObj.ST, 
                          features = VariableFeatures(object = SeuratObj.ST))
  
  DefaultAssay(SeuratObj.ST) <- "SCT"
  
  return(SeuratObj.ST)
}


#
# S_list <- lapply(X = S_list, FUN = function(x) {
#   x <- Seurat.STnorm.pca(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 8000)
# })
# 

bm <- Seurat.STnorm.pca(bm)
S.list <- list(bm, Sobj, Sobj1)
features <- SelectIntegrationFeatures(object.list = S.list, nfeatures = 8000)

S.anchors <- FindIntegrationAnchors(object.list = S.list, 
                                    anchor.features = features,
                                    reduction = "rpca",
                                    k.anchor = 20)

S.combined <- IntegrateData(anchorset = S.anchors)
DefaultAssay(S.combined) <- "integrated"
S.combined <- ScaleData(S.combined, verbose = FALSE)

