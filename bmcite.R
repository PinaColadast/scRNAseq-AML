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
library(cowplot)


# global_variables --------------------------------------------------------
project<-"AML-bmcite"
Sys.setenv(language="en")
if (project == "AML-bmcite"){
  working_dir<-"/home/tjinyan/work_dir/AML/data"
  graph_dir<- "/home/tjinyan/work_dir/AML/data/output/images"
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}

print(working_dir)

#==========================================================================
# helper functions 
#==========================================================================
Seurat.STnorm.pca <- function(SeuratObj){
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  SeuratObj <- CellCycleScoring(SeuratObj,
                                s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, 
                                                    pattern = "^MT-")
  
  
  
  SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SeuratObj.ST <- SCTransform(SeuratObj, assay = "RNA",vars.to.regress =c( "percent.mt", "S.Score", "G2M.Score" ), return.only.var.genes = FALSE)
  SeuratObj.ST <- FindVariableFeatures(SeuratObj.ST, nfeatures = 15000) 
  SeuratObj.ST <- ScaleData(SeuratObj.ST)
  SeuratObj.ST <- RunPCA(SeuratObj.ST)
  
  DefaultAssay(SeuratObj.ST) <- 'ADT'
  # we will use all ADT features for dimensional reduction
  # we set a dimensional reduction name to avoid overwriting the 
  # VariableFeatures(SeuratObj.ST) <- rownames(SeuratObj.ST[["ADT"]])
  SeuratObj.ST <- NormalizeData(SeuratObj.ST, assay = "ADT",normalization.method = 'CLR', margin = 2) %>% 
    ScaleData() %>% RunPCA(reduction.name = 'apca')
  return(SeuratObj.ST)
}



#==========================================================================
# data loading 
#==========================================================================
SeuratData::InstallData("bmcite")
bm <- LoadData(ds = "bmcite")


#==========================================================================
# SCT transformation  
#==========================================================================

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


bm <- Seurat.STnorm.pca(bm)
saveRDS(bm, paste(working_dir, "bmcite.SCT.rds", sep = "/"))





