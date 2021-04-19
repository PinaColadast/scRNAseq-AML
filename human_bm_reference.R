rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape)
library(corrplot)
library(Matrix)
library(Seurat)
# library(SingleR)
library(tidyverse)
library(cowplot)
library(sctransform)
library(SeuratDisk)
library(SeuratData)
library(patchwork)


InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
bm <- ScaleData(bm, assay = 'RNA')
bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')

bm <- FindNeighbors(
  object = bm,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)


data_dir <- "/home/tjinyan/work_dir/AML/data/classifier/mut/AMLnature_train.rds"
#==========================================================================
#read in mtx 
#==========================================================================
# Read in `matrix.mtx`
# counts <- readMM(paste(data_dir,
#                        "matrix.mtx.gz",
#                        sep= "/")
# )

# # Read in `genes.tsv`
# genes <- read_tsv(paste(data_dir,
#                         "features.tsv.gz",
#                         sep= "/"), col_names = FALSE)
# gene_ids <- genes$X2

# Read in `barcodes.tsv`
# cell_ids <- read_tsv(paste(data_dir, "barcodes.tsv.gz", sep = "/"),
#                      col_names = FALSE)$X1

# rownames(counts) <- gene_ids
# colnames(counts) <- cell_ids
SeuratObj <- readRDS(data_dir)

# DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()




#recode Azimuth names to cell ontology

# p1 = DimPlot(SeuratObj, reduction = "ref.umap", group.by = "cell.type", label = TRUE, label.size = 3, repel = TRUE)# + NoLegend()
# p2 = DimPlot(SeuratObj, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3 ,repel = TRUE) #+ NoLegend()
# p1 + p2

#saveRDS()




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

SeuratObj<- Seurat.STnorm.pca(SeuratObj)

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
# DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap") 



           
            


anchors <- FindTransferAnchors(
  reference = bm,
  query = SeuratObj,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
  reference.neighbors = "spca.annoy.neighbors"
)
# for non-multimodal data

SeuratObj <- MapQuery(
  anchorset = anchors,
  query = SeuratObj,
  reference = bm,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

saveRDS(SeuratObj, "~/work_dir/AML/data/classifier/mut/AMLnature_train.SCT.rds")
