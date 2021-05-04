rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

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
library(optparse)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
working_dir <- args[1]
filename <- args[2]

data_dir <- paste(working_dir, filename, sep = "/")
SeuratObj <- readRDS(data_dir)

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


bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
# DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap") 



           
            


anchors <- FindTransferAnchors(
  reference = bm,
  query = SeuratObj,
#   normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
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

SeuratObj<- Seurat.STnorm.pca(SeuratObj)
SeuratObj <- RECODE_Ontology(SeuratObj)

saveRDS(SeuratObj, paste(working_dir, "/SCT_", filename, sep = ""))

