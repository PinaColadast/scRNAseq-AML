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

# global_variables --------------------------------------------------------
project<-"AML-bmcite"
Sys.setenv(language="en")
if (project == "AML-bmcite"){
  working_dir<-"C:/Users/jtao/work_dir/AML/data"
  graph_dir<- "C:/Users/jtao/work_dir/AML/data/images/"
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}

print(working_dir)
setwd(working_dir)
getwd()

targets <- c("CD70", "CD33", "IL3RA", "CD123", "GP67")
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


my_DimPlot <- function(SeuratObject, reduction){
  pdim<- DimPlot(SeuratObject, reduction = reduction,
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


#==========================================================================
# clustering 
#==========================================================================
bm<- readRDS(paste(working_dir, "/bmcite6000.SCT.rds", sep = ""))

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- RunUMAP(bm, features = VariableFeatures(object = bm))
bm <- RunTSNE(bm, nn.name = "weighted.nn", reduction.name = "wnn.tsne",
              reduction.key = "wnnTSNE_")

Seurat::DefaultAssay(bm) <- "SCT"
#for multi-model data, use function

bm<- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

bm <- FindNeighbors(bm)
?FindNeighbors
bm <- FindClusters(bm, graph.name = "SCT_snn", resolution = 0.5,
                   algorithm = 3,
                             verbose = TRUE, )
?FindClusters
bm <- FindClusters(bm, graph.name = "wsnn",
                   algorithm = 3, resolution = 0.5,
                   verbose =TRUE)


Idents(bm) <-bm$seurat_clusters
cluster_no <- length(unique(bm$seurat_clusters))
cluster.markers.list <- list()
for (i in 0:(cluster_no-1)){
  j = i+1
  cluster.markers.list[[j]] <- FindMarkers(bm, ident.1 = i, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
}

cluster.markers.list

# l1 <- as.data.frame(table(bm@meta.data[["celltype.l1"]]))
# l2 <- as.data.frame(table(bm@meta.data[["celltype.l2"]]))
# 
# pl1 <- ggplot(data = l2, aes(x = Var1, y = Freq)) + 
#   geom_bar(stat="identity", width = 0.5,)+ coord_flip() +
#   # scale_fill_manual(values = c("grey59", "grey59", "grey59","skyblue3")) +
#   theme (legend.position = "none", aspect.ratio = 1.5/1, axis.text = element_text(size = 12),
#   ) + labs(y = "Cell Counts", x = "cell type by Seurat") + 
#   theme(panel.grid.minor = element_blank(),
#         panel.background = element_blank()) 


par(mar = c(1, 10, 4, 3))
Idents(bm) <- bm$SCT_snn_res.0.5
p1 <- my_DimPlot(bm, reduction = "umap")
Idents(bm) <- bm$cell_types_CL
p2 <- my_DimPlot(bm, reduction = "umap")

plot_grid(p1, p2)

cell_typel2 <- as.vector(unique(bm$celltype.l2))
cell_ontology <- as.data.frame(cell_typel2, cell_name)
write.table(cell_ontology, "cell_ontology.txt", sep = ",")
cell_typel2
cell_name <- c("progenitor RBC", "gamma-delta T cell", "CD4 T cell",
               "CD4 T cell", "monocyte(myeloid cell)", "B cell",
               "CD8 T cell", "regulatory T cell", "CD8 T cell", "NK cell",
               "granulocyte monocyte progenitor cell(myeloid cell)", 
               "CD8 T cell", "monocyte(myeloid cell)", "Plasmacytoid dendritic cell", 
               "CD8 T cell", "MAIT", "B cell", "dendritic cell",
               "NK cell", "pro-B cell", "megakaryocyte progenitor cell(myeloid cell)",
               "CD8 T cell", "Plasmablast", "hematopoietic stem cell",
               "early lymphoid progenitor", "macrophage dendritic cell progenitor(myeloid cell)",
               "pro-B cell")

s_cell_name <- c("progenitor RBC", "gamma-delta T cell", "CD4 T cell",
                "CD4 T cell", "myeloid cell", "B cell",
                "CD8 T cell", "regulatory T cell", "CD8 T cell", "NK cell",
                "myeloid cell", 
                "CD8 T cell", "myeloid cell", "Plasmacytoid dendritic cell", 
                "CD8 T cell", "MAIT", "B cell", "dendritic cell",
                "NK cell", "pro-B cell", "myeloid cell",
                "CD8 T cell", "Plasmablast", "hematopoietic stem cell",
                "early lymphoid progenitor", "myeloid cell",
                "pro-B cell" )

match_cell_name <- function(celltypel2, c_cell_typel2, cell_names){
  idx <- match (celltypel2, c_cell_typel2)
  cell_type <- cell_names[[idx]]
  
  return(cell_type) 
}

cell_types_curated <- c()
for (i in bm$celltype.l2){
  cell_types_curated <- append(cell_types_curated, match_cell_name(i, cell_typel2, cell_name))
  }
  

cell_concise <- c()
for (i in bm$celltype.l2){
  cell_concise <- append(cell_concise, match_cell_name(i, cell_typel2, s_cell_name))
}

bm$cell_types_general <- cell_concise
bm$cell_types_CL <- cell_types_curated
Idents(bm) <- bm$cell_types_CL

data<- c()
#=========================================================================
#data visualization 
#=========================================================================
Idents(bm) <- bm$cell_types_general
myeloid <- subset(bm, ident = c("myeloid cell"))
DefaultAssay(HSC.bm) <- "RNA"
heatmap <- DoHeatmap(HSC.bm, 
          features = targets,
          ,size = 3.5, angle = 30, slot = "data")+
          labs(fills = "Cluster", colours = "Expression") + 
          theme(legend.position= "bottom", 
                legend.text = element_text(size = 9))+ 
          scale_fill_distiller(palette = "RdBu")
heatmap
save_png(heatmap, "heatmap_adt_ml", 12, 5)

heatmap
rownames(myeloid@assays[["ADT"]]@data)
cd123adt <- as.vector(myeloid@assays[["ADT"]]@data["CD123", ])
cd123sct <- as.vector(myeloid@assays[["SCT"]]@data["IL3RA", ])
df <- data.frame(cd123adt, cd123sct)
colnames(df) <- c("protein expression", "gene expression")

scatter <- ggscatter(df, x= "protein expression", 
                                       y= "gene expression",
                                       add = "reg.line", cor.coef = TRUE, cor.method = "pearson",
                                       main = "CD123/IL3RA expression",
                                       xlab = "protein expression",
                                       ylab = "gene expression")

scatter_scaled <- ggscatter(df, x= "protein expression", 
                y= "gene expression",
                add = "reg.line", cor.coef = TRUE, cor.method = "pearson",
                main = "CD123/IL3RA expression",
                xlab = "scaled protein expression",
                ylab = "scaled gene expression")

combined<- plot_grid(scatter, scatter_scaled)
save_png(combined, "scatter_cd123_scaled_unscaled", 10, 8)
scatter



