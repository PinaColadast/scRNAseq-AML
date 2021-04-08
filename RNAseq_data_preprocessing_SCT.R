#raw data processing 
#based on RNA transcriptomic matrix
args = commandArgs(trailingOnly = TRUE)

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
library(patchwork)
library(cowplot)
library(optparse)

# option_list <- list(
#   make_option(c("-d", "--datadir"), action = "store",
#               type = "character",
#               help="matrix object filepath"),
#   make_option(c("-f", "--filename"), action= "store",
#               type = "character",
#               help = "file name")
# )

# parser <- OptionParser(option_list=option_list)

# args <- parse_args(parser, positional_arguments = 1)
# working_dir <- args$datadir
# file <- args$filename
# file_dir <- paste(working_dir, file, sep = "/")
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
working_dir <- args[1]
file <- args[2]
file_dir <- paste(working_dir, file, sep = "/")


if (file.exists(file_dir)==FALSE){
  print(file_dir)
  stop("input file doesn't exist, please show correct file")
}

# global_variables --------------------------------------------------------
setwd(working_dir)
working_dir <- getwd()
Sys.setenv(language="en")


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
  SeuratObj.ST <- FindVariableFeatures(SeuratObj.ST, nfeatures = 8000) 
  SeuratObj.ST <- ScaleData(SeuratObj.ST)
  SeuratObj.ST <- RunPCA(SeuratObj.ST)
  
  # DefaultAssay(SeuratObj.ST) <- 'ADT'
  # # we will use all ADT features for dimensional reduction
  # # we set a dimensional reduction name to avoid overwriting the 
  # # VariableFeatures(SeuratObj.ST) <- rownames(SeuratObj.ST[["ADT"]])
  # SeuratObj.ST <- NormalizeData(SeuratObj.ST, assay = "ADT",normalization.method = 'CLR', margin = 2) %>% 
  #   ScaleData() %>% RunPCA(reduction.name = 'apca')
  return(SeuratObj.ST)
}

#=============================================================================
# load expression matrix, create Seurat object
#=============================================================================
mtx <- read.table(file_dir, row.names=1, header = TRUE, sep = ",")
Sobj <- CreateSeuratObject(mtx, min.cells = 3, min.features = 200)

#==========================================================================
# SCT transformation  
#==========================================================================

Sobj <- Seurat.STnorm.pca(Sobj)
out_name <- strsplit(file, ".", fixed = TRUE)[[1]][[1]]
out_dir <- paste(working_dir,"/output", sep ="")

if (dir.exists(out_dir)){
  out_dir <- out_dir
} else {dir.create(out_dir)}

saveRDS(Sobj, paste(working_dir,"/output/", "SCT_", out_name, ".rds", sep =""))
