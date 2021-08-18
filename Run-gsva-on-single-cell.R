
rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape)
library(corrplot)
library(Matrix)
library(singlecell)
library(Seurat)
library(SingleR)
library(tidyverse)
library(cowplot)
library(GSVA)
library(GSEABase)
library(msigdbr)


#library(celldex)

setwd("C:/data/10x datasets/Seurat objects")
source("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/SingleCellSeuratBased/scRNAseq-unimodal-Seurat-based-functions.R")

SeuratObj <- readRDS("C:/data/10x datasets/Seurat objects/GSE126310_CTCL_multimodal-Seurat-celltypes.rds")

data.mat <- as.matrix(SeuratObj[["SCT"]]@data)

signatureList <- getGmt('/home/tjinyan/work_dir/AML/data/Signatures-for-AML.gmt')
gsva_es <- gsva(data.mat,gset.idx.list=signatureList)
gsvadf <- data.frame(Cell_id=colnames(gsva_es), t(gsva_es))
rownames(gsvadf) <-gsvadf$Cell_id
df.annotation <- SeuratObj@meta.data
dfGSVA <- merge(df.annotation, gsvadf, by=0)
rownames(dfGSVA)<- dfGSVA$Row.names
SeuratObj<- AddMetaData(SeuratObj, dfGSVA[,20:93])
saveRDS(SeuratObj,"C:/data/10x datasets/Seurat objects/GSE126310_CTCL_multimodal-Seurat-celltypes-gsva.rds")
