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




SeuratObj <- readRDS("/home/tjinyan/work_dir/AML/data/SCT_BM_Tcell_8000.rds")

data.mat <- as.matrix(SeuratObj[["SCT"]]@data)

signatureList <- getGmt('/home/tjinyan/work_dir/AML/data/genesets-T-cells.gmt')
gsva_es <- gsva(data.mat,gset.idx.list=signatureList)
gsvadf <- data.frame(Cell_id=colnames(gsva_es), t(gsva_es))
rownames(gsvadf) <-gsvadf$Cell_id
df.annotation <- SeuratObj@meta.data
dfGSVA <- merge(df.annotation, gsvadf, by=0)
rownames(dfGSVA)<- dfGSVA$Row.names
SeuratObj<- AddMetaData(SeuratObj, dfGSVA)
saveRDS(SeuratObj,"/home/tjinyan/work_dir/AML/data/classifier/AMLnature/GSEA_SCT_BM_Tcell_8000.rds")





signatureList <- getGmt('C:/Users/jtao/work_dir/AML/data/GSEA/genesets-T-cells-modified.gmt')
