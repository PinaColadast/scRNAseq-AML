# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GSEABase")
# BiocManager::install("GSVAdata")
# BiocManager::install("GSVA")
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
library(ggpubr)
library(ggplot2)
library(comprehenr)
# gmt <- GSA.read.gmt('C:/Users/jtao/work_dir/AML/data/GSEA/Signatures-for-AML.gmt')
# gmt <- getGmt('C:/Users/jtao/work_dir/AML/data/GSEA/Signatures-for-AML.gmt')
# 
# gmt
# 
# install.packages("GSVA")
# install.packages("GSEABase")
# install.packages("GSVAdata")
# 
# install_github("GSEA-MSigDB/GSEA_R")
# gmt[1]
# # source(system.file('extdata', 'Run.GSEA.R', package = 'GSEA'))
#===============================================================================
#helper function
#===============================================================================
my_DimPlot <- function(SeuratObject){
  pdim<- DimPlot(SeuratObject, reduction = "umap",
                 label = TRUE, 
                 # cols = MyColorsCellTypes,
                 pt.size = 1, 
                 repel = TRUE,
                 label.size = 4)+theme_void()+theme( legend.position="bottom")
  return (pdim)}

#===============================================================================
#global variables
#===============================================================================

project<-"GSEA"
Sys.setenv(language="en")
if (project == "GSEA"){
  working_dir<-"C:/Users/jtao/work_dir/AML/data/GSEA"
  graph_dir<- "C:/Users/jtao/work_dir/AML/data/GSEA/images/"
  setwd(working_dir)
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}

AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                  "granulocyte monocyte progenitor", "hematopoietic stem cell",
                  "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")

sobj_gsea <- readRDS("C:/Users/jtao/work_dir/AML/data/GSEA/GSVA_Tcell_SCT_72_78.rds")
sobj_gsea <- RECODE_Ontology(sobj_gsea)
sobj_gsea <- add.mut.CellType(sobj_gsea, reference = TRUE)

Idents(sobj_gsea) <- sobj_gsea$l2CellType.mut
Idents(sobj_gsea) <- sobj_gsea$predicted.celltype.l2
Sobj_stem <- subset(sobj_gsea, ident = AML_category)


DefaultAssay(sobj_gsea) <- "integrated"
sobj_gsea <- RunUMAP(sobj_gsea, features = VariableFeatures(sobj_gsea))
umap <- as.data.frame(sobj_gsea@reductions[["umap"]]@cell.embeddings)
pro_scatter <- ggplot (umap, 
                       aes(x= UMAP_1, y=UMAP_2, 
                       )) 



LSC17 <- pro_scatter + labs(colour = "cytotoxic") + 
  theme_void()+ 
  theme(legend.position = "bottom") +
  geom_point(size = 0.5, aes(color = sobj_gsea$Cytotoxicity)) + 
  scale_color_gradientn(colours = rainbow(5))

cyto <- FeaturePlot(sobj_gsea, features = "Cytotoxicity", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(title = "cytotoxicity GSEA")

Idents(sobj_gsea) <- sobj_gsea$predicted.celltype.l2
cyto
Idents(Sobj_stem) <- Sobj_stem$l2CellType.mut
LSC17

dimplot <- my_DimPlot(sobj_gsea)
dimplot
save_png(dimplot, "dimplot_50_stem", 11, 8)
plot_grid(LSC17, dimplot)
LSC17
save_png(LSC17, "24006_HIGHLSCvsHSC", 8,8)
dimplot

fontsize <- 7
df_meta<- Sobj_stem@meta.data
df_score<- df_meta[,c("GSE24006_HighLSCvsHSC", "l2CellType.mut")]
df_score
df_score$ident <- df_meta$predicted.celltype.l2
df_score$mut <- to_vec(for(i in df_score$l2CellType.mut) if (str_detect(i, "malig")) "malignant" else "normal")
df_cnt <- as.data.frame(table(df_score$l2CellType.mut))
df_cnt

head(df_score, 1, 10)
names(df_score) <- c("HighLSCvsHSC score", "Cell.mut", "Cell", "malignancy")

Gscore <- as.vector(df_score$`HighLSCvsHSC score`)
lab_y <- max(Gscore) * (max(Gscore)/(max(Gscore) - min(Gscore)) - 0.05)


give.n <- function(x){
  return(data.frame(y = max(x)-0.5, label = paste0("n = ", length(x))))
}

p1 <- ggboxplot(df_score, x = "malignancy", y = "HighLSCvsHSC score",
                color = "malignancy", size  = 0.1,
                add = "jitter", palette=  c("#00AFBB", "#E7B800"),
                ylab = "GSEA score", xlab = "",
                main = paste("Subject:", 508084, "\n", "GSE24006_HighLSCvsHSC GSEA score", sep = ""), 
                facet.by ="Cell", short.panel.labs = FALSE,
                font.label = list(size = 14, face = "bold")) +
  stat_compare_means(label = "p.format", size = fontsize-2,
                     family = "mono", 
                     label.x = 1.3, label.y = lab_y )+
  stat_summary(fun.data = give.n, geom = "text", size = fontsize-2) + 
  theme(text = element_text(size = fontsize + 3),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"),
                                    face = "bold", size = fontsize ),
        axis.text.x  = element_text(size = fontsize),
        legend.text = element_text(size = fontsize)) + 
  font("legend.title",  face = "bold", size = fontsize+2)+
  font("legend.text", size = fontsize +2)+
  font("xy.text", size = fontsize +2 )+
  font("title", size = fontsize +2 )


p1 + facet(p1, facet.by = "Cell", 
           panel.labs.font = list(size = fontsize +2))+ 
  facet_wrap( ~ variable, scales="free")

p1
