library(Seurat)
library(stringr)
# for a1
a1 <- read.table("328_707_expression_mat.txt", sep = ",", header = TRUE,
                 row.names = 1)
a1.anno <- read.table("328_707_annotation.txt", sep = ",",
                      header = TRUE, row.names = 2)





rownames(a1.anno) <- str_replace_all(rownames(a1.anno), "-", ".")


S1 <- CreateSeuratObject(a1)
S1@meta.data <- transform(merge(S1@meta.data, a1.anno, by = 0, all= TRUE),
                          row.names = Row.names, Row.names = NULL)


bm.anno <- read.table("~/work_dir/AML/data/classifier/BM_annotation.txt",
                      sep = "\t", header = TRUE, row.names = 2)
bm <- read.table("~/work_dir/AML/data/classifier/BM_expression_mat.txt",
                 sep = "\t", header = TRUE, row.names = 1)
Sbm <- CreateSeuratObject(bm)

rownames(bm.anno) <- str_replace_all(rownames(bm.anno), "-", ".")
Sbm@meta.data <- transform(merge(Sbm@meta.data, bm.anno, by = 0, all= TRUE),
                          row.names = Row.names, Row.names = NULL)


# str_replace <- function(string){
#   new_str <- str_replace(string, "-", ".")
#   return(new_str)
#   
# }
# rownames <- unlist(lapply(rownames(a1.anno), str_repalce))

#for a2
a2 <- read.table("D0_expression_mat.txt", sep = ",",
                 header = TRUE, row.names = 1)

a2.anno <- read.table("D0_annotation.txt",  sep = ",",
                      header = TRUE, row.names = 2)

rownames(a2.anno) <- str_replace_all(rownames(a2.anno), "-", ".")
S2 <- CreateSeuratObject(a2)
S2@meta.data <- transform(merge(S2@meta.data, a2.anno, by = 0, all= TRUE),
                          row.names = Row.names, Row.names = NULL)


S37_bm <- merge(S1, Sbm, add.cell.ids = c("AML", "bm_norm"))
SDO_bm <- merge(S2, Sbm, add.cell.ids= c("AML_D0", "bm_norm"))