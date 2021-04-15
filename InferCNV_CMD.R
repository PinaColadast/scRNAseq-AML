args = commandArgs(trailingOnly = TRUE)


library(infercnv)
library(biomaRt)
library(dplyr)
library(Seurat)
library(optparse)


#=====================================================================#
# arguments parser set-up
#=====================================================================#

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
working_dir <- args[1]
mat_file <- args[2]
anno_file <- args[3]
mat_dir <- paste(working_dir, mat_file, sep = "/")
anno_dir <- paste(working_dir, anno_file, sep = "/")


if (file.exists(mat_dir)==FALSE){
  print(file_dir)
  stop("input file doesn't exist, please show correct file")
}
#set working directory
setwd(working_dir)
working_dir <- getwd()
Sys.setenv(language="en")
#----------------------------------------------------------------------------
# #Create Infercnv Object # -----------------------------------------------
#----------------------------------------------------------------------------

# 1. raw_counts matrix 
# extract from Seurat object 
# Read rds dataset
data <- read.table(mat_dir, sep = "\t", header = TRUE, row.names=1, quote="")

raw_counts_matrix <- as.matrix(data)

# 2. cell annotation files

#made already

#------------------------------------------------------
# 3. gene order file (order on Chromosomes)
#------------------------------------------------------
all.genes <- c(rownames(data))

df_gene <- as.data.frame(table(all.genes), row.names = NULL)
colnames(df_gene) <- c("hgnc_symbol", "freq")
head(df_gene,10)
# retrieve chromosomes positions given a list of genes 
query_genes = df_gene
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results_id <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'chromosome_name',
                                   'start_position', 'end_position'),
                    filters = "hgnc_symbol", 
                    values = df_gene$hgnc_symbol, 
                    mart = ensembl)

chromo_list <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15", "16", "17", "18", "19",
                 "20", "21", "22", "X", "Y")

results <- results_id %>% filter(chromosome_name %in% chromo_list)
results <- results %>% select(hgnc_symbol, chromosome_name, start_position, end_position)

#check if any duplicates in gene position
rep_gene <- data.frame(table(results$hgnc_symbol))
results[results$hgnc_symbol %in% rep_gene[rep_gene$Freq>1, ]$Var1, ]

#clear replicates
results_unique <- results[!duplicated(results$hgnc_symbol), ]

# write table of gene notations
write.table(results_unique, paste(getwd(), "/output/gene_chromopos.txt", sep = ""),
            col.names = FALSE, row.names = FALSE,  sep = "\t", quote = FALSE)

# filter the counts matrix according to results of chromosome positions
counts_matrix <- raw_counts_matrix[c(results$hgnc_symbol), ]
# write.table(counts_matrix, file = paste(getwd(), "data/output/cnt_matrix", sep = "/")
#              , sep = "\t", col.names= FALSE, row.names = FALSE)


#-------------------------------------------------------------------------------
# Create InferCNV object and run ------------------------------------------
#-------------------------------------------------------------------------------
out_dir <- paste(getwd(), "/output/InferCNV/", sep = "")

if (dir.exists(out_dir)){
  out_dir <- out_dir
} else {dir.create(out_dir)}

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=anno_dir,
                                    delim="\t",
                                    gene_order_file= paste(getwd(), "output/gene_chromopos.txt", sep = "/"),
                                    ref_group_names=c("normal")) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

print(paste("files are written in:", out_dir, sep = " "))
      
