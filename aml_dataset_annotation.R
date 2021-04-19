library(Seurat)
library(dplyr)



# Sobj1<- readRDS("C:/Users/jtao/work_dir/AML/data/AMLnature/809653.seurat.rds")
# Sobj1 <- UpdateSeuratObject(Sobj1)
# 
# 
# mutate_file_path <- "C:/Users/jtao/work_dir/AML/data/AMLnature/scrna_mutations/scrna_mutations/cb_sniffer_results/809653_CellCoverage_Coding_180813.txt"


#=========================================================================
#helper function
#=========================================================================

mutation_function <- function(cell_ID, mutate_cell){
  if (cell_ID %in% mutate_cell){
    mutate_status = "mutation detected"
    
  }else{mutate_status = "Unknown"}
  return (mutate_status)
}

gene_function <- function(cell_ID, mutate_cell, mutation){
  if(cell_ID %in% mutate_cell){
     cell<- mutation %>% filter(Cell == cell_ID)
     if (dim(cell)[[1]] >=2){
       mutate_gene <- paste(cell$row.names, collapse = ",")
     }else{mutate_gene <- cell$row.names}
  }else {mutate_gene <- "no mutation detected"}
  return (mutate_gene)
}

#======================================================
# main functions 
#======================================================

adding_mutation_info <- function(Sobj1, alignment){
  mutation <- alignment %>% filter(AltReads > 0)
  mutate_cell <- as.vector(unique((mutation$Cell)))
  Sobj1$cell.ID <- rownames(Sobj1@meta.data)
  # mutate_cell
  cell.ID <- as.vector(Sobj1$cell.ID)
  Sobj1$mutation_annotation <- unlist(lapply(cell.ID, mutation_function,
                                             mutate_cell = mutate_cell))
  Sobj1$mutated_gene <- unlist(lapply(cell.ID, gene_function,
                                      mutate_cell = mutate_cell,
                                      mutation = mutation))
  
  return (Sobj1)
}


# Sobj1$mutation_annotation <- unlist(lapply(cell.ID, mutation_function))
# Sobj1$mutated_gene <- unlist(lapply(cell.ID, gene_function))
# mutation <- alignment %>% filter(AltReads > 0)
# mutate_cell <- as.vector(unique((mutation$Cell)))
# Sobj1$cell.ID <- rownames(Sobj1@meta.data)
# # mutate_cell
# cell.ID <- as.vector(Sobj1$cell.ID)

# Sobj1_mut <- subset(Sobj1, ident = mmutate_cell)
# saveRDS(Sobj1_mut, "C:/Users/jtao/work_dir/AML/data/AMLnature/mut/809653_mut.seurat.rds") 

sample_ID <- c("508084", "548327", "782328", "809653", "721214")
seurat_list <- list()
for (i in 1:length(sample_ID)){
  rds_path <- paste("/home/tjinyan/work_dir/AML/data/classifier/AMLnature/",
                    sample_ID[i], ".seurat.rds", sep = "")
  mut_anno_path <- paste("/home/tjinyan/work_dir/AML/data/classifier/AMLnature/cb_sniffer_results/",
                         sample_ID[i], "_CellCoverage_Coding_180813.txt", sep = "")
  output_path <- paste("/home/tjinyan/work_dir/AML/data/classifier/AMLnature/",
                       sample_ID[i], ".anno.rds", sep = "")
  
  alignment <- read.table(mut_anno_path, row.names = NULL)
  Sobj1 <- readRDS(rds_path)
  Sobj1 <- UpdateSeuratObject(Sobj1)
  
  Sobj1<- adding_mutation_info(Sobj1, alignment)
  saveRDS(Sobj1, output_path)
  seurat_list[[i]] <- Sobj1
}

seurat.big <- merge(unlist(seurat_list), add.cell.ids = sample_ID, 
                    project = "AML.mutation.annotation")

saveRDS(seurat.big, "/home/tjinyan/work_dir/AML/data/classifier/AMLnature/AML.mutation.annotation.rds")

