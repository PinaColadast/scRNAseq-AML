library("pivottabler")
library("reshape2")
library("tidyverse")
library("scales")
library("ggpubr")
# co-expression analysis (compact version)


#===============================================================================
#take data
#===============================================================================
sample_name <- c("508084"
                 , "548327", "721214", "782328")
AML_mut_category <- c("granulocyte monocyte progenitor-malignant", 
                  "hematopoietic stem cell-malignant",
                  "lymphoid-primed multipotent progenitor-malignant", 
                  "megakaryocyte progenitor cell-malignant",
                  "erythrocyte progenitor cell-malignant")
AML_category <- c("common dendritic progenitor", "erythroid progenitor cell",
                  "granulocyte monocyte progenitor", "hematopoietic stem cell",
                  "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")
AML_sub_category <- c("granulocyte monocyte progenitor", "hematopoietic stem cell",
                      "lymphoid-primed multipotent progenitor", "megakaryocyte progenitor cell")
co_markers <- c("CD33", "IL3RA")
# num <- 1
# for (i in sample_name){
# 
#   Sobj_path <- paste0("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_",
#                       i, ".anno.rds")
#   AML_nature <- readRDS(Sobj_path)
#   AML_nature <- RECODE_Ontology(AML_nature)
#   AML_nature <- add.mut.CellType(AML_nature, reference = TRUE)
# 
#   df_meta <- AML_nature@meta.data
#   cell_UMI <- rownames(df_meta %>% filter(l2CellType.mut %in% AML_mut_category))
#   co_expression_matrix <- as.data.frame(t(AML_nature@assays[["SCT"]]@counts[co_markers, cell_UMI]))
#   rownames(co_expression_matrix) <-c()
#   co_expression_matrix$Cell <- df_meta[cell_UMI, "l2CellType.mut"]
#   co_expression_matrix$sample <- i
# 
#   if (num == 1){
#     data_frame <- co_expression_matrix
#   }else{
#     data_frame <- rbind(data_frame, co_expression_matrix)
#   }
#   num <- num+1
#   AML_nature <- c()
#   }
# 
# # co_expression_matrix
# # AML_nature <- read.RDS("C:/Users/jtao/work_dir/AML/data/AMLnature/SCT_508084.anno.rds")
# dim(data_frame)
# write.table(data_frame, "data_co_expression.txt")

plot_co_expre<- function(Sobj){
  
  Idents(Sobj) <- Sobj$predicted.celltype.l2
  mutual_celltype <- intersect(AML_category, unique(Sobj$predicted.celltype.l2))
  Sobj <- subset(Sobj, ident = mutual_celltype)
  
  # Idents(Sobj) <- Sobj$l2CellType.mut
  # Sobj <- subset(Sobj, ident = grep("malignant", unique(Sobj$l2CellType.mut), value = TRUE))
  # 
  subject <- unique(Sobj$orig.ident)
  co_markers<- c("CD33", "IL3RA")
  df_meta <- Sobj@meta.data
# cell_UMI <- rownames(df_meta %>% filter(l2CellType.mut == celltype))
  df<- Sobj@assays[["SCT"]]@counts[co_markers, ]

  df_aa <- data.frame("CD33" = as.vector(df["CD33", ]),
                    "IL3RA" = as.vector(df["IL3RA", ]),
                    "Cell" = as.vector(df_meta$predicted.celltype.l2))

# df_aa <- df_aa[df_aa[, "CD33"] >0 & df_aa[, "IL3RA"] >0, ]
  df_aa$color <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "orchid",
                      ifelse(df_aa$CD33>0.5, "skyblue3"
                             , ifelse(df_aa$IL3RA>0.5, "indianred3", "grey59")) )
  df_aa$label <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "CD33 > 0.5, IL3RA > 0.5",
                      ifelse(df_aa$CD33>0.5 & df_aa$IL3RA<0.5, "CD33 > 0.5, IL3RA < 0.5"
                             , ifelse(df_aa$IL3RA>0.5, "CD33 < 0.5, IL3RA > 0.5", "CD33 < 0.5, IL3RA < 0.5")) )

# head(df_aa, 5,5)

  pivot_table <- dcast(df_aa, label ~ Cell, fun.aggregate = length, value.var = "label")
  melted <- melt(pivot_table)
  colnames(melted) <- c("Expression", "Cell", "Cell Counts")

  fontsize = 8

  p1 <- ggbarplot(melted, x = "Expression",y= "Cell Counts", fill = "Expression",
                color = "Expression", size  = 0.1,
                palette=c("grey59", "indianred3", "skyblue3", "orchid"),
                ylab = "Counts", xlab = "",
                main = paste("Subject: ", subject, "\n", "coexpression of CD33 & IL3RA", sep = ""), 
                short.panel.labs = FALSE,
                font.label = list(size = 10, face = "bold"))+
  guides(fill = guide_legend(nrows = 2, byrow= TRUE))

  p2 <- facet(p1, facet.by = "Cell", panel.labs.font.x = list(size= 9, face = "bold")) +
  theme(text = element_text(size = fontsize +2),
        # strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"),
        # face = "bold", size = fontsize),
        legend.position = "bottom",
        # axis.text.x  = element_text(size = fontsize+2, angle = 90),
        axis.text.x  = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(size = fontsize+2)) + 
  font("legend.title",  face = "bold", size = fontsize+2)+
  font("legend.text", size = fontsize+2)+
  font("xy.text", size = fontsize+2)+
  font("title", size = fontsize +2) + rremove("x.text")
  return(p2)
}


plexp <- plot_co_expre(hcabm40k)
plexp 
save_png(plexp, "hcabm40k_CD33_IL3RA", 10, 6)
# bp <- ggplot(df_aa, aes(x=label, fill = color)) + geom_bar(width = 0.5, stat = "count") +
#   
#   scale_fill_manual(values = c("grey59", "indianred3", "orchid","skyblue3"),
#                     labels = c("0" = "Absent", "1" = "Present")) +
#   
#   theme (legend.position = "none", aspect.ratio = 1.5/1, axis.text = element_text(size = 10),
#          axis.title.y = element_blank()
#   ) + theme(aspect.ratio =1) + labs(y = "Cell Counts") +facet_wrap(~Cell)
# bp
unique(AML_sub$l2CellType.mut)

"IL3RA" %in% rownames(hcabm40k@assays[["SCT"]]@counts)

unique(hcabm40k@assays[["SCT"]]@counts["IL3RA", ])

#===============================================================================
# co_expression_multi_sample
#===============================================================================
coexp <- read.table("data_co_expression.txt")
coexp$color <- ifelse(coexp$CD33>0.5 & coexp$IL3RA >0.5, "orchid",
                      ifelse(coexp$CD33>0.5, "skyblue3"
                             , ifelse(coexp$IL3RA>0.5, "indianred3", "grey59")) )
coexp$label <- ifelse(coexp$CD33>0.5 & coexp$IL3RA >0.5, "CD33 > 0.5, IL3RA > 0.5",
                      ifelse(coexp$CD33>0.5 & coexp$IL3RA<0.5, "CD33 > 0.5, IL3RA < 0.5"
                             , ifelse(coexp$IL3RA>0.5, "CD33 < 0.5, IL3RA > 0.5", "CD33 < 0.5, IL3RA < 0.5")) )
coexp

# df_a <- coexp %>% filter(sample == "508084")
# df_a$label
plot_coexpression <- function(df, sample) {
  df_aa <- df %>% filter(sample == sample)
  sp <- ggplot(df_aa) 
  geom_point(aes(x = IL3RA, y = CD33, colour = label), size  = 2) + 
    scale_color_manual(values =c("grey59", "indianred3","skyblue3","orchid")) + 
    geom_hline(yintercept= 0.5, linetype = "dotted", color = "red", 
               size = 1)  +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size = 1) +
    labs(x = 'IL3RA Expression',y = "CD33 Expression", colour = "Expression Group", title = "") + 
    theme(axis.title= element_text(size = 13),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 11),
          aspect.ratio = 1/3) +coord_fixed()
    
  return(sp)
}

# aaa <- plot_coexpression(coexp, "508084")


sp <- ggscatter(df_a, x= "IL3RA", y ="CD33", color = "label")  + 
  scale_color_manual(values =c("grey59", "indianred3","skyblue3","orchid")) + 
  geom_hline(yintercept= 0.5, linetype = "dotted", color = "red", 
             size = 1)  +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size = 1) +
  labs(x = 'IL3RA Expression',y = "CD33 Expression", colour = "Expression Group", title = "") + 
  theme(axis.title= element_text(size = 8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        aspect.ratio = 1/1) +coord_fixed()

ssp <- facet(sp, facet.by = "Cell", fnrow =1)
save_png(ssp, "508084.coexp", 10, 5)
cell.order <- c("granulocyte monocyte progenitor-malignant",
                "hematopoietic stem cell-malignant",
                "lymphoid-primed multipotent progenitor-malignant",
                "megakaryocyte progenitor cell-malignant")

sample.order <-c("508084", "548327", "721214", "782328")
var_width = 30

proportion <- c()
n_cells <-c()
for (i in sample.order){
  for (j in cell.order){
    sample_cell <- coexp %>% filter(Cell== j , sample == i)
    if (nrow(sample_cell %>% filter(label == "CD33 > 0.5, IL3RA > 0.5")) < 1){
      freq <- 0
    }else{
      freq <- prop.table(table(sample_cell$label))[["CD33 > 0.5, IL3RA > 0.5"]]
    }
    n_cells <- append(n_cells, nrow(sample_cell))
    proportion <- append(proportion, freq)
  }
}
n_cells
proportion
give.n <- function(x){
  return(data.frame(y = max(x)-0.5, label = paste0("n = ", length(x))))
}

coexp <- mutate(coexp, Cell.Name = str_wrap(Cell, width = var_width))
dat_sum <- data.frame(
  num_of_cells = n_cells,
  proportion = proportion,
  Cell.Name = c(cell.order,cell.order,cell.order,cell.order),
  sample = c("508084", "508084", "508084", "508084",
  "548327", "548327","548327","548327",
  "721214", "721214", "721214","721214",
  "782328","782328","782328","782328")
  
)
dat_sum <- mutate(dat_sum, Cell.Name = str_wrap(Cell.Name, width = var_width))

dat_sum$text <- sprintf("number of cells = %s\nCD33++IL3RA++ %s", format(dat_sum$num_of_cells, nsmall = 0), 
                                                  percent(dat_sum$proportion))


gsp <- ggplot(coexp)  + 
geom_point(aes(x = IL3RA, y = CD33, colour = label), size  = 2) + 
  scale_color_manual(values =c("grey59", "indianred3","skyblue3","orchid")) + 
  geom_hline(yintercept= 0.5, linetype = "dotted", color = "red", 
             size = 1)  +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size = 1) +
  labs(x = 'IL3RA Expression',y = "CD33 Expression", colour = "Expression Group", title = "") + 
  theme(axis.title= element_text(size = 8+4),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8+4),
        aspect.ratio = 1/1.2) +coord_fixed()+
  facet_grid(cols = vars(Cell.Name), rows = vars(sample), margins = FALSE,
             )+
  geom_text(check_overlap = TRUE,
    size    = 4.5,
    data    = dat_sum,
    mapping = aes(x = Inf, y = Inf, label = text)
    ,hjust   = 1.05,
    vjust   = 1.5
  )+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 13, face = "bold"),
        strip.text.y = element_text(size = 11, face = "bold")
        )
save_png(gsp, "facet", 18, 16)
gsp




#===============================================================================
#for bm
#===============================================================================

cell.order <- c("granulocyte monocyte progenitor",
                "hematopoietic stem cell",
                "lymphoid-primed multipotent progenitor",
                "megakaryocyte progenitor cell")
unique(hcabm40k$orig.ident)
hcabm40k <- RECODE_Ontology(hcabm40k)
mutual_celltype <- intersect(AML_sub_category, unique(hcabm40k$predicted.celltype.l2))
Idents(hcabm40k) <- hcabm40k$predicted.celltype.l2
bm_sub <- subset(hcabm40k, ident = mutual_celltype)

head(coexp)
df_meta <- bm_sub@meta.data

co_expression_matrix <- as.data.frame(t(bm_sub@assays[["RNA"]]@counts[co_markers,]))
rownames(co_expression_matrix) <-c()
co_expression_matrix$Cell <- df_meta$predicted.celltype.l2
df_aa<-co_expression_matrix
df_aa$color <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "orchid",
                      ifelse(df_aa$CD33>0.5, "skyblue3"
                             , ifelse(df_aa$IL3RA>0.5, "indianred3", "grey59")) )
df_aa$label <- ifelse(df_aa$CD33>0.5 & df_aa$IL3RA >0.5, "CD33 > 0.5, IL3RA > 0.5",
                      ifelse(df_aa$CD33>0.5 & df_aa$IL3RA<0.5, "CD33 > 0.5, IL3RA < 0.5"
                             , ifelse(df_aa$IL3RA>0.5, "CD33 < 0.5, IL3RA > 0.5", "CD33 < 0.5, IL3RA < 0.5")) )
var_width_bm <- 20
df_aa <- mutate(df_aa, Cell.Name = str_wrap(Cell, width = var_width_bm))

proportion_bm <- c()
n_cells_bm <-c()
for (i in cell.order){
  sample_cell <- df_aa %>% filter(Cell==i)
  if (nrow(sample_cell %>% filter(label == "CD33 > 0.5, IL3RA > 0.5")) < 1){
    freq <- 0
    }else{
      freq <- prop.table(table(sample_cell$label))[["CD33 > 0.5, IL3RA > 0.5"]]
    }
    n_cells_bm <- append(n_cells_bm, nrow(sample_cell))
    proportion_bm <- append(proportion_bm, freq)
    }


dat_sum_bm <- data.frame(
  num_of_cells = n_cells_bm,
  proportion = proportion_bm,
  Cell.Name = cell.order
)

dat_sum_bm <- mutate(dat_sum_bm, Cell.Name = str_wrap(Cell.Name, width = var_width_bm))

dat_sum_bm$text <- sprintf("number of cells = %s\nCD33++IL3RA++ %s", format(dat_sum_bm$num_of_cells, nsmall = 0), 
                        percent(dat_sum_bm$proportion))



gbmp <- ggplot(df_aa)  + 
  geom_point(aes(x = IL3RA, y = CD33, colour = label), size  = 2) + 
  scale_color_manual(values =c("grey59", "indianred3","skyblue3","orchid")) + 
  geom_hline(yintercept= 0.5, linetype = "dotted", color = "red", 
             size = 1)  +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size = 1) +
  labs(x = 'IL3RA Expression',y = "CD33 Expression", colour = "Expression Group", title = "") + 
  theme(axis.title= element_text(size = 8+4),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8+4),
        aspect.ratio = 1/1.2) +coord_fixed()+
  facet_grid(cols = vars(Cell.Name), margins = FALSE)+
  facet_wrap(~ Cell.Name, ncol=2)+
  geom_text(check_overlap = TRUE,
            size    = 4.5,
            data    = dat_sum_bm,
            mapping = aes(x = Inf, y = Inf, label = text)
            ,hjust   = 1.05,
            vjust   = 1.5
  )+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 12, face = "bold"),
  )
gbmp
save_png(gbmp, "RNA_hcabm40k", 10, 10)

