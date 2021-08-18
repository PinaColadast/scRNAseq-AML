library(Seurat)
library(caret)
bmcite_pred <- read.table(file=paste0("bmcite", ".nonbi.predict.txt"), sep = "\t")
pred_nonbi <- read.table(file = paste0("bmcite", ".nonbi.predict.txt"),
                          sep = "\t")
bmcite_pred

table(bmcite_pred$celltype.l2)
aa<- sort(table(bmcite_pred$RF2.class))
table(bmcite_pred$RF2.class)

table(bmcite_pred$celltype.l2)
aa<- barplot(aa)
aa


#cell name curation

classifier_name <- c("progenitor RBC", "T", "T",
               "T", "Mono", "B",
               "CTL", "T", "CTL", "NK",
               "GMP", 
               "CTL", "Mono", "pDC", 
               "CTL", "MAIT", "B", "cDC",
               "NK", "ProB", "megakaryocyte progenitor cell",
               "CTL", "Plasma", "HSC",
               "early lymphoid progenitor", "ProDC",
               "ProB")

cell_typel2 <- as.vector(unique(bmcite_pred$celltype.l2))
match_cell_name <- function(celltypel2, c_cell_typel2, cell_names){
  idx <- match (celltypel2, c_cell_typel2)
  cell_type <- cell_names[[idx]]
  
  return(cell_type) 
}

cell_types_classifier <- c()
for (i in bmcite_pred$celltype.l2){
  cell_types_classifier <- append(cell_types_classifier, match_cell_name(i, cell_typel2, classifier_name))
}

bmcite_pred$cell_type.classifier <- cell_types_classifier
write.table(bmcite_pred, file=paste0("bmcite", ".nonbi.predict_2.txt"), sep = "\t")

cell_typel2
confusionMatrix(bmcite_pred$RF2.class, bmcite_pred$cell_type.classifier)





#significance 
marker_importance <- as.data.frame(rf.all.inner$importance)
sort_mark <- marker_importance[order(marker_importance$MeanDecreaseGini),]
barplot(height = des$MeanDecreaseGini[1:50],
        names = rownames(des)[1:50],
        )

rownames(rf.all.inner$importance)[1:100]
des <- marker_importance %>% arrange(desc(MeanDecreaseGini))

desc(marker_importance, MeanDecreaseGini )

des$name <- rownames(des)
ggplot(des[1:50, ], aes(x= reorder(name, -MeanDecreaseGini), y=MeanDecreaseGini)) +
  geom_bar(stat = "identity") +
  coord_flip()
