options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(randomForest)


## load classifier
load("RF1.nonbi.RData", verbose = TRUE)   # as generated using the RFtrain script
load("RF2.nonbi.RData", verbose = TRUE)


## load data
# f <- commandArgs(trailingOnly=TRUE)[1]  # input RData file (as generated using the prepareBackSpinKNN script)


X_test <- X_test[genes, ]
X_test <- X
setdiff(rownames(rf.all.inner$importance), rownames(X_test))
# X1<- readRDS(paste(working_dir, "/SCT_BM_expression_mat.rds", sep = ""))
# X1 <- X1@assays[["SCT"]]@data
# intersect(rownames(X), intersect(rownames(X1), rownames(Y)))
## predict, within seconds
# RF1
X.predict.all <- predict(rf.all.inner, newdata = t(X_test[rownames(rf.all.inner$importance), ]), type="prob")
X.predict.all.max <- apply(X.predict.all, 1, max)
X.predict.all.class <- colnames(X.predict.all)[apply(X.predict.all, 1, which.max)]
table(X.predict.all.class)

# RF2
X.predict.tumor <- predict(rf.tumor.inner, newdata = t(X_test[rownames(rf.tumor.inner$importance), ]), type="prob")
X.predict.tumor.max <- apply(X.predict.tumor, 1, max)
X.predict.tumor.class <- colnames(X.predict.tumor)[apply(X.predict.tumor, 1, which.max)]
table(X.predict.tumor.class)


## write output
Sobj_Xtest <- readRDS("C:/Users/jtao/work_dir/AML/data/bmcite6000.SCT.rds")
D.stats<- Sobj_Xtest@meta.data
D.stats<-a

d.out <- data.frame(D.stats[names(X.predict.tumor.max), ],
                    RF1.class=X.predict.all.class, RF1.score=X.predict.all.max, RF1=X.predict.all, 
                    RF2.class=X.predict.tumor.class, RF2.score=X.predict.tumor.max, RF2=X.predict.tumor)
write.table(d.out, file=paste0("bmcite.nonbi", ".predict.txt"), quote = FALSE, sep = "\t")

save(list=c("X.predict.all", "X.predict.tumor"), file=paste0("bmcite.nonbi", ".predict.RData"))

dev.off()

p1 <- hist(d.out$RF2.score, main = "prediction score distribution on training set")
p1
