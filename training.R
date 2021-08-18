# load("BM_6915cell.RData", verbose = TRUE)  # as generated using the prepareBackSpinKNN script
# install.packages("randomForest")
library("randomForest")


project<-"AML-classifier"
Sys.setenv(language="en")
if (project == "AML-classifier"){
  working_dir<-"C:/Users/jtao/work_dir/AML/data/AMLPaper"
  setwd(working_dir)
  graph_dir<- "C:/Users/jtao/work_dir/AML/data/images/"
  if (dir.exists(graph_dir)){
    graph_dir <- graph_dir
  } else {dir.create(graph_dir)}
}


# load BackSPIN cell type annotation
a <- read.table(paste(working_dir, "/BM_6915cells.BackSPIN.txt", sep = ""), row.names=1, header=TRUE)  # cell annotations available on GEO
X<- readRDS(paste(working_dir, "/SCT_BM_expression_mat.rds", sep = ""))

X<- X@assays[["SCT"]]@data
Y<- readRDS(paste(working_dir, "/SCT_AMLmut_expression_mat.rds", sep = ""))
Y<- Y@assays[["SCT"]]@data

X_test <- readRDS("C:/Users/jtao/work_dir/AML/data/bmcite6000.SCT.rds")
X_test <- X_test@assays[["SCT"]]@data
genes <- as.vector(intersect(rownames(X_test),
                             intersect(rownames(Y), rownames(X))))

#vectorize the dataset

# X<- as.matrix((X>0) + 0)
# X_test<- as.matrix((X_test>0)+0)
# Y<- as.matrix((Y>0) + 0)
# X <- log(X+1)  # filtered gene set
# X<-t(X)

X <- as.matrix(X[genes, ])
Y <- as.matrix(Y[genes, ])
X_test <- as.matrix(X_test[genes,])
# table(c(X))
# setdiff(rownames(X), rownames(Y))
# setdiff(colnames(X), rownames(a))

unique(a$cluster)
sample<- intersect(rownames(a), colnames(X))
a_sample <- a[rownames(a) %in% sample, ]

a_sample <- data.frame(a_sample)
rownames(a_sample) <- sample
colnames(a_sample) <- "cluster"
### 2) Train RF classifier #1 (normal BM data only)
pdf("RF1-bina_replaced.pdf", width = 6, height = 6)   # updated for cv heatmap
par(mar=c(4, 4, 4, 4))

## outer classifier (variable selection)
rf.all.expr <- names(which(rowMeans(exp(X)-1) > 0.01))
length(rf.all.expr)
set.seed(123)
rf.all.outer <- randomForest(x = t(X[rf.all.expr, rownames(a_sample)]),
                             y = factor(a_sample$cluster, unique(a_sample$cluster)),
                             sampsize = rep(300, length(unique(a_sample$cluster))),
                             replace = TRUE,
                             ntree = 1000,
                             do.trace=TRUE)

plot(sort(rf.all.outer$importance[, 1], decreasing = TRUE), type="l")
abline(v=1000)

cols <- apply(colorRamp(c("#FFFFFF", "red"))(seq(0, 1, 0.001)), 1, function(rr) rgb(rr[1], rr[2], rr[3], maxColorValue = 255))
image(as.matrix(seq(0, 1, 0.001)), col=cols, main="scale")

z <- rf.all.outer$confusion[, rownames(rf.all.outer$confusion)] / rowSums(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM, outer classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)]==0, NA, rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])))

1 - sum(diag(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])) / sum(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])

# select 1000 most informative genes (genes chosen most often in the outer classifier)
rf.all.outer.used1k <- names(sort(table(rownames(rf.all.outer$importance)[rf.all.outer$forest$bestvar[rf.all.outer$forest$bestvar!=0]]), decreasing = TRUE)[1:1000])
plot(rf.all.outer$importance[rf.all.outer.used1k, 1])  # correlates to importance measure


## inner classifier
set.seed(123)
rf.all.inner <- randomForest(x = t(X[rf.all.outer.used1k, rownames(a_sample)]),
                             y = factor(a_sample$cluster, unique(a_sample$cluster)),
                             sampsize = rep(30, length(unique(a_sample$cluster))),
                             ntree = 1000,
                             do.trace=TRUE)

z <- rf.all.inner$confusion[, rownames(rf.all.inner$confusion)] / rowSums(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM, inner classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)]==0, NA, rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])))

1 - sum(diag(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])) / sum(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])

## inner classifier, 5-fold cross validation
y <- factor(a_sample$cluster, unique(a_sample$cluster))
names(y) <- rownames(a_sample)
x <- rownames(a_sample)
cv <- split(x, rep(1:5, 1E6)[1:length(y)])
lengths(cv)

rf.all.cv <- lapply(cv, function(y2) {
  set.seed(123)
  randomForest(x = t(X[rf.all.outer.used1k, setdiff(x, y2)]),
               y = y[setdiff(x, y2)],
               sampsize = rep(10, length(unique(y))),
               ntree = 1000,
               do.trace=TRUE)
})

# setcv[1]
#predict the sets that were not used for training
rf.all.cv.predict.prob <- lapply(rf.all.cv, function(rf) {
  predict(rf, t(X[rownames(rf$importance), setdiff(colnames(X), names(rf$y))]), type = "prob")
})
rf.all.cv.predict <- lapply(rf.all.cv.predict.prob, function(x) {
  y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))  # max
  names(y) <- rownames(x)
  y
})

# generate CV confusion matrix
rf.all.cv.predict.matrix <- table(unlist(rf.all.cv.predict), y[unlist(lapply(rf.all.cv.predict, names))])
rf.all.cv.predict.matrix <- rf.all.cv.predict.matrix[, rownames(rf.all.cv.predict.matrix)]
sum(rf.all.cv.predict.matrix)
write.table(rf.all.cv.predict.matrix, sep="\t", quote=FALSE, file="RF1.cv.txt")

z <- t(t(rf.all.cv.predict.matrix) / colSums(rf.all.cv.predict.matrix))

heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM, inner classifier, CV confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.all.cv.predict.matrix==0, NA, rf.all.cv.predict.matrix)))

1 - sum(diag(rf.all.cv.predict.matrix)) / sum(rf.all.cv.predict.matrix)
heatmap
# save_png(heatmap, "BM, inner classifier CV", 10, 10)
dev.off()

rf.all.cv.predict.matrix
## export RData object for RF1
save(list=c("rf.all.outer", "rf.all.inner", "rf.all.cv"), file="RF1_nonbi.RData")


### 3) Classify cells from tumor samples with detected mutations
# load tumor data (cell mutation status available on GEO)
# load("AMLmut_923cells.RData", verbose = TRUE)    # as generated using the prepareBackSpinKNN script


## classify
pdf("RF1.AMLclassify.nonbinarized.pdf", width = 6, height = 6)
par(mar=c(4, 4, 4, 4))
Y[rownames(rf.all.inner$importance), ]

setdiff(rownames(rf.all.inner$importance), rownames(Y))
rf.all.inner.mut <- predict(rf.all.inner, newdata = t(Y[rownames(rf.all.inner$importance), ]), type = "prob")
rf.all.inner.mut.max <- factor(colnames(rf.all.inner.mut)[apply(rf.all.inner.mut, 1, which.max)], colnames(rf.all.inner.mut))
names(rf.all.inner.mut.max) <- rownames(rf.all.inner.mut)
barplot(table(rf.all.inner.mut.max), las=2)

# reclassify Prog with highest score for HSC to HSC to have 65 cells in total
rf.all.inner.mut.prog2hsc <- names(sort(rf.all.inner.mut[rf.all.inner.mut.max == "Prog", "HSC"], decreasing = TRUE)[1:(65-table(rf.all.inner.mut.max)["HSC"])])
plot(rf.all.inner.mut[, "Prog"], rf.all.inner.mut[, "HSC"], col=ifelse(rownames(rf.all.inner.mut) %in% rf.all.inner.mut.prog2hsc, "red", "black"), xlim=c(0, 0.5), ylim=c(0, 0.5), pch=16)
abline(0, 1)

# reclassify earlyEry with higher score for Prog than for lateEry to Prog
rf.all.inner.mut.ery2prog <- names(which(rf.all.inner.mut[rf.all.inner.mut.max=="earlyEry", "Prog"] > rf.all.inner.mut[rf.all.inner.mut.max=="earlyEry", "lateEry"]))
plot(rf.all.inner.mut[, "Prog"], rf.all.inner.mut[, "earlyEry"], xlim=c(0, 0.5), ylim=c(0, 0.5), col=ifelse(rownames(rf.all.inner.mut) %in% rf.all.inner.mut.ery2prog, "red", "black"), pch=16)
plot(rf.all.inner.mut[, "Prog"], rf.all.inner.mut[, "lateEry"], xlim=c(0, 0.5), ylim=c(0, 0.5), col=ifelse(rownames(rf.all.inner.mut) %in% rf.all.inner.mut.ery2prog, "red", "black"), pch=16)

# create list
a.mut <- split(names(rf.all.inner.mut.max), rf.all.inner.mut.max)[c("HSC", "Prog", "GMP", "ProMono", "Mono", "cDC")]
lengths(a.mut)
a.mut$HSC <- c(a.mut$HSC, rf.all.inner.mut.prog2hsc)
a.mut$Prog <- setdiff(a.mut$Prog, rf.all.inner.mut.prog2hsc)
a.mut$Prog <- c(a.mut$Prog, rf.all.inner.mut.ery2prog)
lengths(a.mut)

dev.off()


### 4) Train RF classifier #2 (normal BM and AML cells with mutations)
pdf("RF2_nonbinarized.pdf", width = 6, height = 6)
par(mar=c(4, 4, 4, 4))

a2 <- c(split(rownames(a_sample), a_sample$cluster)[unique(a_sample$cluster)], Tumor=a.mut)
lengths(a2)

XY <- cbind(X, Y[rownames(X), ])[, unlist(a2)]
dim(XY)

## outer classifier (variable selection)
rf.tumor.expr <- names(which(rowMeans(XY[, unlist(a2)]) >= 0.01))
length(rf.tumor.expr)

set.seed(123)
rf.tumor.outer <- randomForest(x = t(XY[rf.tumor.expr, unlist(a2)]),
                               y = factor(rep(names(a2), lengths(a2)), names(a2)),
                               sampsize = rep(30, length(a2)),
                               ntree = 1000,
                               do.trace=TRUE)

z <- rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)] / rowSums(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM+mut, outer classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)]==0, NA, rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])))

1 - sum(diag(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])) / sum(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])

rf.tumor.outer.used1k <- names(sort(table(rownames(rf.tumor.outer$importance)[rf.tumor.outer$forest$bestvar[rf.tumor.outer$forest$bestvar!=0]]), decreasing = TRUE)[1:1000])
plot(rf.tumor.outer$importance[rf.tumor.outer.used1k, 1])

## inner classifer
set.seed(123)
rf.tumor.inner <- randomForest(x = t(XY[rf.tumor.outer.used1k, unlist(a2)]),
                               y = factor(rep(names(a2), lengths(a2)), names(a2)),
                               sampsize = rep(30, length(a2)),
                               ntree = 1000,
                               do.trace=TRUE)

z <- rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)] / rowSums(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM+mut, inner classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)]==0, NA, rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])))

1 - sum(diag(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])) / sum(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])

## inner, 5-fold cv
y <- factor(rep(names(a2), lengths(a2)), names(a2))
names(y) <- unlist(a2)
x <- unlist(a2)
cv <- split(x, rep(1:5, 1E6)[1:length(y)])

rf.tumor.cv <- lapply(cv, function(y2) {
  set.seed(123)
  randomForest(x = t(XY[rf.tumor.outer.used1k, setdiff(x, y2)]),
               y = y[setdiff(x, y2)],
               sampsize = rep(20, length(unique(y))),
               ntree = 1000,
               do.trace=TRUE)
})

rf.tumor.cv.predict.prob <- lapply(rf.tumor.cv, function(rf) {
  predict(rf, t(XY[rownames(rf$importance), setdiff(colnames(XY), names(rf$y))]), type = "prob")
})
rf.tumor.cv.predict <- lapply(rf.tumor.cv.predict.prob, function(x) {
  y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))
  names(y) <- rownames(x)
  y
})

rf.tumor.cv.predict.matrix <- table(unlist(rf.tumor.cv.predict), y[unlist(lapply(rf.tumor.cv.predict, names))])
rf.tumor.cv.predict.matrix <- rf.tumor.cv.predict.matrix[, rownames(rf.tumor.cv.predict.matrix)]
sum(rf.tumor.cv.predict.matrix)
write.table(rf.all.cv.predict.matrix, sep="\t", quote=FALSE, file="RF2.cv.nonbi.txt")

z <- t(t(rf.tumor.cv.predict.matrix) / colSums(rf.tumor.cv.predict.matrix))
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM+mut, inner classifier, CV confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.tumor.cv.predict.matrix==0, NA, rf.tumor.cv.predict.matrix)))

1 - sum(diag(rf.tumor.cv.predict.matrix)) / sum(rf.tumor.cv.predict.matrix)

dev.off()



## export RData object for RF2
save(list=c("rf.tumor.outer", "rf.tumor.inner", "rf.tumor.cv"), file="RF2.nonbi.RData")


