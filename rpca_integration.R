
library(Seurat)

#data integration 
DefaultAssay(X_test)<-"SCT"
DefaultAssay(X) <- "SCT"
Seurat_list <- list(X, X_test)
features <- SelectIntegrationFeatures(object.list = Seurat_list)

Seurat_list <- lapply(X = Seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


immune.anchors <- FindIntegrationAnchors(object.list = Seurat_list, anchor.features = features, reduction = "rpca")
combined <- IntegrateData(anchorset = immune.anchors)


#data integration with SCT

Seurat_list <- lapply(X = Seurat_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = Seurat_list, nfeatures = 6000)
Seurat_list <- PrepSCTIntegration(object.list = Seurat_list, anchor.features = features)
Seurat_list <- lapply(X = Seurat_list, FUN = RunPCA, features = features)
Seurat.anchors <- FindIntegrationAnchors(object.list = Seurat_list, normalization.method = "SCT", 
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 30)
Seurat.combined.sct <- IntegrateData(anchorset = Seurat.anchors, normalization.method = "SCT", dims = 1:30)
Seurat.combined.sct <- RunPCA(Seurat.combined.sct, verbose = FALSE)
Seurat.combined.sct <- RunUMAP(Seurat.combined.sct, reduction = "pca", dims = 1:30)