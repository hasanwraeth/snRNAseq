library(Seurat)

hepcom <- readRDS(file = "hep_harmony_annotated.rds")
p1 <- DimPlot(hepcom, reduction = "umap", group.by = "assay_type")
p2 <- DimPlot(hepcom, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2

o=subset(hepcom, subset = assay_type == "single_nuc")
p2 <- DimPlot(o, reduction = "umap", label = TRUE, repel = TRUE)
p2
p2 <- DimPlot(o, reduction = "pca", label = TRUE, repel = TRUE)
p2
p2 <- DimPlot(o, reduction = "harmony", label = TRUE, repel = TRUE)
p2
