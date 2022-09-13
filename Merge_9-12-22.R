library(Seurat)

#h5ad not possible see python instructions
#For h5 files
#pbmc8k.data = Read10X_h5("pbmc8k_raw_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
#pbmc8k=CreateSeuratObject(counts = pbmc8k.data, project = "PBMC8K")
#pbmc8k

pbmc4k.data = Read10X("filtered_gene_bc_matrices/GRCh38/")
pbmc4k=CreateSeuratObject(counts = pbmc4k.data, project = "PBMC4K")
pbmc4k


pbmc8k.data = Read10X("filtered_gene_bc_matrices 2/GRCh38/")
pbmc8k=CreateSeuratObject(counts = pbmc8k.data, project = "PBMC8K")
pbmc8k

pbmc.combined <- merge(pbmc4k, y = pbmc8k, add.cell.ids = c("4K", "8K"), project = "PBMC12K")
pbmc.combined

head(colnames(pbmc.combined))
table(pbmc.combined$orig.ident)

pbmc.list <- SplitObject(pbmc.combined, split.by = "orig.ident")
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = pbmc.list)
immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 1.2) 
#resolution 0.4 to 1.2 for 3k cells controls cluster no., higher for more cells


p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident")

cluster2.markers <- FindMarkers(immune.combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

cluster21.markers <- FindMarkers(immune.combined, ident.1 = 21, min.pct = 0.25)
head(cluster21.markers, n = 5)
