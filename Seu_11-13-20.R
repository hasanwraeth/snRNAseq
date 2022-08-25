library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)


#Simple 2700 PBMCs -------------------
#loading data
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#QC 
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

#Feature Selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge=0, ynudge=0)
plot1+ plot2

#Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "./pbmc_tutorial.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                               "FCGR3A", "LYZ", "PPBP", "CD8A"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T",
                     "FCGR3A+ Mono","NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "./pbmc3k_final.rds")

#Multiple Dataset Integration----------------------
#Data Preprocessing
InstallData("panc8")

data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)}

#Integration of 3 datasets
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

#Cell type classification using an integrated reference
pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype, 
                            dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

table(pancreas.query$predicted.id)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

#Mouse Cell Atlas ---------------------
mca.matrix <- readRDS(file = "./MCA/MCA_merged_mat.rds")
mca.metadata <- read.csv(file = "./MCA/MCA_All-batch-removed-assignments.csv", row.names = 1)

mca <- CreateSeuratObject(counts = mca.matrix, meta.data = mca.metadata, project = "MouseCellAtlas")
# Only keep annotated cells
mca <- subset(mca, cells = names(which(!is.na(mca$ClusterID))))
# Leaves us with 242k cells
mca

#Data Preprocessing
mca <- NormalizeData(mca, normalization.method = "LogNormalize", scale.factor = 10000)
mca <- FindVariableFeatures(mca)
mca[["percent.mt"]] <- PercentageFeatureSet(mca, pattern = "^mt-")
mca <- ScaleData(mca, vars.to.regress = "percent.mt")

#Dimensional Reduction (PCA)
mca <- RunPCA(mca, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
ElbowPlot(mca, ndims = 100)
DimHeatmap(mca, dims = c(1:3, 70:75), cells = 500, balanced = TRUE)

#Graph-based clustering
mca <- FindNeighbors(mca, reduction = "pca", dims = 1:75, nn.eps = 0.5)
mca <- FindClusters(mca, resolution = 3, n.start = 10)

#isntall fftw and fast_tsne before running
#Visualization (t-SNE/UMAP)
mca <- RunTSNE(mca, dims = 1:75, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
mca <- RunUMAP(mca, dims = 1:75, min.dist = 0.75)

library(ggplot2)
p1 <- DimPlot(mca, reduction = "tsne", pt.size = 0.1) + ggtitle(label = "FIt-SNE")
p2 <- DimPlot(mca, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP")
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
(p1 + p2) & NoLegend()

p3 <- FeaturePlot(mca, features = c("S100a9", "Sftpc"), reduction = "umap", pt.size = 0.1, combine = FALSE)
p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
wrap_plots(p3)









