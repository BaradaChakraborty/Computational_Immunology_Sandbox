library(Seurat)
library(dplyr)
library(ggplot2)

print("Initializing BioGrademy R Pipeline...")

counts <- read.csv("data/raw/GSE115978_counts.csv.gz", row.names = 1)
seu_obj <- CreateSeuratObject(counts = counts, project = "Melanoma_ICB", min.cells = 3, min.features = 200)

print("Executing Quality Control filtering...")
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")

seu_obj <- subset(seu_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

print("Applying LogNormalization...")
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)

print("Identifying highly variable features and scaling data...")
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj, features = rownames(seu_obj))

print("Running PCA, UMAP, and t-SNE algorithms...")
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))

seu_obj <- FindNeighbors(seu_obj, dims = 1:20)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)

seu_obj <- RunUMAP(seu_obj, dims = 1:20)
seu_obj <- RunTSNE(seu_obj, dims = 1:20)

dir.create("figures", showWarnings = FALSE)
p1 <- DimPlot(seu_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP: Tumor Microenvironment")
p2 <- DimPlot(seu_obj, reduction = "tsne", label = TRUE) + ggtitle("t-SNE: Tumor Microenvironment")
ggsave("figures/01_UMAP_Clustering.png", plot = p1, width = 7, height = 5)
ggsave("figures/02_tSNE_Clustering.png", plot = p2, width = 7, height = 5)

print("Mapping canonical T-Cell exhaustion markers...")
exhaustion_features <- c("CD8A", "PDCD1", "HAVCR2", "LAG3", "CTLA4")

p3 <- FeaturePlot(seu_obj, features = exhaustion_features, reduction = "umap", ncol = 3)
ggsave("figures/03_Exhaustion_Signatures.png", plot = p3, width = 12, height = 8)

dir.create("data/processed", showWarnings = FALSE)
saveRDS(seu_obj, file = "data/processed/seu_obj_processed.rds")
print("R Pipeline Complete. Object exported for Python ICB validation.")