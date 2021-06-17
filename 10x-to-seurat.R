# import libraries
library(Seurat)

### setup path for folder 
# path for cellranger outs folder
path = "/data/data_becavin/pig_multiome/data/outs/"
# path for dimensionality reduction files
reddim_path = paste0(path,"analysis/dimensionality_reduction/")
# path for clustering files
clust_path = paste0(path,"analysis/clustering/")

### setup path for datasets
# dimred datasets
h5_file = paste0(path,"filtered_feature_bc_matrix.h5")
gex_umap = paste0(reddim_path,"gex/umap_projection.csv")
gex_tsne = paste0(reddim_path,"gex/tsne_projection.csv")
gex_pca = paste0(reddim_path,"gex/pca_projection.csv")
atac_umap = paste0(reddim_path,"atac/umap_projection.csv")
atac_tsne = paste0(reddim_path,"atac/tsne_projection.csv")
atac_lsa = paste0(reddim_path,"atac/lsa_projection.csv")
# cluster datasets
gex_clust = paste0(clust_path,"gex/graphclust/clusters.csv")
atac_clust = paste0(clust_path,"atac/graphclust/clusters.csv")

seurat_file = paste0(path,"filtered_feature_bc_matrix.rds")

# set genome name
genome_name = "susScr11"

# Load h5 file
seurat_h5 <- Read10X_h5(h5_file)

# Create seuratobject
seurat_object <- CreateSeuratObject(counts = seurat_h5$`Gene Expression`)
seurat_object[['ATAC']] <- CreateChromatinAssay(
  counts = seurat_h5$`Peaks`,
  sep = c(":", "-"),
  genome = genome_name,
)

# Add dimred tables
# UMAP
dimred_table = read.table(gex_umap, sep=",", row.names = 1, header = T)
dimred <- CreateDimReducObject(as.matrix(dimred_table), key = "UMAP_", assay = "RNA")
seurat_object[["umap.rna"]] = dimred
dimred_table = read.table(atac_umap, sep=",", row.names = 1, header = T)
dimred <- CreateDimReducObject(as.matrix(dimred_table), key="UMAP_", assay = "ATAC")
seurat_object[["umap.atac"]] = dimred
# T-SNE
dimred_table = read.table(gex_tsne, sep=",", row.names = 1, header = T)
dimred <- CreateDimReducObject(as.matrix(dimred_table), key = "TSNE_", assay = "RNA")
seurat_object[["tsne.rna"]] = dimred
dimred_table = read.table(atac_tsne, sep=",", row.names = 1, header = T)
dimred <- CreateDimReducObject(as.matrix(dimred_table), key="TSNE_", assay = "ATAC")
seurat_object[["tsne.atac"]] = dimred
# PCA
dimred_table = read.table(gex_pca, sep=",", row.names = 1, header = T)
dimred <- CreateDimReducObject(as.matrix(dimred_table), key = "PCA_", assay = "RNA")
seurat_object[["pca.rna"]] = dimred
dimred_table = read.table(atac_lsa, sep=",", row.names = 1, header = T)
dimred <- CreateDimReducObject(as.matrix(dimred_table), key="PCA_", assay = "ATAC")
seurat_object[["lsa.atac"]] = dimred

# add graph clusters
clust_table = as.matrix(read.table(gex_clust, sep=",", row.names = 1, header = T))
seurat_object[["gex_graph_based"]] <- clust_table
Idents(seurat_object) <- "gex_graph_based"
clust_table = as.matrix(read.table(atac_clust, sep=",", row.names = 1, header = T))
seurat_object[["atac_graph_based"]] <- clust_table

# verify your umap for gex
DimPlot(seurat_object, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("UMAP RNA")
DimPlot(seurat_object, reduction = "tsne.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("t-SNE RNA")
Idents(seurat_object) <- "atac_graph_based"
DimPlot(seurat_object, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("UMAP ATAC")
DimPlot(seurat_object, reduction = "tsne.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("t-SNE ATAC")

# verify with feature plot
FeaturePlot(object = seurat_object, features = "nCount_RNA")
FeaturePlot(object = seurat_object, features = "nCount_ATAC")
FeaturePlot(object = seurat_object, features = "CFTR")

# normalize seurat data
seurat_object <- SCTransform(seurat_object)
DefaultAssay(seurat_object) <- 'ATAC'
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 10)
seurat_object <- RunTFIDF(seurat_object)

# verify with feature plot
DefaultAssay(seurat_object) <- "RNA"
FeaturePlot(object = seurat_object, features = "nCount_RNA")
FeaturePlot(object = seurat_object, features = "nCount_ATAC")
FeaturePlot(object = seurat_object, features = "CFTR")
DefaultAssay(seurat_object) <- "ATAC"
FeaturePlot(object = seurat_object, features = "18-28760563-28761524")
DefaultAssay(seurat_object) <- "RNA"


# save to rds
saveRDS(seurat_object, seurat_file)
