# read loom file with R
library(SeuratDisk)
library(Seurat)
library(anndata)
library(rhdf5)
library(Signac)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
cbmc <- readRDS('/home/chengyue/data/multi-omics/GSE214979_revised_seurat.rds')

cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.1, verbose = FALSE)
Idents(cbmc) <- 'wsnn_res.0.1'
cbmc <- RunUMAP(cbmc, dims = 1:10)
cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3"
)

p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.1', label = TRUE, cols = cluster_colors)


DefaultAssay(cbmc)<-'RNA'
my_color <- c("#CE3D32FF", "#5DB1DDFF", "#749B58FF", "#E4AF69FF")
VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)#, cols = my_color)
# # dev.off()
my_color1 <- c("#CE3D32FF", "#5DB1DDFF")
VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE, cols = my_color1)

allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)
Idents(cbmc) <- 'wsnn_res.0.2'
cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f",
  "#a1d99b"
)
DimPlot(seurat_object, reduction = "umap", cols = cluster_colors)
cbmc <- RunUMAP(cbmc, dims = 1:10)
p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.2', label = TRUE, cols = cluster_colors)


DefaultAssay(cbmc)<-'RNA'
Idents(cbmc) <- 'wsnn_res.0.2'
# my_color <- c("#CE3D32FF", "#5DB1DDFF", "#749B58FF", "#E4AF69FF")
VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)#, cols = my_color)

# my_color1 <- c("#CE3D32FF", "#5DB1DDFF")

VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE)#, cols = my_color1)

allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)

DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:10)
cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f",
  "#a1d99b", "#9ecae1", "#fdae6b", "#c994c7", "#d9d9d9",
  "#bcbd22", "#17becf", "#c7e9b4"
)

p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.3', label = TRUE, cols = cluster_colors)

Idents(cbmc) <- 'wsnn_res.0.3'
DefaultAssay(cbmc)<-'RNA'
# my_color <- c("#CE3D32FF", "#5DB1DDFF", "#749B58FF", "#E4AF69FF")
VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)#, cols = my_color)

my_color1 <- c("#CE3D32FF", "#5DB1DDFF")
VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE)#, cols = my_color1)
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])

pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)



cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.4, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:10)

cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f",
  "#a1d99b", "#9ecae1", "#fdae6b", "#c994c7", "#d9d9d9",
  "#bcbd22", "#17becf", "#c7e9b4", "#fdcdac", "#b2abd2"
)

p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.4', label = TRUE, cols = cluster_colors)

DefaultAssay(cbmc)<-'RNA'
Idents(cbmc) <- 'wsnn_res.0.4'
VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)
VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE)#, cols = my_color1)

allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)




cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:10)
cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f",
  "#a1d99b", "#9ecae1", "#fdae6b", "#c994c7", "#d9d9d9",
  "#bcbd22", "#17becf", "#c7e9b4", "#fdcdac", "#b2abd2",
  "#fbb4ae", "#decbe4"
)

p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.5', label = TRUE, cols = cluster_colors)
DefaultAssay(cbmc)<-'RNA'
Idents(cbmc) <- 'wsnn_res.0.5'

VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)
VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE)#, cols = my_color1)

allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.6, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:10)
cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f",
  "#a1d99b", "#9ecae1", "#fdae6b", "#c994c7", "#d9d9d9",
  "#bcbd22", "#17becf", "#c7e9b4", "#fdcdac", "#b2abd2",
  "#fbb4ae", "#decbe4", "#ccebc5", "#fed9a6", "#ffffcc"
)

p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.6', label = TRUE, cols = cluster_colors)
ggsave(p, '/home/chengyue/data/multi-omics/GSE214979_seurat_0.6_col.pdf', width = 8, height = 6)
DefaultAssay(cbmc)<-'RNA'
Idents(cbmc) <- 'wsnn_res.0.6'
my_color <- c("#CE3D32FF", "#5DB1DDFF", "#749B58FF", "#E4AF69FF")
VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)#, cols = my_color)
my_color1 <- c("#CE3D32FF", "#5DB1DDFF")
VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE)#, cols = my_color1)
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)



cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.7, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:10)
cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f",
  "#a1d99b", "#9ecae1", "#fdae6b", "#c994c7", "#d9d9d9",
  "#bcbd22", "#17becf", "#c7e9b4", "#fdcdac", "#b2abd2",
  "#fbb4ae", "#decbe4", "#ccebc5", "#fed9a6", "#ffffcc",
  "#8c510a"
)
p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.7', label = TRUE, cols = cluster_colors)
DefaultAssay(cbmc)<-'RNA'
Idents(cbmc) <- 'wsnn_res.0.7'
my_color <- c("#CE3D32FF", "#5DB1DDFF", "#749B58FF", "#E4AF69FF")
VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)#, cols = my_color)
my_color1 <- c("#CE3D32FF", "#5DB1DDFF")
VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE)#, cols = my_color1)
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)
allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:10)
cluster_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#cab2d6", "#6a3d9a",
  "#ff7f00", "#b15928", "#fdbf6f", "#b3de69", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3b3b3", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f",
  "#a1d99b", "#9ecae1", "#fdae6b", "#c994c7", "#d9d9d9",
  "#bcbd22", "#17becf", "#c7e9b4", "#fdcdac", "#b2abd2",
  "#fbb4ae", "#decbe4", "#ccebc5", "#fed9a6", "#ffffcc",
  "#8c510a", "#01665e", "#542788", "#bf812d", "#35978f",  
  "#dfc27d"
)

p<-DimPlot(cbmc, reduction = "umap", group.by = 'wsnn_res.0.8', label = TRUE, cols = cluster_colors)
DefaultAssay(cbmc)<-'RNA'
Idents(cbmc) <- 'wsnn_res.0.8'
my_color <- c("#CE3D32FF", "#5DB1DDFF", "#749B58FF", "#E4AF69FF")
VlnPlot(cbmc, features = c("DOCK8", "ARHGAP24", "ST6GAL1", "ARHGAP15"), stack = TRUE)#, cols = my_color)
my_color1 <- c("#CE3D32FF", "#5DB1DDFF")
VlnPlot(cbmc, features = c("ATP10A", "ABCB1"), stack = TRUE)#, cols = my_color1)
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Genes", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)


DefaultAssay(cbmc)<-'ATAC'
allmarker <- FindAllMarkers(cbmc , only.pos = TRUE,logfc.threshold = 0.25, min.pct = 0.25)

allmarker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

allmarker %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) -> top6
top_gene<-top6$gene
top_gene<-unique(top_gene)
cluster.averages <- AverageExpression(cbmc , 
                                      assays = "ATAC", 
                                      return.seurat = TRUE,
                                      slot = "data")
mat <- GetAssayData(cluster.averages, slot = "data")
mat <- as.matrix(mat[top_gene, ])
pheatmap(mat, scale = "row", cellwidth = 14, cellheight = 10, main="Peaks", cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE)#, color = colors)
