# read loom file with R
library(SeuratDisk)
library(Seurat)
library(rhdf5)
library(Signac)
library(Matrix)
library(EnsDb.Hsapiens.v86)
source("r_utils.R")

matrix_rna = './GSE214979/paired_RNA/matrix_pairedrna.mtx'
matrix_atac = './GSE214979/paired_ATAC/matrix_pairedatac.mtx'
feature_rna = './GSE214979/paired_RNA/gene.tsv'
feature_atac = './GSE214979/paired_ATAC/peak.tsv'
barcode = './GSE214979/paired_RNA/barcodes.csv'
frag_path = './GSE214979/fragments.tsv.gz'

start <- proc.time()

cbmc = read_rna(matrix_rna, feature_rna, barcode)
atac = read_atac(matrix_atac, feature_atac, barcode, frag_path)

cbmc[["ATAC"]] <- atac
DefaultAssay(cbmc) <- "RNA"

cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)

DefaultAssay(cbmc) <- "ATAC"
VariableFeatures(cbmc) <- rownames(cbmc[["ATAC"]])
cbmc <- NormalizeData(cbmc)
cbmc <- ScaleData(cbmc) 
cbmc <- RunPCA(cbmc, reduction.name = 'apca')
cbmc<- FindNeighbors(cbmc, dims = 1:20)
cbmc <- FindMultiModalNeighbors(
  cbmc, reduction.list = list("pca", "apca"), #, "apca"
  dims.list = list(1:30, 1:9), modality.weight.name = "RNA.weight"#, 1:16
)
cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 0.4, verbose = FALSE)
# seurat cluster based on louvain/leiden amd resolution effect cluster number
cbmc <- RunUMAP(cbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pred <- as.vector(cbmc$seurat_clusters)

write.csv(pred, './GSE214979/seurat-output/pred.csv')
write.csv(embedding, './GSE214979/seurat-output/z.csv', row.names = FALSE)

end <- proc.time()
elapsed <- end - start
print(elapsed)