library(Seurat)
library(Signac)
library(Matrix)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
source("r_utils.R")
start <- proc.time()
matrix_rna = './GSE214979/paired_RNA/matrix_pairedrna.mtx'
matrix_atac = './GSE214979/paired_ATAC/matrix_pairedatac.mtx'
feature_rna = './GSE214979/paired_RNA/gene.tsv'
feature_atac = './GSE214979/paired_ATAC/peak.tsv'
barcode = './GSE214979/paired_RNA/barcodes.csv'
frag_path = './GSE214979/fragments.tsv.gz'
out_dir = './GSE214979/figr-output/'


paired_rna = read_rna(matrix_rna, feature_rna, barcode)
paired_atac = read_atac(matrix_atac, feature_atac, barcode, frag_path)

dataset_vec <- rep(c("Multiome-RNA","Multiome-ATAC"),
                   c(ncol(paired_rna),
                     ncol(paired_atac)))
paired_atac@meta.data$dataset <- "Multiome-ATAC"
paired_rna@meta.data$dataset <- "Multiome-RNA"

names(dataset_vec) <- c(paste0("prna_",colnames(paired_rna)),
                        paste0("patac_",colnames(paired_atac)))
print(table(dataset_vec))

paired_rna <- RenameCells(paired_rna,add.cell.id = "prna",for.merge = FALSE)
paired_atac <- RenameCells(paired_atac,add.cell.id = "patac",for.merge = FALSE)

DefaultAssay(paired_rna) <- "RNA"
paired_rna <- NormalizeData(paired_rna)
paired_rna <- FindVariableFeatures(paired_rna, nfeatures = 5000)


DefaultAssay(paired_atac) <- "ATAC"
paired_atac <- RunTFIDF(paired_atac)
paired_atac <- FindTopFeatures(paired_atac, min.cutoff = 50)
paired_atac <- RunSVD(paired_atac)
paired_atac <- RunUMAP(paired_atac, reduction = "lsi", dims = 1:15)

gene.activities <- GeneActivity(paired_atac)
paired_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(paired_atac) <- "ACTIVITY"
paired_atac <- NormalizeData(paired_atac)
paired_atac <- FindVariableFeatures(paired_atac, nfeatures = 5000)
paired_atac <- ScaleData(paired_atac)

gene.use <- intersect(VariableFeatures(paired_rna), 
                      VariableFeatures(paired_atac))

paired_rna[["RNA"]] <- as(object = paired_rna[["RNA"]], Class = "Assay")
paired_atac[['RNA']]<- as(object = paired_atac[["ACTIVITY"]], Class = "Assay")
cca_res <- RunCCA(object1 = paired_rna, 
                  object2 = paired_atac,
                  assay1 = "RNA",
                  assay2 = "RNA",
                  num.cc = 30,
                  features = gene.use,
                  renormalize = FALSE,
                  rescale = TRUE)
half_index <- nrow(cca_res@reductions[["cca"]]@cell.embeddings) %/% 2
embed <- cca_res@reductions[["cca"]]@cell.embeddings[1:half_index, ]

write.table(embed, "./GSE214979/figr-output/embedding.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

end <- proc.time()
elapsed <- end - start
print(elapsed)