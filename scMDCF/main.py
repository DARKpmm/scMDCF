import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import argparse
import torch
#from sklearn.metrics import silhouette_samples, silhouette_score

from utils import read_data, normalize#, get_adj, eva, read_dataset, CLR_normalization, GetCluster
from layer import scDEFR_multi
from train import pre_train, alt_train
from time import time
import warnings
warnings.filterwarnings('ignore')


def parameter_setting(): 
    #spleen_lymph GSE128639_BMNC inhouse_rna Pbmc10k-ATAC human_brain_3k_rna pbmc_10x_rna_public GSM4949911_tea_rna Multiome_nextgem_Chromium_rna
    parser = argparse.ArgumentParser(description='train')
    parser.add_argument('--file_path1', default='/home/chengyue/data/multi-omics/GSE128639_BMNC.h5')#peripheral_blood_rna.h5ad  pbmc_spector.h5 Pbmc10k-RNA
    parser.add_argument('--file_path2', default='/home/chengyue/data/multi-omics/human_pbmc_3k_atac.h5ad')
    parser.add_argument('--label_file', default='/home/chengyue/data/multi-omics/GSE128639_BMNC_ypred.txt')#'/home/chengyue/data/multi-omics/peripheral_blood_label.tsv'
    parser.add_argument('--save_results', default='False', type=bool)
    parser.add_argument('--file_type', default='h5', type=str)
    parser.add_argument('--model_file', default='/home/chengyue/data/multi-omics/test1.pth.tar')
    parser.add_argument("--highly_genes", default = 2500, type = int)#
    parser.add_argument("--lr_pre", default = "1e-2", type = float)#2
    parser.add_argument("--lr_alt", default = "1e-3", type = float)#3
    parser.add_argument("--epoch_pre", default = "200", type = int)#2
    parser.add_argument("--epoch_alt", default = "200", type = int)#2
    parser.add_argument("--beginkl", default = "200", type = int)
    parser.add_argument("--device", default='cuda:2', type=str)
    parser.add_argument("--enc1", default = "512", type = int)#512
    parser.add_argument("--enc2", default = "64", type = int)#64
    parser.add_argument("--zdim", default = "8", type = int)#
    parser.add_argument("--alpha", default = "1.", type = float)#
    parser.add_argument("--tau", default = "1", type = float)#
    parser.add_argument("--tol", default = "1e-3", type = float)#
    parser.add_argument("--batch_size", default = "256", type = int)
    parser.add_argument("--gamma", default = "1.", type = float)
    # parser.add_argument("--sigma1", default = "2.5", type = float)
    # parser.add_argument("--sigma2", default = "0.1", type = float)
    parser.add_argument("--lamb", default = "10", type = float)# 10 for atac; 0.5 for adt
    #alpha
    return parser
#import scanpy as sc
# data = sc.read_h5ad('/home/chengyue/data/multi-omics/spleen_lymph_206.h5ad')
# print(data.obsm['protein_expression'].shape[1])
# print(data.X.shape[1])

parser=parameter_setting()
args = parser.parse_args()
adata_RNA, adata_ATAC, cluster_number, y = read_data(args.file_path1, args.file_path2, args.file_type, args.label_file)
adata_RNA = normalize(adata_RNA, highly_genes=args.highly_genes, normalize_input=True)
adata_ATAC = normalize(adata_ATAC, highly_genes=args.highly_genes, normalize_input=True)
#adata_ATAC.X = CLR_normalization(adata_ATAC.X)#.toarray()
#adata_RNA.X = CLR_normalization(adata_RNA.X)
print(adata_ATAC)
print(adata_RNA)
args.RNA_input = adata_RNA.X.shape[1]
print(args.RNA_input)
print(y.shape[0])
args.ATAC_input = adata_ATAC.X.shape[1]
args.n_cell = adata_RNA.X.shape[0]
args.n_clusters = cluster_number
print(cluster_number)

args.layere_view = [adata_RNA.X.shape[1], args.enc1, args.enc2, args.zdim]#for atac adata_RNA.X.shape[1], args.enc1, args.enc2
args.layere_adt_view = [adata_ATAC.X.shape[1], args.enc1, args.enc2, args.zdim]#for atac adata_ATAC.X.shape[1], args.enc1, args.enc2
args.layerd_view = [args.zdim, args.enc2, args.enc1, adata_RNA.X.shape[1]]
args.layerd_adt_view = [args.zdim, args.enc2, args.enc1, adata_ATAC.X.shape[1]]

model = scDEFR_multi(args).to(args.device)
print(model)
x_rna, x_atac = torch.from_numpy(adata_RNA.X).to(args.device).float(), torch.from_numpy(adata_ATAC.X).to(args.device).float() #for human pbmc and brain 10xmalt and bmnc .float()
t0=time()
pre_train(args, model, x_rna, x_atac, y)

alt_train(args, model, x_rna, x_atac, y)
print('Total time:{:.4f} seconds. Cell number:{}'.format(time()-t0, args.n_cell))
# y_old = model.y_pred

# silhouette_avg = silhouette_score(model.z.data.cpu().numpy(), y_old)

# print(silhouette_avg)
