import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import argparse
import torch

from utils import read_data, normalize
from layer import scMDCF
from train import pre_train, alt_train
from time import time
import warnings
warnings.filterwarnings('ignore')


def parameter_setting(): 
    
    parser = argparse.ArgumentParser(description='train')
    parser.add_argument('--file_name1', default='pbmc_10x_rna_public.h5ad')#peripheral_blood_rna.h5ad  pbmc_spector.h5 Pbmc10k-RNA
    parser.add_argument('--file_name2', default='pbmc_10x_atac_public.h5ad')
    parser.add_argument('--label_file', default=None)#'/home/chengyue/data/multi-omics/GSE128639_BMNC_ypred.txt')#'/home/chengyue/data/multi-omics/peripheral_blood_label.tsv'
    parser.add_argument('--save_results', default='False', type=bool)
    parser.add_argument('--file_type', default='h5ad', type=str)
    parser.add_argument('--model_file', default='/home/chengyue/data/multi-omics/test1.pth.tar')
    parser.add_argument("--highly_genes", default = 2500, type = int)#SNARE-seq 2500; CITE-seq 1000
    parser.add_argument("--lr_pre", default = "1e-2", type = float)
    parser.add_argument("--lr_alt", default = "1e-3", type = float)
    parser.add_argument("--epoch_pre", default = "200", type = int)
    parser.add_argument("--epoch_alt", default = "200", type = int)
    parser.add_argument("--device", default='cuda:2', type=str)
    parser.add_argument("--enc1", default = "512", type = int)
    parser.add_argument("--enc2", default = "64", type = int)
    parser.add_argument("--zdim", default = "8", type = int)
    parser.add_argument("--alpha", default = "1.", type = float)
    parser.add_argument("--gamma", default = "1.", type = float)
    parser.add_argument("--lamb", default = "10", type = float)# 10 for atac; 0.5 for adt
    
    return parser

parser=parameter_setting()
args = parser.parse_args()
file_path1 = '/home/chengyue/scMDCF/dataset/'+args.file_name1
file_path2 = '/home/chengyue/scMDCF/dataset/'+args.file_name2
adata_RNA, adata_ATAC, cluster_number, y = read_data(file_path1, file_path2, args.file_type, args.label_file)
adata_RNA = normalize(adata_RNA, highly_genes=args.highly_genes, normalize_input=True)
adata_ATAC = normalize(adata_ATAC, highly_genes=args.highly_genes, normalize_input=True)

print(adata_RNA)
print(adata_ATAC)
args.RNA_input = adata_RNA.X.shape[1]
args.ATAC_input = adata_ATAC.X.shape[1]

args.n_cell = adata_RNA.X.shape[0]
args.n_clusters = cluster_number

args.layere_omics1_view = [adata_RNA.X.shape[1], args.enc1, args.enc2]#for atac adata_RNA.X.shape[1], args.enc1, args.enc2; for adt adata_RNA.X.shape[1], args.enc1, args.enc2, args.zdim
args.layere_omics2_view = [adata_ATAC.X.shape[1], args.enc1, args.enc2]#for atac adata_ATAC.X.shape[1], args.enc1, args.enc2; adata_ATAC.X.shape[1], args.enc1, args.enc2, args.zdim
args.layerd_omics1_view = [args.zdim, args.enc2, args.enc1, adata_RNA.X.shape[1]]
args.layerd_omics2_view = [args.zdim, args.enc2, args.enc1, adata_ATAC.X.shape[1]]
args.fusion_layer = [2*args.enc2, 32, args.zdim]

model = scMDCF(args).to(args.device)
print(model)
x_rna, x_atac = torch.from_numpy(adata_RNA.X).to(args.device).float(), torch.from_numpy(adata_ATAC.X).to(args.device).float() #for human pbmc and brain 10xmalt and bmnc .float()
t0=time()
pre_train(args, model, x_rna, x_atac, y)

alt_train(args, model, x_rna, x_atac, y)
print('Total time:{:.4f} seconds. Cell number:{}'.format(time()-t0, args.n_cell))