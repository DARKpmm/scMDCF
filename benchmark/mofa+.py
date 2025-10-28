import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import torch
print(torch.version.cuda)
import scanpy as sc
import muon as mu
import argparse
from muon import MuData
import numpy as np
from sklearn.cluster import KMeans
from utils import read_data, eva, read_data_nolabel, GetCluster, eva_nolabel
from sklearn.metrics import silhouette_score
from time import time
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
import resource

t0=time()
def parameter_setting():
    
    parser = argparse.ArgumentParser(description='train')
    parser.add_argument('--file_path1', default='./GSE214979_rna.h5ad')
    parser.add_argument('--file_path2', default='./GSE214979_atac.h5ad')
    parser.add_argument('--label_file', default='./id.csv')
    parser.add_argument('--file_type', default='h5ad', type=str)
    return parser

parser=parameter_setting()
args = parser.parse_args()
adata_RNA, adata_ATAC, cluster_number, y = read_data(args.file_path1, args.file_path2, args.file_type, args.label_file)

mdata = MuData({'rna': adata_RNA, 'atac': adata_ATAC})
print(mdata)
mu.tl.mofa(mdata)
print(mdata)
latent = mdata.obsm['X_mofa']
np.savetxt('./GSE214979/mofa-output/z.txt', latent)
kmeans = KMeans(n_clusters = cluster_number, n_init=20)
y_pred = kmeans.fit_predict(latent)   
np.savetxt('./GSE214979/mofa-output/ypred.txt', y_pred)
nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred)
print('z for clustering, NMI:{:.4f}, ARI:{:.4f}, AMI:{:.4f}, FMI:{:.4f}, HOM:{:.4f}, COM:{:.4f}, V:{:.4f}'.format(nmi, ari, ami, fmi, hom, com, v))

db, ch, asw = eva_nolabel(latent, y_pred)
print('z for clustering, DB:{:.4f}, CH:{:.4f}, ASW:{:.4f}'.format(db, ch, asw))  

t1=time()
t=t1-t0
print('time: ', t)
max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(f"max_mem: , {max_mem/1024/1024:.2f}GB")