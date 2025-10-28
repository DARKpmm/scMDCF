import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'#https://blog.csdn.net/weixin_46713695/article/details/130134955 OpenBLAS warning
import mudata as md
import scvi
from scvi.model import TOTALVI
import argparse
from utils import read_data, eva, read_data_nolabel
from sklearn.cluster import KMeans
from time import time
scvi.settings.seed = 1000
#conda activate scvi-env
def parameter_setting():
    
    parser = argparse.ArgumentParser(description='train')
    parser.add_argument('--file_path1', default='./inhouse_rna.h5ad')
    parser.add_argument('--file_path2', default='./inhouse_adt.h5ad')
    parser.add_argument('--label_file', default=None)
    parser.add_argument('--file_type', default='h5ad', type=str)
    #parser.add_argument('--cluster_number', default='31', type=int)
    return parser

parser=parameter_setting()
args = parser.parse_args()
adata_RNA, adata_ATAC, cluster_number, y = read_data(args.file_path1, args.file_path2, args.file_type, args.label_file)

t0=time()
mdata = md.MuData({"rna": adata_RNA, "protein": adata_ATAC})
TOTALVI.setup_mudata(mdata,
                                rna_layer=None,
                                protein_layer=None,
                                modalities={"rna_layer": "rna",
                                            "protein_layer": "protein"})
vae = TOTALVI(mdata)
vae.train()
latent = vae.get_latent_representation()
kmeans = KMeans(n_clusters = cluster_number, n_init=20)
y_pred = kmeans.fit_predict(latent)   
from sklearn.metrics import silhouette_score
silhouette_avg = silhouette_score(latent, y_pred)
print(silhouette_avg)
import numpy as np
np.savetxt('./results/spleen_lymph_z_multivi.txt', latent)
np.savetxt('./results/spleen_lymph_ypred_multivi.txt', y_pred)
nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred) 
print('z for clustering, NMI:{:.4f}, ARI:{:.4f}, AMI:{:.4f}, FMI:{:.4f}, HOM:{:.4f}, COM:{:.4f}, V:{:.4f}'.format(nmi, ari, ami, fmi, hom, com, v))
t1=time()
t=t1-t0
print('time is:{:.4f}'.format(t))