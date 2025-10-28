import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import scanpy as sc
from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt
import os
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from time import time
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, adjusted_mutual_info_score
from sklearn.metrics import fowlkes_mallows_score, homogeneity_score, completeness_score, v_measure_score
t0=time()
def eva(y_true, y_pred):
    nmi = nmi_score(y_true, y_pred)
    ari = ari_score(y_true, y_pred)
    ami = adjusted_mutual_info_score(y_true, y_pred)
    fmi = fowlkes_mallows_score(y_true, y_pred)
    hom = homogeneity_score(y_true, y_pred)
    com = completeness_score(y_true, y_pred)
    v = v_measure_score(y_true, y_pred)
    return nmi, ari, ami, fmi, hom, com, v

def GetCluster(X, res, n):
    adata0=sc.AnnData(X)
    if adata0.shape[0]>200000:
       np.random.seed(adata0.shape[0])#set seed 
       adata0=adata0[np.random.choice(adata0.shape[0],200000,replace=False)] 
    sc.pp.neighbors(adata0, n_neighbors=n, use_rep="X")
    sc.tl.louvain(adata0,resolution=res)
    Y_pred_init=adata0.obs['louvain']
    Y_pred_init=np.asarray(Y_pred_init,dtype=int)
    if np.unique(Y_pred_init).shape[0]<=1:
        #avoid only a cluster
        exit("Error: There is only a cluster detected. The resolution:"+str(res)+"is too small, please choose a larger resolution!!")
    else: 
        print("Estimated n_clusters is: ", np.shape(np.unique(Y_pred_init))[0]) 
    return(np.shape(np.unique(Y_pred_init))[0])


rna = sc.read_h5ad('./multiome_bmmc_site1_or_donor1_RNA.h5ad')
atac = sc.read_h5ad('./multiome_bmmc_site1_or_donor1_ATAC_pmat.h5ad')
count_rna = rna.X
rna.obs['barcodes'] = rna.obs_names
barcodes_rna = rna.obs['barcodes']
print(rna.var)
features_rna = rna.var_names
print(features_rna)
snare_mrna = SingleData("GeneExpr", "SNARE-seq", features_rna, count_rna, barcodes_rna)

count_atac = atac.X
atac.obs['barcodes'] = atac.obs_names
barcodes_atac = atac.obs['barcodes']

features_atac = atac.var_names
snare_atac = SingleData("ChromAccess", "SNARE-seq", features_atac, count_atac, barcodes_atac)


multi_dt = MultiomicDataset.from_singledata(snare_atac, snare_mrna)
print(multi_dt)

model = Cobolt(dataset=multi_dt, n_latent=10, lr=1e-4)
model.train(num_epochs=5)

model.calc_all_latent()
latent = model.get_all_latent()

cell_name = atac.obs['cell_type']
cell_type, y = np.unique(cell_name, return_inverse=True)
n_clusters = int(max(y) - min(y) + 1)  

kmeans = KMeans(n_clusters = n_clusters, n_init=20)
latent = latent[0]
print(latent)
y_pred_z = kmeans.fit_predict(latent) 

nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
print('z for clustering, NMI:{:.4f}, ARI:{:.4f}, AMI:{:.4f}, FMI:{:.4f}, HOM:{:.4f}, COM:{:.4f}, V:{:.4f}'.format(nmi, ari, ami, fmi, hom, com, v))


t1=time()
t=t1-t0
print('time:{:.4f}'.format(t))
