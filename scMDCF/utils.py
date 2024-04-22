from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, adjusted_mutual_info_score
from sklearn.metrics import fowlkes_mallows_score, homogeneity_score, completeness_score, v_measure_score, silhouette_score
import scanpy as sc
import numpy as np
import pandas as pd
import h5py

def eva(y_true, y_pred):
    nmi = nmi_score(y_true, y_pred)
    ari = ari_score(y_true, y_pred)
    ami = adjusted_mutual_info_score(y_true, y_pred)
    fmi = fowlkes_mallows_score(y_true, y_pred)
    hom = homogeneity_score(y_true, y_pred)
    com = completeness_score(y_true, y_pred)
    v = v_measure_score(y_true, y_pred)
    return nmi, ari, ami, fmi, hom, com, v

def eva_nolabel(data, y):
    db = davies_bouldin_score(data, y)
    ch = calinski_harabasz_score(data, y)
    asw = silhouette_score(data, y)
    return db, ch, asw

def read_data(file_path1, file_path2, file_type, label_file):
    if file_type=='h5ad':
        atac = sc.read_h5ad(file_path2)
        rna = sc.read_h5ad(file_path1)
        atac_X = np.array(atac.X.toarray())
        rna_X = np.array(rna.X.toarray())
        if label_file==None:
            cell_name = np.array(atac.obs["cell_type"])
        else:
            cell_name = pd.read_csv(label_file, usecols=[1])
        cell_type, y = np.unique(cell_name, return_inverse=True)
        print(y)
        cluster_number = int(max(y) - min(y) + 1)   
        adata_RNA = sc.AnnData(rna_X)
        adata_ATAC = sc.AnnData(atac_X)
    elif file_type=='h5':
        data_mat = h5py.File(file_path1)
        rna_X = np.array(data_mat['X1'])
        atac_X = np.array(data_mat['X1'])
        y = np.array(data_mat['Y'])
        data_mat.close()
        cluster_number = int(max(y) - min(y) + 1) 
        adata_RNA = sc.AnnData(rna_X)
        adata_ATAC = sc.AnnData(atac_X)
    elif file_type=='loom':
        adata_RNA=sc.read_loom(file_path1)
        adata_ATAC=sc.read_loom(file_path2)
        cell_name = pd.read_csv(label_file, usecols=[1])
        cell_type, y = np.unique(cell_name, return_inverse=True)
        cluster_number = int(max(y) - min(y) + 1)   

    return adata_RNA, adata_ATAC, cluster_number, y

def read_data_nolabel(file_path1, file_path2, file_type):
    if file_type=='h5ad':
        atac = sc.read_h5ad(file_path2)
        rna = sc.read_h5ad(file_path1)# Pbmc10k Chen-2019
        #atac_X = np.array(atac.X.toarray())#.toarray()
        #rna_X = np.array(rna.X.toarray())#.toarray()
        #adata_RNA = sc.AnnData(rna_X)
        #adata_ATAC = sc.AnnData(atac_X)
    elif file_type=='h5':
        data_mat = h5py.File(file_path1)
        rna_X = np.array(data_mat['X1'])
        atac_X = np.array(data_mat['X1'])
        data_mat.close()
        adata_RNA = sc.AnnData(rna_X)
        adata_ATAC = sc.AnnData(atac_X)
    elif file_type=='loom':
        atac = sc.read_loom(file_path2)
        rna = sc.read_loom(file_path1)
    return adata_RNA, adata_ATAC

def normalize( adata, filter_min_counts=True, size_factors=True, highly_genes=None,
               normalize_input=False, logtrans_input=True):

    if filter_min_counts:
        sc.pp.filter_genes(adata, min_counts=1)#1
        sc.pp.filter_cells(adata, min_counts=1)#1
    
    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata
    
    if size_factors:
        sc.pp.normalize_per_cell(adata)
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0
    
    if logtrans_input:
        sc.pp.log1p(adata)
    
    if highly_genes != None:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = highly_genes, subset=True)

    if normalize_input:
        sc.pp.scale(adata)

    return adata

def GetCluster(adata0, res, n):
    #adata0=sc.AnnData(X)
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