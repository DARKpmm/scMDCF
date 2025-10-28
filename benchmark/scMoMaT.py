#!/usr/bin/python

import sys 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))

#==== method specific ==== 
import scmomat
import torch

#from matplotlib import rcParams
from anndata import AnnData
import anndata as ad
import scipy
import numpy as np
import pandas as pd
import scipy.io as sio
import os
import scanpy as sc
from copy import deepcopy
#from utils_eval import read_mtx_folder, write_adata
import timeit

# ------ Helper Functions -----
# copied from https://github.com/PeterZZQ/scMoMaT/blob/main/scmomat/utils.py
# slight adjustment to preprocess() line 91 of this document. An error occured due to third dimension introduced in the previous calculations. 
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, adjusted_mutual_info_score
from sklearn.metrics import fowlkes_mallows_score, homogeneity_score, completeness_score, v_measure_score, silhouette_score
from sklearn.cluster import KMeans
import anndata as ad
def read_mtx_folder(dir_path,mod_key,var_list,obs_list):
    mtx_path = []
    tsv_path = []
    for file in os.listdir(dir_path):
        if file.endswith(".mtx"):
            mtx_path.append(os.path.join(dir_path, file))
    if(len(mtx_path)==0):
        print("no .mtx file found, exiting function")
        return
    counts = sio.mmread(mtx_path[0])
    #if type(counts) != "scipy.sparse.csr.csr_matrix": counts = counts.tocsr()
    if not isinstance(counts, scipy.sparse.csr_matrix): 
        counts = counts.tocsr()
    
    var_v = []
    for i in var_list:
        var_v.append(pd.read_csv(os.path.join(dir_path,i+".tsv"), sep = '\t'))
    var_df = pd.concat(var_v, axis=1)
    #var_df = var_df.set_axis(var_list, axis='columns')
    
    obs_v = []
    for i in obs_list:
        obs_v.append(pd.read_csv(os.path.join(dir_path,i+".csv"), sep = '\t'))
    obs_df = pd.concat(obs_v, axis=1)
    obs_df = obs_df.set_axis(obs_list, axis='columns')
    if(counts.shape[0]!=obs_df.shape[0]): counts=deepcopy(counts.transpose())
    adata = AnnData(counts,obs=obs_df,var=var_df)
    adata.obs = adata.obs.set_axis(list(adata.obs["barcodes"]),axis="index")
    adata.var= adata.var.set_axis(list(adata.var['feature_name']),axis="index")
    adata.var["modality"] = mod_key
    return(adata)

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

def quantile_norm(X):
    from scipy import stats
    """Normalize the columns of X to each have the same distribution.
    Given an expression matrix (microarray data, read counts, etc) of M genes
    by N samples, quantile normalization ensures all samples have the same
    spread of data (by construction).
    The data across each row are averaged to obtain an average column. Each
    column quantile is replaced with the corresponding quantile of the average
    column.
    Parameters
    ----------
    X : 2D array of float, shape (M, N)
        The input data, with M rows (genes/features) and N columns (samples).
    Returns
    -------
    Xn : 2D array of float, shape (M, N)
        The normalized data.
    """
    # compute the quantiles
    quantiles = np.mean(np.sort(X, axis=0), axis=1)

    # compute the column-wise ranks. Each observation is replaced with its
    # rank in that column: the smallest observation is replaced by 1, the
    # second-smallest by 2, ..., and the largest by M, the number of rows.
    ranks = np.apply_along_axis(stats.rankdata, 0, X)

    # convert ranks to integer indices from 0 to M-1
    rank_indices = ranks.astype(int) - 1

    # index the quantiles for each rank with the ranks matrix
    Xn = quantiles[rank_indices]

    return(Xn)

def quantile_norm_log(X, log = True):
    if log:
        logX = np.log1p(X)
    else:
        logX = X
    logXn = quantile_norm(logX)
    return logXn


def preprocess(counts, modality = "RNA", log = True):
    """\
    Description:
    ------------
        Preprocess the dataset, for count, interaction matrices
    
    Parameters:
    -------------
        counts: count matrix
        modality: modality that the matrix belong to, can be ``ATAC'', ``RNA'', ``protein''
        log: log transform the data or not
    """
    if modality == "ATAC":
        # make binary, maximum is 1
        counts = (counts > 0).astype(float) 

    else:
        # other cases, e.g. Protein, RNA, etc
        counts = quantile_norm_log(counts, log = log)[:,:,0]
        counts = counts/np.max(counts)

    return counts
# ------ Helper Functions End -----
def run_scmomat_fn(in_dir, out_dir, label_dir=None):#in_dir, 
    start = timeit.default_timer()
    
    rna_adata = read_mtx_folder(os.path.join(in_dir,"paired_RNA/"),
                                       "Gene Expression",
                                       ["gene"],
                                       ["barcodes"])
    
    atac_adata = read_mtx_folder(os.path.join(in_dir,"paired_ATAC/"),
                                       "Peaks",
                                       ["peak"],
                                       ["barcodes"])
    
    cell_bc = np.concatenate([rna_adata.obs.index.values.squeeze()])
    rna_adata.obs['dataset'] = 'multiomeRNA'
    atac_adata.obs['dataset'] = 'multiomeATAC'
    
    # find the hvg for RNA modality
    rna_adata.layers["counts"] = rna_adata.X.copy()
    atac_adata.layers["counts"] = atac_adata.X.copy()
    sc.pp.normalize_total(rna_adata, target_sum=1e4)
    sc.pp.log1p(rna_adata)
    sc.pp.highly_variable_genes(rna_adata,n_top_genes=5000)
    hvg_prna = rna_adata.var['highly_variable'].index[rna_adata.var['highly_variable']].tolist()
    hvg_rna_common = hvg_prna
    
    ## calculate the conversion between ATAC to pseudoRNA counts
    A_long = pd.read_csv('./GSE214979/pbmc_peak_GxR.csv')
    A_long_sel = A_long[A_long.loc[:, 'gene.name'].isin(hvg_rna_common)]
    A_sel = pd.crosstab(A_long_sel.loc[:, 'peak'], A_long_sel.loc[:, 'gene.name'])

    # filter for peaks that are within 2000bp of TSS or along the gene body, as well as genes that have at least one peak
    gene_final = A_sel.columns.values.squeeze()
    peak_final = A_sel.index.values.squeeze()
    
    rna = rna_adata[:,rna_adata.var_names.isin(gene_final)].layers['counts'].todense()

    rna = preprocess(rna, modality = "RNA", log = True) 

    atac = atac_adata[:,atac_adata.var_names.isin(peak_final)].X.todense()

    atac = preprocess(atac, modality = "ATAC")
    
    counts = {"rna":[rna], "atac": [atac] }
    
    feats_name = {"rna": gene_final, "atac": peak_final}
    
    counts["feats_name"] = feats_name
    
    counts["nbatches"] = 1

    # Running scMoMaT model
    device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
    lamb = 0.001
    batchsize = 0.1
    # running seed
    seed = 0
    # number of latent dimensions
    K = 30
    interval = 1000
    T = 4000
    lr = 1e-2
    # 1st stage training, learning cell factors
    model = scmomat.scmomat_model(counts = counts, K = K, device = device, lr=lr, lamb=lamb, interval=interval)
    losses = model.train_func(T = T)

    # extract cell factors/latent representations
    
    zs = model.extract_cell_factors()
    res_df = pd.DataFrame(np.concatenate(zs, axis=0),
                          index=cell_bc)
    
    
    res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")
    res_df = res_df.squeeze()
    print(res_df)

    cell_name = pd.read_csv(label_dir, usecols=[1])
    cell_type, y = np.unique(cell_name, return_inverse=True)
    n_clusters = int(max(y) - min(y) + 1)   

    kmeans = KMeans(n_clusters = n_clusters, n_init=20)
    y_pred_z = kmeans.fit_predict(res_df)   
    nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
    db, ch, asw = eva_nolabel(res_df, y_pred_z)
    print('nmi, ari, ami, fmi, hom, com, v:', nmi, ari, ami, fmi, hom, com, v)
    print('db, ch, asw:', db, ch, asw)
    
    
    
    # save latent representation and model
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir,"scmomat"), exist_ok=True)
    os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    
    csv_out = os.path.join(out_dir, "scmomat","scmomat_result.csv")
    y_out = os.path.join(out_dir, "scmomat","scmomat_y.csv")
    res_df.to_csv(csv_out)
    y_pred_z = pd.DataFrame(y_pred_z)
    y_pred_z.to_csv(y_out)
    stop = timeit.default_timer()
    
    print('Time(s): ', stop - start)  
    # record time 
    runtime_out = os.path.join(out_dir,"runtime","scmomat_runtime.txt")
    print(stop - start,  file=open(runtime_out, 'w'))
    print("------ Done ------")

print("argument 1:",sys.argv[1])
print("argument 2:",sys.argv[2])

run_scmomat_fn(sys.argv[1],sys.argv[2])#,sys.argv[3])