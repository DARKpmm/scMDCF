#!/usr/bin/python

import sys 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))

#==== method specific ==== 
import networkx as nx
import scglue
from itertools import chain
import seaborn as sns
from matplotlib import rcParams
import h5py

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

from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, adjusted_mutual_info_score
from sklearn.metrics import fowlkes_mallows_score, homogeneity_score, completeness_score, v_measure_score, silhouette_score
from sklearn.cluster import KMeans
import anndata as ad

os.environ["CUDA_VISIBLE_DEVICES"] = "1"
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
        var_v.append(pd.read_csv(os.path.join(dir_path,i+".csv"), sep = '\t'))
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

def run_glue_fn(in_dir, out_dir):
    start = timeit.default_timer()

    rna = read_mtx_folder(os.path.join(in_dir,"paired_RNA/"),
                                       "Gene Expression",
                                       ["gene"],
                                       ["barcodes"])
    
    atac = read_mtx_folder(os.path.join(in_dir,"paired_ATAC/"),
                                       "Peaks",
                                       ["peak"],
                                       ["barcodes"])

    cell_name = pd.read_csv('./GSE214979/id.csv')
    cell_type, y = np.unique(cell_name, return_inverse=True)
    cluster_number = int(max(y) - min(y) + 1)   
    
    rna.obs['dataset'] = 'multiomeRNA'
    atac.obs['dataset'] = 'multiomeATAC'

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir,"glue"), exist_ok=True)
    os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    
    # preprocessing of scRNA
    rna.layers["counts"] = rna.X.copy()
    sc.pp.filter_genes(rna, min_cells=3)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")
    
    # preprocessing of snATAC
    sc.pp.filter_genes(atac,min_counts=1)
    sc.pp.normalize_total(atac)
    sc.pp.log1p(atac)
    sc.pp.scale(atac)
    sc.tl.pca(atac, n_comps=100, svd_solver="auto")
    scglue.data.lsi(atac, n_components=100, n_iter=15)

    # build graph
    scglue.data.get_gene_annotation(
        # this works for human hg38 genome-build 
        rna, gtf="./gencode.v47.annotation.gtf.gz", 
        gtf_by="gene_name"
    )
    rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
    
    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    
    atac_chrs = atac.var['chrom'].value_counts().index.tolist()
    row_keep = rna.var_names[rna.var['chrom'].isin(atac_chrs).tolist()]
    row_keep = list(set(row_keep))
    #rna = rna.loc[:, row_keep].copy()
    rna.var_names_make_unique()
    rna = rna[:,row_keep].copy()
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    scglue.graph.check_graph(guidance, [rna, atac])
    
    # prepare for training 
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi"#X_lsi
    )

    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()
    
    # GLUE training 
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        fit_kws={"directory": os.path.join(out_dir,"glue")}
    )
    
    dx = scglue.models.integration_consistency(
        glue, {"rna": rna, "atac": atac}, guidance_hvf
    )
    
    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    combined = ad.concat([rna, atac])

    # extract latent representation
    res_df = pd.DataFrame(combined.obsm['X_glue'],index=combined.obs.index)
    # set column names as latent_x 
    res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")
    res_df['dataset'] = combined.obs['dataset']
    res_df['dataset'] = res_df['dataset'].astype("string")
    
    # save latent representation and model
    
    csv_out = os.path.join(out_dir, "glue","glue_result.csv")
    res_df.to_csv(csv_out)
    model_out = os.path.join(out_dir,"glue","glue.dill")
    glue.save(model_out)
    stop = timeit.default_timer()
    
    n_clusters = cluster_number

    kmeans = KMeans(n_clusters = n_clusters, n_init=20)
    y_pred_z = kmeans.fit_predict(rna.obsm["X_glue"])   
    rna_z = os.path.join(out_dir,"glue","rna_z.csv")
    rna_pred = os.path.join(out_dir,"glue","rna_pred.csv")
    np.savetxt(rna_pred, y_pred_z)
    np.savetxt(rna_z, rna.obsm["X_glue"])
    print(y_pred_z)
    nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
    
    db, ch, asw = eva_nolabel(rna.obsm["X_glue"], y_pred_z)
    print('RNA: nmi, ari, ami, fmi, hom, com, v:', nmi, ari, ami, fmi, hom, com, v)
    print('RNA: db, ch, asw:', db, ch, asw)

    kmeans = KMeans(n_clusters = n_clusters, n_init=20)
    y_pred_z = kmeans.fit_predict(atac.obsm["X_glue"])   
    atac_z = os.path.join(out_dir,"glue","atac_z.csv")
    atac_pred = os.path.join(out_dir,"glue","atac_pred.csv")
    np.savetxt(atac_pred, y_pred_z)
    np.savetxt(atac_z, atac.obsm["X_glue"])
    nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
    
    db, ch, asw = eva_nolabel(atac.obsm["X_glue"], y_pred_z)
    print('ATAC: nmi, ari, ami, fmi, hom, com, v:', nmi, ari, ami, fmi, hom, com, v)
    print('ATAC: db, ch, asw:', db, ch, asw)
    print('Time(s): ', stop - start)  
    # record time 
    runtime_out = os.path.join(out_dir,"runtime","glue_runtime.txt")
    print(stop - start,  file=open(runtime_out, 'w'))
    print("------ Done ------")
    

run_glue_fn(sys.argv[1],sys.argv[2])