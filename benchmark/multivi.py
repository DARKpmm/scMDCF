#!/usr/bin/python
import sys 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '6'

from anndata import AnnData
import mudata as md
from scvi.model import TOTALVI
import anndata as ad
import scipy
import resource
import numpy as np
import pandas as pd
import scipy.io as sio
import scanpy as sc
from copy import deepcopy
import timeit
from utils import *
from sklearn.cluster import KMeans
import h5py
from scipy.sparse import csr_matrix

#==== method specific ==== 
import scvi

scvi.settings.seed = 420

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
        var_v.append(pd.read_csv(os.path.join(dir_path,i+".tsv"), sep = '\t'))#csv
    var_df = pd.concat(var_v, axis=1)
    #var_df = var_df.set_axis(var_list, axis='columns')
    
    obs_v = []
    for i in obs_list:
        obs_v.append(pd.read_csv(os.path.join(dir_path,i+".csv"), sep = '\t', header = None))#csv
    obs_df = pd.concat(obs_v, axis=1)
    obs_df = obs_df.set_axis(obs_list, axis='columns')
    if(counts.shape[0]!=obs_df.shape[0]): counts=deepcopy(counts.transpose())
    adata = AnnData(counts,obs=obs_df,var=var_df)
    adata.obs = adata.obs.set_axis(list(adata.obs["barcodes"]),axis="index")
    adata.var= adata.var.set_axis(list(adata.var['feature_name']),axis="index")
    adata.var["modality"] = mod_key
    return(adata)

def run_multivi_fn(in_dir, out_dir,label_file=None):
    
    # save latent representation and model 
    start = timeit.default_timer()
    scvi.settings.seed = 420
    adata_prna = read_mtx_folder(os.path.join(in_dir,"paired_RNA/"),
                                       "Gene Expression",
                                       ["gene"],
                                       ["barcodes"])
    
    adata_patac = read_mtx_folder(os.path.join(in_dir,"paired_ATAC/"),
                                       "Peaks",
                                       ["peak"],
                                       ["barcodes"])
    cell_name = pd.read_csv('/data1/chengyue/gm/fig2/GSE214979/id.csv', usecols=[1])
    cell_type, y = np.unique(cell_name, return_inverse=True)
    adata_prna.var_names_make_unique()
    adata_prna.var_names = ["rna_" + str(gene) for gene in adata_prna.var_names]
    adata_patac.var_names = ["atac_" + str(gene) for gene in adata_patac.var_names]
    adata_patac.obs['modality']='Peaks'
    adata_patac.var['modality']='Peaks'
    
    adata_prna.obs['modality']='Gene Expression'
    adata_prna.var['modality']='Gene Expression'
    # Horizontally stack two modalities of paired dataset 
    adata_paired = AnnData(scipy.sparse.hstack((deepcopy(adata_prna.X), 
                                            deepcopy(adata_patac.X)), 
                                           format='csr'),
                           obs = deepcopy(adata_prna.obs),
                           var = pd.concat([deepcopy(adata_prna.var[["modality"]]),deepcopy(adata_patac.var[["modality"]])]))
    # organize_mulitome_anndatats: concatenate paired and two unpaired anndata
    adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired)#, adata_urna, adata_uatac)
    # # gene features need to be before chromatin peaks (algorithm assumption)
    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
    sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
    # setup batch annotation
    scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key='modality')
    # setup model 
    mvi = scvi.model.MULTIVI(
        adata_mvi,
        n_genes=(adata_mvi.var['modality']=='Gene Expression').sum(),
        n_regions=(adata_mvi.var['modality']=='Peaks').sum(), # Peaks
    )
    # train 
    mvi.train()
    os.makedirs(out_dir, exist_ok=True)
    # get latent representation 
    adata_mvi.obsm["MultiVI_latent"] = mvi.get_latent_representation()
    print(adata_mvi.obsm["MultiVI_latent"])
    adata_mvi.obs = adata_mvi.obs.set_axis([s. split("_", 1)[0] for s in adata_mvi.obs.index], axis='index')

    # extract latent representation
    res_df = pd.DataFrame(adata_mvi.obsm['MultiVI_latent'],index=adata_mvi.obs.index)
    
    n_clusters = int(max(y) - min(y) + 1)   
    
    kmeans = KMeans(n_clusters = n_clusters, n_init=20)
    y_pred_z = kmeans.fit_predict(res_df)   
    csv_out = os.path.join(out_dir, "multivi_result.csv")
    res_df.to_csv(csv_out)
    pred_out = os.path.join(out_dir, "multivi_pred.csv")
    np.savetxt(pred_out, y_pred_z)
    nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
    db, ch, asw = eva_nolabel(res_df, y_pred_z)
    print('nmi, ari, ami, fmi, hom, com, v:', nmi, ari, ami, fmi, hom, com, v)
    print('db, ch, asw:', db, ch, asw)
    
    # save latent representation and model 
    os.makedirs(os.path.join(out_dir,"multivi"), exist_ok=True)
    os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    
    model_out = os.path.join(out_dir,"multivi","trained_multivi")
    mvi.save(model_out, overwrite=True)
    stop = timeit.default_timer()
    print('Time(s): ', stop - start)  
    # record time 
    runtime_out = os.path.join(out_dir,"runtime","multivi_runtime.txt")
    print(stop - start,  file=open(runtime_out, 'w'))
    print("------ Done ------")
    
    stop = timeit.default_timer()
    max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"max_mem: , {max_mem/1024/1024:.2f}GB")
    print('Time(s): ', stop - start)  

print("argument 1:",sys.argv[1])
print("argument 2:",sys.argv[2])
#print("argument 3:",sys.argv[3])

mvi = run_multivi_fn(sys.argv[1],sys.argv[2])#,sys.argv[3])#,sys.argv[4])

# igraph leidenbase