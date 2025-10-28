import anndata
import scanpy as sc
import sccross
import pandas as pd
from matplotlib import rcParams
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
import numpy as np
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, adjusted_mutual_info_score
from sklearn.metrics import fowlkes_mallows_score, homogeneity_score, completeness_score, v_measure_score, silhouette_score
from sklearn.cluster import KMeans
import timeit

start = timeit.default_timer()
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

rcParams["figure.figsize"] = (4, 4)
rna = anndata.read_h5ad("./GSE214979/rna_preprocessed_n.h5ad")
atac = anndata.read_h5ad("./GSE214979/atac_preprocessed_n.h5ad")
atac2rna = anndata.read_h5ad("./GSE214979/atac2rna_n.h5ad")#5000
label_dir = './GSE214979/id.csv'
atac.uns['gene'] = atac2rna

# label_dir = ""
rna.obs['domain'] = 'scRNA-seq'
atac.obs['domain'] = 'scATAC-seq'

sccross.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer = 'counts',
     use_rep="X_pca"
)

sccross.models.configure_dataset(
    atac, "NB", use_highly_variable=False,
    use_rep="X_lsi"
)
# = rna.obsm['X_pca'].iloc[:,-50:]
sccross.data.mnn_prior([rna,atac])

cross = sccross.models.fit_SCCROSS(
    {"rna": rna, "atac": atac}, 
    fit_kws={"directory": "sccross"}
)

rna.obsm["X_cross"] = cross.encode_data("rna", rna)
atac.obsm["X_cross"] = cross.encode_data("atac", atac)

combined = anndata.concat([rna, atac])


latent = combined.obsm['X_cross']

cell_name = pd.read_csv(label_dir, usecols=[1])
cell_type, y = np.unique(cell_name, return_inverse=True)
n_clusters = int(max(y) - min(y) + 1)   

kmeans = KMeans(n_clusters = n_clusters, n_init=20)

y_pred_z = kmeans.fit_predict(latent)   
nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
db, ch, asw = eva_nolabel(latent, y_pred_z)
print('Combined: nmi, ari, ami, fmi, hom, com, v:', nmi, ari, ami, fmi, hom, com, v)
print('Combined: db, ch, asw:', db, ch, asw)
np.savetxt('./sccross-output/combine-pred.csv', y_pred_z)
np.savetxt('./sccross-output/combine-latent.csv', latent)

latent_rna = rna.obsm["X_cross"]
y_pred_z = kmeans.fit_predict(latent_rna)   
nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
db, ch, asw = eva_nolabel(latent_rna, y_pred_z)
print('RNA: nmi, ari, ami, fmi, hom, com, v:', nmi, ari, ami, fmi, hom, com, v)
print('RNA: db, ch, asw:', db, ch, asw)
np.savetxt('./hpap/sccross-output/rna-pred.csv', y_pred_z)
np.savetxt('./hpap/sccross-output/rna-latent.csv', latent_rna)

latent_atac = atac.obsm["X_cross"]
y_pred_z = kmeans.fit_predict(latent_atac)   
nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred_z)
db, ch, asw = eva_nolabel(latent_atac, y_pred_z)
print('ATAC: nmi, ari, ami, fmi, hom, com, v:', nmi, ari, ami, fmi, hom, com, v)
print('ATAC: db, ch, asw:', db, ch, asw)
np.savetxt('./hpap/sccross-output/atac-pred.csv', y_pred_z)
np.savetxt('./hpap/sccross-output/atac-latent.csv', latent_atac)

stop = timeit.default_timer()
print('Time(s): ', stop - start)  