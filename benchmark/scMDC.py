from time import time
import math, os
os.environ['OPENBLAS_NUM_THREADS'] = '2'
from sklearn import metrics
from sklearn.cluster import KMeans
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.nn import Parameter
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import silhouette_score
from scMDC import scMultiCluster
import numpy as np
import collections
import h5py
import scanpy as sc
from preprocess import read_dataset, normalize, clr_normalize_each_cell
from utils import *
import timeit
import resource
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, adjusted_mutual_info_score
from sklearn.metrics import fowlkes_mallows_score, homogeneity_score, completeness_score, v_measure_score, silhouette_score


def eva(y_true, y_pred):
    nmi = nmi_score(y_true, y_pred)
    ari = ari_score(y_true, y_pred)
    ami = adjusted_mutual_info_score(y_true, y_pred)
    fmi = fowlkes_mallows_score(y_true, y_pred)
    hom = homogeneity_score(y_true, y_pred)
    com = completeness_score(y_true, y_pred)
    v = v_measure_score(y_true, y_pred)
    return nmi, ari, ami, fmi, hom, com, v
if __name__ == "__main__":
    start = timeit.default_timer()
    # setting the hyper parameters
    import argparse
    parser = argparse.ArgumentParser(description='train',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_clusters', default=27, type=int)#27
    parser.add_argument('--cutoff', default=0.5, type=float, help='Start to train combined layer after what ratio of epoch')
    parser.add_argument('--batch_size', default=256, type=int)
    parser.add_argument('--data_file', default='./spleen_lymph.h5') #SMAGESeq;GSE128639_BMNC
    parser.add_argument('--maxiter', default=5000, type=int)
    parser.add_argument('--pretrain_epochs', default=400, type=int)#400
    parser.add_argument('--gamma', default=.1, type=float,
                        help='coefficient of clustering loss')
    parser.add_argument('--tau', default=1., type=float,
                        help='fuzziness of clustering loss')                    
    parser.add_argument('--phi1', default=0.001, type=float,
                        help='coefficient of KL loss')
    parser.add_argument('--phi2', default=0.001, type=float,
                        help='coefficient of KL loss')
    parser.add_argument('--update_interval', default=1, type=int)
    parser.add_argument('--tol', default=0.001, type=float)
    parser.add_argument('--lr', default=1., type=float)
    parser.add_argument('--ae_weights', default=None)
    parser.add_argument('--save_dir', default='results/')
    parser.add_argument('--ae_weight_file', default='AE_weights_1.pth.tar')
    parser.add_argument('--resolution', default=0.2, type=float)
    parser.add_argument('--n_neighbors', default=30, type=int)
    parser.add_argument('--embedding_file', action='store_true', default=False)
    parser.add_argument('--prediction_file', action='store_true', default=False)
    parser.add_argument('-el','--encodeLayer', nargs='+', default=[256,64,32,16])
    parser.add_argument('-dl1','--decodeLayer1', nargs='+', default=[16,64,256])
    parser.add_argument('-dl2','--decodeLayer2', nargs='+', default=[16,20])
    parser.add_argument('--sigma1', default=2.5, type=float)
    parser.add_argument('--sigma2', default=1.5, type=float)
    parser.add_argument('--f1', default=1000, type=float, help='Number of mRNA after feature selection')
    parser.add_argument('--f2', default=2000, type=float, help='Number of ADT/ATAC after feature selection')
    parser.add_argument('--filter1', action='store_true', default=False, help='Do mRNA selection')
    parser.add_argument('--filter2', action='store_true', default=False, help='Do ADT/ATAC selection')
    parser.add_argument('--run', default=1, type=int)
    parser.add_argument('--device', default='cuda')
    args = parser.parse_args()
    print(args)
    
    atac = sc.read_h5ad('./GSE214979_atac.h5ad')#inhouse_ADT
    rna = sc.read_h5ad('./GSE214979_rna.h5ad')
    x2 = np.array(atac.X.toarray())
    x1 = np.array(rna.X.toarray())
    
    import pandas as pd
    cell_name = pd.read_csv('./GSE214979/id.csv', usecols=[1])
    cell_type, y = np.unique(cell_name, return_inverse=True)
    
    args.n_clusters=int(max(y)-min(y)+1)
    print(args.n_clusters)
    #Gene filter
    if args.filter1:
        importantGenes = geneSelection(x1, n=args.f1, plot=False)
        x1 = x1[:, importantGenes]
    if args.filter2:
        importantGenes = geneSelection(x2, n=args.f2, plot=False)
        x2 = x2[:, importantGenes]
        
    # preprocessing scRNA-seq read counts matrix
    adata1 = sc.AnnData(x1)
    adata1 = normalize(adata1,
                      size_factors=True,
                      normalize_input=True,
                      logtrans_input=True)
    
    adata2 = sc.AnnData(x2)
    input_size1 = adata1.n_vars
    input_size2 = adata2.n_vars

    encodeLayer = list(map(int, args.encodeLayer))
    decodeLayer1 = list(map(int, args.decodeLayer1))
    decodeLayer2 = list(map(int, args.decodeLayer2))
    
    model = scMultiCluster(input_dim1=input_size1, input_dim2=input_size2, tau=args.tau,
                        encodeLayer=encodeLayer, decodeLayer1=decodeLayer1, decodeLayer2=decodeLayer2,
                        activation='elu', sigma1=args.sigma1, sigma2=args.sigma2, gamma=args.gamma, 
                        cutoff = args.cutoff, phi1=args.phi1, phi2=args.phi2, device=args.device).to(args.device)
    
    
    if not os.path.exists(args.save_dir):
            os.makedirs(args.save_dir)
            
    t0 = time()
    if args.ae_weights is None:
        model.pretrain_autoencoder(X1=adata1.X, X_raw1=adata1.raw.X, sf1=adata1.obs.size_factors,
                X2=adata2.X, X_raw2=adata2.X, sf2=adata1.obs.size_factors, batch_size=args.batch_size, 
                epochs=args.pretrain_epochs, ae_weights=args.ae_weight_file)###有改
    else:
        if os.path.isfile(args.ae_weights):
            print("==> loading checkpoint '{}'".format(args.ae_weights))
            checkpoint = torch.load(args.ae_weights)
            model.load_state_dict(checkpoint['ae_state_dict'])
        else:
            print("==> no checkpoint found at '{}'".format(args.ae_weights))
            raise ValueError
    
    print('Pretraining time: %d seconds.' % int(time() - t0))
    
    #get k
    latent = model.encodeBatch(torch.tensor(adata1.X).to(args.device), torch.tensor(adata2.X).to(args.device))
    latent = latent.cpu().numpy()
    if args.n_clusters == -1:
       n_clusters = GetCluster(latent, res=args.resolution, n=args.n_neighbors)
    else:
       print("n_cluster is defined as " + str(args.n_clusters))
       n_clusters = args.n_clusters

    y_pred, _, _, _, _ = model.fit(X1=adata1.X, X_raw1=adata1.raw.X, sf1=adata1.obs.size_factors, 
        X2=adata2.X, X_raw2=adata2.X, sf2=adata1.obs.size_factors, y=y,
        n_clusters=n_clusters, batch_size=args.batch_size, num_epochs=args.maxiter, 
        update_interval=args.update_interval, tol=args.tol, lr=args.lr, save_dir=args.save_dir)
    print('Total time: %d seconds.' % int(time() - t0))
    
    
    final_latent = model.encodeBatch(torch.tensor(adata1.X).to(args.device), torch.tensor(adata2.X).to(args.device))
    final_latent = final_latent.cpu().numpy()
    #np.savetxt(args.save_dir + "/" + str(args.run) + "_embedding.csv", final_latent, delimiter=",")
    np.savetxt("./GSE214979/scmdc-output/z.csv", final_latent, delimiter=',')

    y_pred_ = best_map(y, y_pred)
    y_pred_ = best_map(y, y_pred) - 1
    np.savetxt("./GSE214979/scmdc-output/predy.csv", y_pred, delimiter=',')
    nmi_z, ari_z, ami_z, fmi_z, hom_z, com_z, v_z = eva(y, y_pred_)
    print('z for clustering, NMI:{:.4f}, ARI:{:.4f}, AMI:{:.4f}, FMI:{:.4f}, HOM:{:.4f}, COM:{:.4f}, V:{:.4f}'.format(nmi_z, ari_z, ami_z, fmi_z, hom_z, com_z, v_z))

    stop = timeit.default_timer()
    print('time: ', stop-start)
    max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"max_mem: , {max_mem/1024/1024:.2f}GB")