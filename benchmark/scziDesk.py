from preprocess import *
from network import *
from utils import *
import argparse
import os
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics import fowlkes_mallows_score, homogeneity_score, completeness_score, v_measure_score, adjusted_mutual_info_score
from preprocess import normalize
import time
import numpy as np
import scanpy as sc
import tensorflow as tf

def cluster_acc(y_true, y_pred):
    y_true = y_true.astype(np.int64)
    assert y_pred.size == y_true.size
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    from sklearn.utils.linear_assignment_ import linear_assignment
    ind = linear_assignment(w.max() - w)
    return sum([w[i, j] for i, j in ind]) * 1.0 / y_pred.size

if __name__ == "__main__":
    start=time.time()

    random_seed = [1111, 2222, 3333, 4444, 5555, 6666, 7777, 8888, 9999, 10000]

    parser = argparse.ArgumentParser(description="train", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dataname", default = "Qx_Limb_Muscle", type = str)
    parser.add_argument("--distribution", default = "ZINB")
    parser.add_argument("--self_training", default = True)
    parser.add_argument("--dims", default = [500, 256, 64, 32])
    parser.add_argument("--highly_genes", default = 500)
    parser.add_argument("--alpha", default = 0.001, type = float)
    parser.add_argument("--gamma", default = 0.001, type = float)
    parser.add_argument("--learning_rate", default = 0.0001, type = float)
    parser.add_argument("--random_seed", default = random_seed)
    parser.add_argument("--batch_size", default = 256, type = int)
    parser.add_argument("--update_epoch", default = 10, type = int)
    parser.add_argument("--pretrain_epoch", default = 1000, type = int)
    parser.add_argument("--funetrain_epoch", default = 2000, type = int)
    parser.add_argument("--t_alpha", default = 1.0)
    parser.add_argument("--noise_sd", default = 1.5)
    parser.add_argument("--error", default = 0.001, type = float)
    parser.add_argument("--gpu_option", default = "0")

    args = parser.parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu_option
    adata = sc.read_h5ad('./multiome_bmmc_site1_or_donor1_RNA.h5ad')
    cell_name = adata.obs['cell_type']
    cell_type, Y = np.unique(cell_name, return_inverse=True)
    Y=list(map(int, Y))
    # # print(Y)
    X = adata.X.astype(np.float32)
    count_X = X
    cluster_number = int(max(Y) - min(Y) + 1)  
    
    adata  = normalize( adata, filter_min_counts=True, highly_genes=args.highly_genes,
                        size_factors=True, normalize_input=False, 
                        logtrans_input=True ) 
    size_factor = np.array(adata.obs.size_factors).reshape(-1, 1).astype(np.float32)
    result = []
    

    for seed in args.random_seed:
        np.random.seed(seed)
        tf.compat.v1.reset_default_graph()
        chencluster = autoencoder(args.dataname, args.distribution, args.self_training, args.dims, cluster_number, args.t_alpha,
                                  args.alpha, args.gamma, args.learning_rate, args.noise_sd)
        chencluster.pretrain(X, count_X, size_factor, args.batch_size, args.pretrain_epoch, args.gpu_option)
        chencluster.funetrain(X, count_X, size_factor, args.batch_size, args.funetrain_epoch, args.update_epoch, args.error)
        
        ARI = np.around(adjusted_rand_score(Y, chencluster.Y_pred), 5)
        NMI = np.around(normalized_mutual_info_score(Y, chencluster.Y_pred), 5)
        AMI = adjusted_mutual_info_score(Y, chencluster.Y_pred)
        FMI = fowlkes_mallows_score(Y, chencluster.Y_pred)
        HOM = homogeneity_score(Y, chencluster.Y_pred)
        COM = completeness_score(Y, chencluster.Y_pred)
        V = v_measure_score(Y, chencluster.Y_pred)
        print('Evaluating cells: NMI= %.4f, ARI= %.4f, AMI= %.4f, FMI= %.4f, HOM= %.4f, COM= %.4f, V= %.4f' % (NMI, ARI, AMI, FMI, HOM, COM, V))
        db, ch, asw = eva_nolabel(chencluster.latent_repre, chencluster.Y_pred)
        print('z for clustering epoch:{:}, DB:{:.4f}, CH:{:.4f}, ASW:{:.4f}'.format(db, ch, asw))
        print("over")
        
    print('Evaluating cells: NMI= %.4f, ARI= %.4f, AMI= %.4f, FMI= %.4f, HOM= %.4f, COM= %.4f, V= %.4f' % (nmi, ari, ami, fmi, hom, com, v))
    end=time.time()
    print(end-start)



#jietu Plasschaert











