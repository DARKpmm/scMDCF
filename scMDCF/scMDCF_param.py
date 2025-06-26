import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import argparse
import torch
#from sklearn.metrics import silhouette_samples, silhouette_score

from utils import read_data, normalize#, get_adj, eva, read_dataset, CLR_normalization, GetCluster
from layer import scDEFR_multi
from train_param import pre_train, alt_train
from time import time
import warnings
warnings.filterwarnings('ignore')
from hyperopt import fmin, tpe, hp,space_eval,rand,Trials,partial,STATUS_OK


if __name__ == "__main__":

    para = {"feature_number":[1000, 2000, 2500, 3000, 4000],
         "lamb":[0.1, 0.5, 1, 5, 10],
         "latent":[8, 16, 32],
         "enc1":[1024, 512, 256],
         "enc2":[256, 128, 64],
         "weight1":[0.1, 1, 5, 10],
         "weight2":[0.1, 1, 5, 10],
         "weight3":[0.1, 1, 5, 10],
         "weight4":[0.1, 1, 5, 10]}

    space = {"feature_number":hp.choice("feature_number", (1000, 2000, 2500, 3000, 4000)),
         "lamb":hp.choice("lamb", (0.1, 0.5, 1, 5, 10)),
         "latent":hp.choice("latent", (8, 16, 32)),
         "enc1":hp.choice("enc1", (1024, 512, 256, 1500, 2000)),
         "enc2":hp.choice("enc2", (256, 128, 64)),
         "weight1":hp.choice("weight1", (0.1, 1, 5, 10)),
         "weight2":hp.choice("weight2", (0.1, 1, 5, 10)),
         "weight3":hp.choice("weight3", (0.1, 1, 5, 10)),
         "weight4":hp.choice("weight4", (0.1, 1, 5, 10))
        }
    
    def hyperpara(argsDict):
        parser = argparse.ArgumentParser(description='train', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('--file_path1', default='/home/chengyue/data/multi-omics/mouse_skin_shareseq_rna_10k.h5ad')#peripheral_blood_rna.h5ad  pbmc_spector.h5 Pbmc10k-RNA
        parser.add_argument('--file_path2', default='/home/chengyue/data/multi-omics/mouse_skin_shareseq_atac_10k.h5ad')
        parser.add_argument('--label_file', default=None)#'/home/chengyue/data/multi-omics/tea/rna_pbmc.csv')
        parser.add_argument('--save_results', default='False', type=bool)
        parser.add_argument('--file_type', default='h5ad', type=str)
        parser.add_argument('--model_file', default='/home/chengyue/data/multi-omics/test1.pth.tar')
        parser.add_argument("--highly_genes", default = 2500, type = int)#####
        parser.add_argument("--lr_pre", default = "1e-2", type = float)#2
        parser.add_argument("--lr_alt", default = "1e-3", type = float)#3
        parser.add_argument("--epoch_pre", default = "200", type = int)#2
        parser.add_argument("--epoch_alt", default = "200", type = int)#2
        parser.add_argument("--beginkl", default = "200", type = int)
        parser.add_argument("--device", default='cuda:2', type=str)
        parser.add_argument("--enc1", default = "512", type = int)#512#####
        parser.add_argument("--enc2", default = "64", type = int)#64#####
        parser.add_argument("--zdim", default = "8", type = int)######
        parser.add_argument("--batch_size", default = "256", type = int)
        parser.add_argument("--gamma", default = "1.", type = float)
        parser.add_argument("--lamb", default = "0.1", type = float)#####
        parser.add_argument("--weight1", default = "0.1", type = float)
        parser.add_argument("--weight2", default = "10", type = float)
        parser.add_argument("--weight3", default = "1", type = float)
        parser.add_argument("--weight4", default = "5", type = float)
        args = parser.parse_args()

        args.highly_genes = argsDict['feature_number']
        args.lamb = argsDict['lamb'] 
        args.zdim = argsDict['latent']
        args.enc1 = argsDict['enc1']
        args.enc2 = argsDict['enc2']
        args.weight1 = argsDict['weight1']
        args.weight2 = argsDict['weight2']
        args.weight3 = argsDict['weight3']
        args.weight4 = argsDict['weight4']

        adata_RNA, adata_ATAC, cluster_number, y = read_data(args.file_path1, args.file_path2, args.file_type, args.label_file)
        adata_RNA = normalize(adata_RNA, highly_genes=args.highly_genes, normalize_input=True)
        adata_ATAC = normalize(adata_ATAC, highly_genes=args.highly_genes, normalize_input=True)

        args.layere_view = [adata_RNA.X.shape[1], args.enc1, args.enc2]
        args.layere_adt_view = [adata_ATAC.X.shape[1], args.enc1, args.enc2]
        args.layerd_view = [args.zdim, args.enc2, args.enc1, adata_RNA.X.shape[1]]
        args.layerd_adt_view = [args.zdim, args.enc2, args.enc1, adata_ATAC.X.shape[1]]

        model = scDEFR_multi(args).to(args.device)
        x_rna, x_atac = torch.from_numpy(adata_RNA.X).to(args.device).float(), torch.from_numpy(adata_ATAC.X).to(args.device).float()
        pre_train(args, model, x_rna, x_atac, y)

        nmi = alt_train(args, model, x_rna, x_atac, y)
    
    algo = partial(tpe.suggest,n_startup_jobs=1)
    best = fmin(hyperpara,space,algo=algo,max_evals=10)
    for key in para.keys():
        para[key]=para[key][best[key]]
    print(best)
    print(para)
    print(hyperpara(para))




