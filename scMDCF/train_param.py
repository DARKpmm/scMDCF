import torch
#import torch.nn as nn
from torch.optim import Adam, Adadelta
import torch.nn.functional as F
from sklearn.cluster import KMeans
import numpy as np
from utils import eva, target_distribution

def pre_train(args, model, X_RNA, X_ATAC, y):
    optimizer = Adam(filter(lambda p: p.requires_grad, model.parameters()), lr = args.lr_pre, amsgrad=True)
    for epoch in range(args.epoch_pre):
        z_RNA, z_ATAC, rec_RNA, rec_ATAC, z, _ = model(X_RNA, X_ATAC)#
        #p = target_distribution(q)
        loss_recrna = F.mse_loss(rec_RNA, X_RNA)
       
        loss_recatac = F.mse_loss(rec_ATAC, X_ATAC)#torch.spmm(adj, X_view2)
        
        cl_loss = model.crossview_contrastive_Loss(z_ATAC, z_RNA, lamb = args.lamb)
        loss = loss_recrna+loss_recatac#+0.1*cl_loss
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch+=1
        if epoch%10==0:
            print('epoch:{}, loss:{:.4f}, loss_RNA:{:.4f}, loss_ATAC:{:.4f}, loss_cl:{:.4f}'.format(epoch, loss, loss_recrna, loss_recatac, cl_loss))
    
    
    # kmeans = KMeans(n_clusters = args.n_clusters, n_init=20)
    # y_pred = kmeans.fit_predict(z.data.cpu().numpy())
    # nmi, ari = eva(y, y_pred)
#    torch.save(model.state_dict(), args.model_file)
    # print('NMI:{:.4f}, ARI:{:.4f}'.format(nmi, ari))
    # return z, y_pred

def alt_train(args, model, X_RNA, X_ATAC, y):
    #model.load_state_dict(torch.load(args.model_file, map_location=args.device))
    with torch.no_grad():
        _, _, _, _, z, _ = model(X_RNA, X_ATAC)#  
    kmeans = KMeans(n_clusters = args.n_clusters, n_init=20)
    y_pred = kmeans.fit_predict(z.data.cpu().numpy())
    nmi, ari, ami, fmi, hom, com, v = eva(y, y_pred)
    model.cluster_layer.data = torch.tensor(kmeans.cluster_centers_).to(args.device)
    print('z for clustering, NMI:{:.4f}, ARI:{:.4f}, AMI:{:.4f}, FMI:{:.4f}, HOM:{:.4f}, COM:{:.4f}, V:{:.4f}'.format(nmi, ari, ami, fmi, hom, com, v))
    
    optimizer = Adadelta(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr_alt, rho=.8)
    
    for epoch in range(args.epoch_alt):
        model.train()
        z_RNA, z_ATAC, rec_RNA, rec_ATAC, z, q = model(X_RNA, X_ATAC)#
        p = target_distribution(q)
        loss_recrna = F.mse_loss(rec_RNA, X_RNA)
        loss_recatac = F.mse_loss(rec_ATAC, X_ATAC)
        
        cl_loss = model.crossview_contrastive_Loss(z_ATAC, z_RNA, lamb = args.lamb)
        loss_clu = model.cluster_loss(args, p, q)
        
        loss=args.weight1*loss_recrna+args.weight2*loss_recatac+args.weight3*loss_clu+args.weight4*cl_loss#5;0.1;5
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch+=1
        if epoch%10==0:
            print('epoch:{}, loss:{:.4f}, loss_RNA:{:.4f}, loss_ATAC:{:.4f}, loss_kl:{:.4f}, loss_cl:{:.4f}'.format(epoch, loss, loss_recrna, loss_recatac, loss_clu, cl_loss))
            with torch.no_grad():
                z, q = encodeZ(model, X_RNA, X_ATAC)
            y_pred_q = torch.argmin(q, dim=1).data.cpu().numpy()
            kmeans = KMeans(n_clusters = args.n_clusters, n_init=20)
            y_pred_z = kmeans.fit_predict(z.data.cpu().numpy())          
            nmi_z, ari_z, ami_z, fmi_z, hom_z, com_z, v_z = eva(y, y_pred_z)
            print('z for clustering epoch:{:}, NMI:{:.4f}, ARI:{:.4f}, AMI:{:.4f}, FMI:{:.4f}, HOM:{:.4f}, COM:{:.4f}, V:{:.4f}'.format(epoch, nmi_z, ari_z, ami_z, fmi_z, hom_z, com_z, v_z))

    model.z = q
    model.y_pred=y_pred_q
    return nmi

def encodeZ(model, X_RNA, X_ATAC):
    model.eval()
    _, _, _, _, z, q = model(X_RNA, X_ATAC)
    return z, q

