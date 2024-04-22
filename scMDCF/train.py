import torch
from torch.optim import Adam, Adadelta
import torch.nn.functional as F
from sklearn.cluster import KMeans
#from utils import eva, target_distribution

def pre_train(args, model, X_RNA, X_ATAC, y):
    optimizer = Adam(filter(lambda p: p.requires_grad, model.parameters()), lr = args.lr_pre, amsgrad=True)
    for epoch in range(args.epoch_pre):
        z_RNA, z_ATAC, rec_RNA, rec_ATAC, z, _ = model(X_RNA, X_ATAC)#
        
        loss_recrna = F.mse_loss(rec_RNA, X_RNA)
       
        loss_recatac = F.mse_loss(rec_ATAC, X_ATAC)
        
        cl_loss = model.crossview_contrastive_Loss(z_ATAC, z_RNA, lamb = args.lamb)
        loss = loss_recrna+loss_recatac+0.1*cl_loss
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch+=1
        if epoch%10==0:
            print('epoch:{}, loss:{:.4f}, loss_RNA:{:.4f}, loss_ATAC:{:.4f}, loss_cl:{:.4f}'.format(epoch, loss, loss_recrna, loss_recatac, cl_loss))

    torch.save(model.state_dict(), args.model_file)
    

def alt_train(args, model, X_RNA, X_ATAC, y):
    
    with torch.no_grad():
        _, _, _, _, z, _ = model(X_RNA, X_ATAC)
    kmeans = KMeans(n_clusters = args.n_clusters, n_init=20)
    y_pred = kmeans.fit_predict(z.data.cpu().numpy())
    
    model.cluster_layer.data = torch.tensor(kmeans.cluster_centers_).to(args.device)
    
    
    optimizer = Adadelta(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr_alt, rho=.8)
    
    for epoch in range(args.epoch_alt):
        model.train()
        z_RNA, z_ATAC, rec_RNA, rec_ATAC, z, q = model(X_RNA, X_ATAC)#
        p = model.target_distribution(q)
        loss_recrna = F.mse_loss(rec_RNA, X_RNA)
        loss_recatac = F.mse_loss(rec_ATAC, X_ATAC)
        
        cl_loss = model.crossview_contrastive_Loss(z_ATAC, z_RNA, lamb = args.lamb)
        
        loss_clu = model.cluster_loss(args, p, q)
        
        loss=0.1*loss_recrna+10*loss_recatac+1*loss_clu+5*cl_loss
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch+=1
        if epoch%10==0:
            print('epoch:{}, loss:{:.4f}, loss_RNA:{:.4f}, loss_ATAC:{:.4f}, loss_kl:{:.4f}, loss_cl:{:.4f}'.format(epoch, loss, loss_recrna, loss_recatac, loss_clu, cl_loss))
            with torch.no_grad():
                z, q = encodeZ(model, X_RNA, X_ATAC)
            
            kmeans = KMeans(n_clusters = args.n_clusters, n_init=20)
            y_pred_z = kmeans.fit_predict(z.data.cpu().numpy())          
   
    model.y_pred = y_pred_z
    model.latent = z

def encodeZ(model, X_RNA, X_ATAC):
    model.eval()
    _, _, _, _, z, q = model(X_RNA, X_ATAC)
    return z, q

