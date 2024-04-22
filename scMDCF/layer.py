import torch
import torch.nn as nn
from torch.nn import Parameter
import sys

class scMDCF_MLP(nn.Module):
    def __init__(self, layers):
        super(scMDCF_MLP, self).__init__()
        self.net=[]
        self.layer=None
        for i, (n_in, n_out) in enumerate(zip(layers[:-1], layers[1:])):
            self.net.append(nn.Linear(n_in, n_out))
            self.net.append(nn.BatchNorm1d(n_out))
            self.net.append(nn.ReLU())
        self.layer=nn.ModuleList(self.net)
        

    def forward(self, feature):
        latent = feature
        for layer in self.net:
            latent = layer(latent)
        
        return latent

class scMDCF(nn.Module):

    def __init__(self, args, v=1.0):
        super(scMDCF, self).__init__()
        self.encoder_RNA = scMDCF_MLP(args.layere_omics1_view)
        self.encoder_ATAC = scMDCF_MLP(args.layere_omics2_view)
        self.decoder_RNA = scMDCF_MLP(args.layerd_omics1_view)
        self.decoder_ATAC = scMDCF_MLP(args.layerd_omics2_view)
        self.fusion_MLP = scMDCF_MLP(args.fusion_layer)
        self.activate_layer = nn.Softmax()
        self.v = v
        self.alpha = args.alpha
        self.gamma = Parameter(torch.zeros(1))
        self.cluster_layer = nn.Parameter(torch.Tensor(args.n_clusters, args.zdim), requires_grad=True)
        torch.nn.init.xavier_normal_(self.cluster_layer.data)
        
    def forward(self, X_RNA, X_ATAC):
        z_RNA = self.encoder_RNA(X_RNA)
        z_ATAC = self.encoder_ATAC(X_ATAC)
        
        z = self.fusion_MLP(torch.cat([z_RNA, z_ATAC], dim=1))
        z_RNA = self.activate_layer(z_RNA)
        z_ATAC = self.activate_layer(z_ATAC)
        rec_ATAC = self.decoder_ATAC(z)
        rec_RNA = self.decoder_RNA(z)#

        q = 1.0 / (1.0 + torch.sum((z.unsqueeze(1) - self.cluster_layer)**2, dim=2) / self.alpha)
        q = q**((self.alpha+1.0)/2.0)
        q = (q.t() / torch.sum(q, dim=1)).t()

        return z_RNA, z_ATAC, rec_RNA, rec_ATAC, z, q
    
    def cal_latent(self, z):
        sum_y = torch.sum(torch.square(z), dim=1)
        num = -2.0 * torch.matmul(z, z.t()) + torch.reshape(sum_y, [-1, 1]) + sum_y
        num = num / self.alpha
        num = torch.pow(1.0 + num, -(self.alpha + 1.0) / 2.0)
        zerodiag_num = num - torch.diag(torch.diag(num))
        latent_p = (zerodiag_num.t() / torch.sum(zerodiag_num, dim=1)).t()
        return num, latent_p

    def crossview_contrastive_Loss(self, view1, view2, lamb=0.1, EPS=sys.float_info.epsilon):
        """Contrastive loss for maximizng the consistency"""
        _, k = view1.size()
        p_i_j = self.compute_joint(view1, view2)
        assert (p_i_j.size() == (k, k))

        p_i = p_i_j.sum(dim=1).view(k, 1).expand(k, k)
        p_j = p_i_j.sum(dim=0).view(1, k).expand(k, k)
        
        # Works with pytorch > 1.2
        p_i_j = torch.where(p_i_j < EPS, torch.tensor([EPS], device = p_i_j.device), p_i_j)
        p_j = torch.where(p_j < EPS, torch.tensor([EPS], device = p_j.device), p_j)
        p_i = torch.where(p_i < EPS, torch.tensor([EPS], device = p_i.device), p_i)

        loss =  -p_i_j * (torch.log(p_i_j) \
                        - (lamb + 1) * torch.log(p_j) \
                        - (lamb + 1) * torch.log(p_i))

        loss = loss.sum()

        return loss
    
    def cluster_loss(self, args, p, q):
        def kld(target, pred):
            return torch.mean(torch.sum(target*torch.log(target/(pred+1e-6)), dim=-1))
        kldloss = kld(p, q)
        return args.gamma*kldloss
    
    def compute_joint(self, view1, view2):
        """Compute the joint probability matrix P"""

        bn, k = view1.size()
        assert (view2.size(0) == bn and view2.size(1) == k)

        p_i_j = view1.unsqueeze(2) * view2.unsqueeze(1)
        p_i_j = p_i_j.sum(dim=0)
        p_i_j = (p_i_j + p_i_j.t()) / 2.  # symmetrise
        p_i_j = p_i_j / p_i_j.sum()  # normalise

        return p_i_j    
    
    def target_distribution(self, Q):
        """
        calculate the target distribution (student-t distribution)
        Args:
            Q: the soft assignment distribution
        Returns: target distribution P
        """
        weight = Q ** 2 / Q.sum(0)
        P = (weight.t() / weight.sum(1)).t()
        return P