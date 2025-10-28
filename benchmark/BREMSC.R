library(BREMSC)
library(anndata)
library(rhdf5)

h5_file= H5Fopen("./GSE128639_BMNC.h5")
rna <- h5_file$X1
atac <- h5_file$X2
testRun <- BREMSC(atac, rna, K = 45, nChains = 2, nMCMC = 10)
write.table(testRun$clusterID, file = './GSE128639_bremscpred.csv', row.names=F, col.names = F)
write.table(testRun$posteriorProb, file = './GSE128639_bremscz.csv', row.names=F, col.names = F)