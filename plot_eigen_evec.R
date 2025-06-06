library(igraph)
library(RSpectra)
library(png)
devtools::load_all('~/SGNprops/')
#library(SGNprops)
results.dir='../sgw-results-test-Dec/'

# Simulate temporal network
tt = 30 # time layers
nd = 100 # node degree
n = 100 # number of nodes
p = 0.1 # probability of edge
corr=0.98  # layer correlation

raw_G = get_conn_ER_ML(p = p, n = n, tt = tt, corr = corr)
supra_adj_G = get_supra_adj(raw_G, omega = 0.01) 
nL = create_norm_laplacian_graph(supra_adj_G)

# get eigenvalues and eigenvectors for the 100 smallest
sRes = eigs_sym(nL, k=100, which='SA')
Lambda = sort(sRes$values,decr=F)
Vecs = sRes$vectors[,order(sRes$values,decreasing=FALSE)]
pdf(paste0(results.dir,'ER_eigen.pdf'), height=10, width=12)
par(mfrow=c(3,3))
cls = c('dodgerblue2','tomato1','goldenrod1','magenta4','palegreen3','turquoise2','deeppink4')
plot(Lambda, ylab='Eigenvalues', main='Eigenvalues', font.main = 2, pch=8, col='blue', cex.main = 1.5, cex.lab = 1.5, cex.axis=1.5 )
cls.full = rep(cls,each=100,times=30)
l.vec = 1:dim(Vecs)[1]
for (i in 1:6){
  plot(l.vec,Vecs[,i], col=cls.full, type='n', ylab='Eigenvector', xlab='Time layers\n represented by nodes', main=paste0('Eigenvector ',i)
       , font.main=2, xaxt='n', ylim=c(min(Vecs[,c(1:6)]),max(Vecs[,c(1:6)])), cex.main = 1.5, cex.lab = 1.5, cex.axis=1.5 )
  segments(l.vec[-length(l.vec)],Vecs[,i][-length(Vecs[,i])],l.vec[-1L],(Vecs[,i])[-1L],col=cls.full)
  axis(1, at=seq(from=50, to=2950, by=100), labels=c(1:30), las=2, cex=1.5)
  
}

plot(c(-1,0,1),c(0,0,0),xlim=c(-20,20),pch=16,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n')

plot(Vecs[,35],col=cls.full,type='n', ylab='Eigenvector', xlab='Time layers\n represented by nodes', main=paste0('Eigenvector ',35)
     , font.main=2, xaxt='n', cex.main = 1.5, cex.lab = 1.5, cex.axis=1.5 )
segments(l.vec[-length(l.vec)],Vecs[,35][-length(Vecs[,35])],l.vec[-1L],(Vecs[,35])[-1L],col=cls.full)
axis(1, at=seq(from=50, to=2950, by=100), labels=c(1:30), las=2, cex=1.5)
dev.off()