library(igraph)
library(RSpectra)
library(png)
devtools::load_all()
library(SGNprops)


library(parallel)

tt = 30 # time layers
n = 100 # number of nodes
corr=0.9 # layer correlation

results.dir='../sgw-results/'
P = c(0.03,0.04,0.05,0.08,0.1,0.3) 
Omega=c(0.01,0.5,1,5) 
Rep=c(1:100)
params = list()
for (p in P){
  for (o in Omega){
    for (r in Rep){
      params[[length(params)+1]]=c(p,o,r)
    }
  }
}


sim.func=function(params){
  p = params[1]
  omega = params[2]
  rep=params[3]
  print(paste0('prob ',p))
    print(paste0('omega ',omega))
    ## for (rep in 1:100){
      print(paste0('rep ',rep))
      raw_G = get_conn_ER_ML(p = p, n = n, tt = tt, corr = corr)
      # get 0 eigenvector for each separate network
      Vec0 = NULL
      for (i in 1:tt){
        QA = create_norm_laplacian_graph(raw_G[[i]])
        res = eigs_sym(QA, k=1, which="SA")$vectors
        Vec0 = cbind(Vec0,res)
      }
      # get 0 eigenvectors for supra laplacian
      supra_adj_G = get_supra_adj(raw_G, omega=omega) 
      nL = create_norm_laplacian_graph(supra_adj_G)
      sRes = eigs_sym(nL, k=100, which='SA')
      sVec = sRes$vectors[,sort.list(sRes$values,decreasing=FALSE)]
      # create 0 padded matrix
      Mat0 = matrix(0,n*tt,tt)
      for (i in 1:tt){
        Mat0[((i-1)*n+1):(i*n),i] = Vec0[,i]
      }
      Err = NULL
      for (j in 1:100){
        beta = as.numeric(lm(sVec[,j]~0+Mat0)$coefficients[-1])
        error = sqrt(sum(t(Mat0 %*% beta-sVec[,j])^2))
        Err = c(Err,error)
      }
      saveRDS(Err,file=paste0(results.dir,'error_p_',p,'_omega_',omega,'_rep_',rep,'.rds'))
      
}

mclapply(params, sim.func, mc.cores=5)
