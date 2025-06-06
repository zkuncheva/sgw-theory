results.dir = "/Users/zhanakuncheva/sgw-results-test-Dec/results_dec23/"

results.list.mean = results.list.sd = list()
Omega = c(0.01,0.5,1,5)
Prob = c(0.03,0.04,0.05,0.08,0.1,0.3)
for (om.ix in 1:length(Omega)){
  omega = Omega[om.ix]
  o.mat.mean = o.mat.sd = NULL
    for (p.ix in 1:length(Prob)){
      p = Prob[p.ix]
      p.mat = NULL
      for (rep in 1:100){
        res = readRDS(file=paste0(results.dir,'error_p_',p,'_omega_',omega,'_rep_',rep,'.rds'))
        p.mat = rbind(p.mat,res)
      }
      o.mat.mean = rbind(o.mat.mean,apply(p.mat,2,mean))
      o.mat.sd = rbind(o.mat.sd,apply(p.mat,2,sd))
    }
  results.list.mean[[om.ix]] = o.mat.mean
  results.list.sd[[om.ix]] = o.mat.sd
}

pdf(paste0(results.dir,'errors-plots.pdf'), height=8, width=10)
mrkrs = c(4,19,6,3,15,20)
cls = c('green','orange','magenta','black','blue','red')
par(family = "serif", mfrow=c(2,2),mai=c(0.3,0.4,0.7,0.3))
for (i in 1:length(Omega)){
  omega = Omega[i]
 for (p in 1:length(Prob))
  if (p==1){
    plot(x=2:100, y=results.list.mean[[i]][p,-1], col=cls[p], type='l', lty=1, ylim=c(0,1.1)
         #, main=expression(Error~epsilon[i]~"of approximating supras-Laplacian eigenvectors")
         , xlab="i", ylab = expression(epsilon[i]), cex.lab=1.3, cex.axis=1.3, cex=1.3)
    points(x=2:100, y=results.list.mean[[i]][p,-1], pch=mrkrs[p], col=cls[p])
    polygon(c(2:100,100:2)
            ,c(results.list.mean[[i]][p,-1]+results.list.sd[[i]][p,-1],rev(results.list.mean[[i]][p,-1]-results.list.sd[[i]][p,-1]))
            , col=adjustcolor(cls[p], alpha.f=0.3), border=NA)
    abline(h=1, col='grey', lty=2)
    
  }else{
    points(x=2:100, y=results.list.mean[[i]][p,-1], pch=mrkrs[p], col=cls[p])
    points(x=2:100, y=results.list.mean[[i]][p,-1], pch=mrkrs[p], col=cls[p], type='l')
    polygon(c(2:100,100:2)
            ,c(results.list.mean[[i]][p,-1]+results.list.sd[[i]][p,-1],rev(results.list.mean[[i]][p,-1]-results.list.sd[[i]][p,-1]))
            , col=adjustcolor(cls[p], alpha.f=0.3), border=NA)
    abline(h=1, col='grey', lty=2)
    
  }

  if (omega==0.01){
    omega.subtitle = expression("Inter-layer weight"~omega*"=0.01")
  }
  if (omega==0.5){
    omega.subtitle = expression("Inter-layer weight"~omega*"=0.5")
  }
  if (omega==1){
    omega.subtitle = expression("Inter-layer weight"~omega*"=1")
  }
  if (omega==5){
    omega.subtitle = expression("Inter-layer weight"~omega*"=5")
  }
  title(omega.subtitle, line=1)
  legend('bottomright',paste0("p=",Prob),col=cls,pch=mrkrs,lty=1)
}
mtext(expression(bold(Error)~bold(epsilon[i])~bold("of approximating supra-Laplacian eigenvectors")), side = 3, line = - 2, outer = TRUE)

dev.off()
