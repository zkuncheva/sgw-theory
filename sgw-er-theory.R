library(igraph)
library(RSpectra)
library(png)
devtools::load_all()
library(SGNprops)
tt = 30 # time layers
nd = 100 # node degree
n = 100 # number of nodes
p = 0.1 # probability of edge
corr=0.98 # layer correlation

raw_G = get_conn_ER_ML(p = p, n = n, tt = tt, corr = corr)
supra_adj_G = get_supra_adj(raw_G) 
nL = create_norm_laplacian_graph(supra_adj_G)


L1 = nL[1:n,1:n]
L12 = nL[1:n,(n+1):(2*n)]
Lambda = NULL
for (k in 0:(tt-1)){
  pp=2*k*pi/tt
  res=eigs(L1+2*cos(pp)*L12,k=dim(L1)[1], , opts = list(retvec = FALSE))
  Lambda = rbind(Lambda,sort(res$values,decr=F))
}

# save small explanatory plot
# png('cos_small_plot.png',height=180,width=300)
# par(mar=c(3.3,5.1,0.1,0.5))
# plot(0:(tt-1),cos(0:(tt-1)*2*pi/tt), pch=16, cex.axis=1.3, las=1, ylab="", xlab="", col='blue')
# title(xlab='k', ylab=expression(cos~"("~italic(k)~frac(italic(2)*italic(pi),italic("T"))~")"),cex.lab=1.3, line=2.25)
# dev.off()

# produce new big eigenvalues plot
#pic <- readPNG('cos_small_plot.png')
pdf('Eigen_CosPlot_fixed_notation_ER.pdf',heigh=7, width=10)
    cls = rainbow(ceiling(tt/2), start=0, end=5/6)
    par(family = "serif")
    for (k in 1:ceiling(tt/2)){
      if (k==1){
        plot(1:100,Lambda[k,], col=cls[k], ylim=c(0,max(Lambda)), las=1, xlim=c(0,120),cex.lab=1.3,cex.axis=1.3,cex=1.3, type='l', lty='dashed', xlab=expression(lambda ~ index), ylab=expression(lambda ~ value)
             , main=expression(100~Smallest~eigenvalues~of~italic(widetilde(L))~+~2~cos~"("~italic(k)~frac(italic(2)*italic(pi),italic("T"))~")"~italic(widetilde(L)[W])~"for"~italic(k)==italic("0,1,...,T-1")))
      }else{
        points(1:100,Lambda[k,], col=cls[k], type='l', lty='dashed')
      }
    }
    legend.text=paste0("k=",(ceiling(tt/2)-1):0,',',ceiling(tt/2):(tt-1))
    legend('topright', lty='dashed',legend=legend.text,col=rev(cls), bty="n",cex=1.3, lwd=2)
    
    par(fig = c(grconvertX(c(35,100), from="user", to="ndc"),
                grconvertY(c(-0.05,0.75), from="user", to="ndc")),
        mar=c(3.3,5.1,0.1,0.5),
        new=TRUE)
    plot(0:(tt-1),cos(0:(tt-1)*2*pi/tt), pch=16, cex.axis=1.3, las=1, ylab="", xlab="", col='blue')
    title(xlab='k', ylab=expression(cos~"("~italic(k)~frac(italic(2)*italic(pi),italic("T"))~")"),cex.lab=1.3, line=2.25)
    #rasterImage(pic,50,0,100,0.7)
dev.off()

# get sin and cos plots
pdf('Cos-Sin-plots.pdf',height=14,width=14)
par(family = "serif", mar=c(5.1,4.5,4.1,0.2))
layout(matrix(c(1,1,2,3,3,1,1,4,3,3,5,5,6,7,7,5,5,8,7,7,9,9,10,11,11,9,9,12,11,11), nrow=6, ncol=5, byrow = TRUE))
for (k in 1:3){
  pp=2*k*pi/tt
  vec0 = eigs_sym(L1+2*cos(pp)*L12, k=1, which="SA")$vectors
  psi_cos=NULL
  psi_sin=NULL
  for (j in 1:tt){
    psi_cos = c(psi_cos,cos(2*pi*(j-1)*k/tt)*vec0)
    psi_sin = c(psi_sin,sin(2*pi*(j-1)*k/tt)*vec0)
  }
    # estimated cos plot
  if (k==1){
    cos.title = expression(cos~"("~frac(italic(2)*italic(pi),italic("T"))*italic(j)*italic(hat(k))~")"*italic(v)~"for"~italic(hat(k))==italic("1")~and~italic(j)==italic("0,1,...,T-1"))
    sin.title= expression(sin~"("~frac(italic(2)*italic(pi),italic("T"))*italic(j)*italic(hat(k))~")"*italic(v)~"for"~italic(hat(k))==italic("1")~and~italic(j)==italic("0,1,...,T-1"))
    small.cos.title = expression(cos~"("~frac(italic(2)*italic(pi),italic("T"))~italic(j)~") for"~italic(k)==italic("1"))
    small.sin.title = expression(sin~"("~frac(italic(2)*italic(pi),italic("T"))~italic(j)~") for"~italic(k)==italic("1"))
  }
  if (k==2){
    cos.title = expression(cos~"("~frac(italic(2)*italic(pi),italic("T"))*italic(j)*italic(hat(k))~")"*italic(v)~"for"~italic(hat(k))==italic("2")~and~italic(j)==italic("0,1,...,T-1"))
    sin.title= expression(sin~"("~frac(italic(2)*italic(pi),italic("T"))*italic(j)*italic(hat(k))~")"*italic(v)~"for"~italic(hat(k))==italic("2")~and~italic(j)==italic("0,1,...,T-1"))
    small.cos.title = expression(cos~"("~frac(italic(2)*italic(pi),italic("T"))~italic(j)~") for"~italic(k)==italic("2"))
    small.sin.title = expression(sin~"("~frac(italic(2)*italic(pi),italic("T"))~italic(j)~") for"~italic(k)==italic("2"))
    
  }
  if (k==3){
    cos.title = expression(cos~"("~frac(italic(2)*italic(pi),italic("T"))*italic(j)*italic(hat(k))~")"*italic(v)~"for"~italic(hat(k))==italic("3")~and~italic(j)==italic("0,1,...,T-1"))
    sin.title= expression(sin~"("~frac(italic(2)*italic(pi),italic("T"))*italic(j)*italic(hat(k))~")"*italic(v)~"for"~italic(hat(k))==italic("3")~and~italic(j)==italic("0,1,...,T-1"))
    small.cos.title = expression(cos~"("~frac(italic(2)*italic(pi),italic("T"))~italic(j)~") for"~italic(k)==italic("3"))
    small.sin.title = expression(sin~"("~frac(italic(2)*italic(pi),italic("T"))~italic(j)~") for"~italic(k)==italic("3"))
    
  }
  fig.txt=2
    plot(psi_cos, col='blue', type='l',xaxt='n', xlab='Time layers represented by nodes', ylab='Eigenvector estimation'
         , main=cos.title
         , cex.lab=fig.txt, cex.axis=fig.txt, cex=fig.txt, cex.main=fig.txt)
    axis(1,at=c((n/2)+n*((1:tt)-1)),label=1:tt, las=2, cex.axis=fig.txt)
    # theoretic cos plot
    if (k==1){
      plot(0:(tt-1),cos(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19,type='l',lty=1, xlab=NA, ylab=NA
           , main = small.cos.title
           , cex.lab=fig.txt, cex.axis=fig.txt, cex=fig.txt, cex.main=fig.txt)
      points(0:(tt-1),cos(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19)
      
    }else{
      plot(0:(tt-1),-cos(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19,type='l',lty=1, xlab=NA, ylab=NA
           , main = small.cos.title
           , cex.lab=fig.txt, cex.axis=fig.txt, cex=fig.txt, cex.main=fig.txt)
      points(0:(tt-1),-cos(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19)
    }
    # estimated sin plot
    plot(psi_sin, col='blue', type='l',xaxt='n', xlab='Time layers represented by nodes', ylab='Eigenvector estimation'
         , main=sin.title
         , cex.lab=fig.txt, cex.axis=fig.txt, cex=fig.txt, cex.main=fig.txt)
    axis(1,at=c((n/2)+n*((1:tt)-1)),label=1:tt, las=2, cex.axis=fig.txt)
    # theoretic sin plot
    if (k==1){
      plot(0:(tt-1),sin(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19,type='l',lty=1, xlab=NA, ylab=NA
           , main = small.sin.title
           , cex.lab=fig.txt, cex.axis=fig.txt, cex=fig.txt, cex.main=fig.txt)
      points(0:(tt-1),sin(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19)
    }else{
      plot(0:(tt-1),-sin(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19,type='l',lty=1, xlab=NA, ylab=NA
           , main = small.sin.title
           , cex.lab=fig.txt, cex.axis=fig.txt, cex=fig.txt, cex.main=fig.txt)
      points(0:(tt-1),-sin(2*pi*c(0:(tt-1)*k/tt)),col='blue',pch=19)
    }

}
dev.off()
