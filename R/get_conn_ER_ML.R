#' get_ER_ML
#' 
#' A function to simulate multilayered slowly changing Erdos-Renyi network
#' @param p edge probability
#' @param n numebr of nodes per layer
#' @param tt number of time layers
#' @param corr pearson correlation between layers
#' @export
get_conn_ER_ML = function(p, n, tt, corr = 0.98){
  if (p<=0 || p>1){stop('Supplied probability p should be between 0 and 1')}
  if (n<=1){stop('Supplied number of nodes n should be greater than 1.')}
  if (n*p<3){stop('The model has very low chance to create a connected graph. Increase n, p or disable conn.')}
  if (tt<1){stop('Supplied number of time layers is less than 1. Increase to tt>=1.')}
  G = list()
  Flag = 0
  while (Flag==0){
    dat = try(conn_gnp(n,p,max.it=100) ,silent = TRUE)
    if (class(dat)=='try-error'){
      Flag = 0
    }else{
      Flag = 1
      G[[1]] = dat
    }
  }
  if (tt>2){
    for (time in 2:tt){
      Flag = 0
      while (Flag==0){
        dat = try(sample_conn_gnp(G[[time-1]], corr=corr, p = p) ,silent = TRUE)
        if (class(dat)=='try-error'){
          Flag = 0
        }else{
          Flag = 1
          G[[time]] = dat
        }
      }
    }
  }
  # create self loops, e.g. diag(G)=1
  #G = lapply(G, FUN=function(g){g[from=V(g),to=V(g)] = 1; g})     
  G = lapply(G, FUN=function(g){as_adjacency_matrix(g, sparse=FALSE)})
  return(G)
}

