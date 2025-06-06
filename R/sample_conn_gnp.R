#' conn_gnp
#' 
#' A function to simulate Erdos-Renyi graph from another graph
#' @param G original graph
#' @param corr Pearson correlation between old and new graph
#' @param p probability of edge from old network 
#' @param max.it maximum number of iterations to try to get a connected network
#' @export
sample_conn_gnp = function(G, corr, p, max.it = 10){
  bool = FALSE
  iter = 1
  while (!bool){
    g=sample_correlated_gnp(G, corr, p)
    bool = is_connected(g)
    if (iter==max.it){
      stop('A connected graph could not be produced for max iterations.')
    }else if(!bool){
      iter = iter+1
    }
  }
  return(g)
}

